#include "lfmc/adaptive/adaptive_estimator.hpp"

#include <algorithm>
#include <chrono>
#include <numeric>

namespace lfmc {

StrategyStats AdaptiveVarianceReduction::run_strategy_batch(const std::string& name,
                                                            const SamplerFn& sampler,
                                                            size_t n_samples, uint64_t seed) {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto samples = sampler(n_samples, seed);
    auto t1 = std::chrono::high_resolution_clock::now();

    const double wall_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    auto [mean, variance] = detail::mean_variance(samples);

    return StrategyStats{
        .name = name,
        .mean = mean,
        .sample_variance = variance,
        .wall_time_ms = wall_ms,
        .n_samples = samples.size(),
    };
}

// Precision weighting: w_k* = (1 / sigma_k^2) / sum_j (1 / sigma_j^2)
std::vector<double>
AdaptiveVarianceReduction::compute_precision_weights(const std::vector<StrategyStats>& stats) {
    std::vector<double> precisions(stats.size());
    double total = 0.0;
    for (size_t k = 0; k < stats.size(); ++k) {
        precisions[k] = (stats[k].sample_variance > 0.0) ? 1.0 / stats[k].sample_variance : 0.0;
        total += precisions[k];
    }

    std::vector<double> weights(stats.size());
    if (total > 0.0) {
        for (size_t k = 0; k < stats.size(); ++k)
            weights[k] = precisions[k] / total;
    } else {
        // Degenerate: all variances are zero - use equal weights
        std::fill(weights.begin(), weights.end(), 1.0 / static_cast<double>(stats.size()));
    }
    return weights;
}

ASVRResult AdaptiveVarianceReduction::run(std::vector<NamedStrategy> strategies,
                                          size_t total_samples, ASVRConfig config) {
    const size_t K = strategies.size();
    assert(K >= 1);

    const size_t hw = std::max(size_t{1}, config.n_threads);

    // Budget split
    const size_t n_explore_each =
        std::max(config.min_exploration_per_strategy,
                 static_cast<size_t>(config.exploration_fraction *
                                     static_cast<double>(total_samples) / static_cast<double>(K)));
    const size_t n_explore_total = K * n_explore_each;
    const size_t n_exploit =
        (n_explore_total < total_samples) ? total_samples - n_explore_total : 0;

    // Phase 1: Parallel exploration — each strategy with seed make_seed(k, 0)
    const size_t explore_thread_cap = std::max(
        size_t{1},
        std::min(K, static_cast<size_t>(config.exploration_fraction * static_cast<double>(hw))));

    std::vector<StrategyStats> exploration_stats(K);

    // Launch exploration in waves of explore_thread_cap threads
    for (size_t wave_start = 0; wave_start < K; wave_start += explore_thread_cap) {
        const size_t wave_end = std::min(wave_start + explore_thread_cap, K);
        std::vector<std::jthread> threads;
        threads.reserve(wave_end - wave_start);
        for (size_t k = wave_start; k < wave_end; ++k) {
            threads.emplace_back([&, k]() {
                exploration_stats[k] = run_strategy_batch(strategies[k].name, strategies[k].sampler,
                                                          n_explore_each, detail::make_seed(k, 0));
            });
        }
        // jthreads join on destruction at end of wave block
    }

    // -- Compute precision weights ------------------------------------------
    auto weights = compute_precision_weights(exploration_stats);

    // -- Allocate exploitation budget ----------------------------------------
    // n_k = floor(w_k* * n_exploit); remainder given to the highest-weight strategy
    std::vector<size_t> exploit_counts(K, 0);
    if (n_exploit > 0) {
        size_t assigned = 0;
        for (size_t k = 0; k < K; ++k) {
            exploit_counts[k] = static_cast<size_t>(weights[k] * static_cast<double>(n_exploit));
            assigned += exploit_counts[k];
        }
        if (assigned < n_exploit) {
            const size_t best_k = static_cast<size_t>(
                std::max_element(weights.begin(), weights.end()) - weights.begin());
            exploit_counts[best_k] += (n_exploit - assigned);
        }
    }

    // Phase 2: Parallel exploitation — seed make_seed(k, 1) independent of phase 0
    const size_t exploit_thread_cap = std::max(size_t{1}, hw);

    std::vector<StrategyStats> exploit_stats(K);
    {
        std::vector<size_t> active;
        for (size_t k = 0; k < K; ++k)
            if (exploit_counts[k] > 0)
                active.push_back(k);

        for (size_t wave_start = 0; wave_start < active.size(); wave_start += exploit_thread_cap) {
            const size_t wave_end = std::min(wave_start + exploit_thread_cap, active.size());
            std::vector<std::jthread> threads;
            threads.reserve(wave_end - wave_start);
            for (size_t i = wave_start; i < wave_end; ++i) {
                const size_t k = active[i];
                threads.emplace_back([&, k]() {
                    exploit_stats[k] =
                        run_strategy_batch(strategies[k].name, strategies[k].sampler,
                                           exploit_counts[k], detail::make_seed(k, 1));
                });
            }
        }
    }

    // Final estimate: mu_ASVR = sum_k w_k* * mu_k^exploit (exploitation samples only)
    // Unbiasedness guaranteed because seeds for phase 1 and 2 are independent
    double mu_exploit = 0.0;
    double w_sum = 0.0;
    for (size_t k = 0; k < K; ++k) {
        if (exploit_counts[k] > 0) {
            mu_exploit += weights[k] * exploit_stats[k].mean;
            w_sum += weights[k];
        }
    }
    if (w_sum > 0.0 && w_sum < 1.0 - 1e-12)
        mu_exploit /= w_sum; // renormalise if some strategies got zero allocation

    const double estimate = (n_exploit > 0) ? mu_exploit : exploration_stats[0].mean;

    // Variance estimation: Var(mu_ASVR) = sum_k w_k*^2 * sigma_k^2 / n_k
    // Note: exploration variance slightly over-estimates CV strategy exploitation variance
    double var_exploit = 0.0;
    for (size_t k = 0; k < K; ++k) {
        if (exploit_counts[k] > 0) {
            var_exploit += weights[k] * weights[k] * exploration_stats[k].sample_variance /
                           static_cast<double>(exploit_counts[k]);
        }
    }
    // Renormalise in the same proportion as the estimate was renormalised
    if (w_sum > 0.0 && w_sum < 1.0 - 1e-12)
        var_exploit /= (w_sum * w_sum);

    const double total_variance = (n_exploit > 0) ? var_exploit
                                                  : exploration_stats[0].sample_variance /
                                                        static_cast<double>(n_explore_each);

    // Plain MC baseline for comparison: sigma_0^2 / N using the first strategy
    // (conventionally plain MC; the caller must ensure strategies[0] is plain MC)
    const double plain_mc_var =
        exploration_stats[0].sample_variance / static_cast<double>(total_samples);
    const double vr_ratio = (total_variance > 0.0) ? plain_mc_var / total_variance : 1.0;

    return ASVRResult{
        .estimate = estimate,
        .estimated_variance = total_variance,
        .estimated_stderr = std::sqrt(total_variance),
        .exploration_stats = exploration_stats,
        .precision_weights = weights,
        .exploitation_counts = exploit_counts,
        .plain_mc_variance_estimate = plain_mc_var,
        .variance_reduction_ratio = vr_ratio,
        .n_exploration = n_explore_total,
        .n_exploitation = n_exploit,
        .total_samples = total_samples,
    };
}

} // namespace lfmc
