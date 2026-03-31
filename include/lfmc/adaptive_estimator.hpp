#pragma once

#include "lfmc/numerical_scheme.hpp"
#include "lfmc/payoff.hpp"
#include "lfmc/stochastic_process.hpp"
#include "lfmc/types.hpp"

#include <cassert>
#include <cstdint>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace lfmc {

// A strategy sampler: given n samples and a seed, returns n iid samples of the
// (possibly variance-reduced) payoff estimator. The seed guarantees RNG independence
// between strategies and between exploration/exploitation phases.
using SamplerFn = std::function<std::vector<double>(size_t n_samples, uint64_t seed)>;

// Per-strategy statistics collected during the exploration phase
struct StrategyStats {
    std::string name;
    double mean;
    double sample_variance; // unbiased (n-1 denominator)
    double wall_time_ms;
    size_t n_samples;
};

// Full result of one ASVR run
struct ASVRResult {
    double estimate;
    double estimated_variance; // Var(mu_ASVR), estimated from exploration data
    double estimated_stderr;   // sqrt(estimated_variance)

    std::vector<StrategyStats> exploration_stats; // one entry per strategy
    std::vector<double> precision_weights;        // w_k* used to allocate exploitation
    std::vector<size_t> exploitation_counts;      // actual n_k samples per strategy

    // Diagnostics — useful for paper tables
    double plain_mc_variance_estimate; // sigma_0^2 / N (first strategy treated as plain MC)
    double variance_reduction_ratio;   // plain_mc_variance / estimated_variance
    size_t n_exploration;
    size_t n_exploitation;
    size_t total_samples;
};

struct ASVRConfig {
    // Fraction of total_samples used for exploration (split equally across K strategies)
    double exploration_fraction = 0.1;
    // Total thread budget; exploration uses up to K threads, exploitation uses the rest
    size_t n_threads = std::thread::hardware_concurrency();
    // Minimum exploration samples per strategy regardless of exploration_fraction.
    // Needs to be large enough for a stable variance estimate (~100+ samples).
    size_t min_exploration_per_strategy = 100;
};

// ─── ASVR Algorithm ──────────────────────────────────────────────────────────
//
// Phase 1 (Exploration): K strategies run in parallel, each generating
//   n_explore_each = max(min_per_strategy, alpha*N/K) samples.
//   Each strategy uses seed make_seed(k, 0) for full RNG independence.
//
// Weight computation: w_k* = (1/sigma_k^2) / sum_j(1/sigma_j^2)
//   (inverse-variance / "precision" weighting)
//
// Phase 2 (Exploitation): strategy k receives n_k = floor(w_k* * n_exploit) samples,
//   each run using seed make_seed(k, 1) — independent of the exploration phase.
//   Strategies with n_k > 0 run in parallel.
//
// Combination (unbiased by independence):
//   mu_explore = (1/K) * sum_k mu_k^explore      (equal-weight pooled)
//   mu_exploit = sum_k w_k* * mu_k^exploit        (precision-weighted)
//   mu_ASVR    = alpha * mu_explore + (1-alpha) * mu_exploit
//
// The key unbiasedness guarantee: w_k* depends only on exploration samples;
// exploitation samples are generated with fresh RNG seeds and are therefore
// independent of the weight selection event.

class AdaptiveVarianceReduction {
  public:
    static ASVRResult run(std::vector<std::pair<std::string, SamplerFn>> strategies,
                          size_t total_samples, ASVRConfig config = {});

    static StrategyStats run_strategy_batch(const std::string& name, const SamplerFn& sampler,
                                            size_t n_samples, uint64_t seed);

  private:
    static std::vector<double> compute_precision_weights(const std::vector<StrategyStats>& stats);
};

// ─── Internal utilities ───────────────────────────────────────────────────────

namespace detail {

// Splitmix64-based seed derivation: strategy k, phase p → unique uint64_t seed.
// Two calls with the same (k, p) give the same seed; different (k, p) pairs are
// statistically independent because they hit different streams of splitmix64.
inline uint64_t make_seed(size_t k, size_t phase) noexcept {
    uint64_t x = (static_cast<uint64_t>(k) * 0x9e3779b97f4a7c15ULL)
               ^ (static_cast<uint64_t>(phase) * 0x6c62272e07bb0142ULL)
               ^ 0xdeadbeefcafe0000ULL;
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

// Generate `n` draws from N(0,1)
inline std::vector<double> gen_normals(size_t n, std::mt19937_64& rng) {
    std::normal_distribution<double> dist{0.0, 1.0};
    std::vector<double> v(n);
    for (auto& x : v)
        x = dist(rng);
    return v;
}

// Correct single-path generator (avoids the pre-allocate + push_back double-length bug
// present in the existing PathGenerator)
template <StochasticProcess SP, NumericalScheme<SP> NS>
Path generate_path(const SP& process, const NS& scheme, const Normals& normals, size_t steps,
                   double T) {
    const double dt = T / static_cast<double>(steps);
    Path path;
    path.reserve(steps + 1);
    double x = process.initial();
    path.push_back(x);
    double t = 0.0;
    for (size_t i = 0; i < steps; ++i) {
        x = scheme.step(process, t, x, dt, normals[i]);
        path.push_back(x);
        t += dt;
    }
    return path;
}

// Unbiased sample mean and variance (Welford would be more numerically stable for large n,
// but two-pass is fine for the batch sizes used here)
inline std::pair<double, double> mean_variance(const std::vector<double>& samples) {
    assert(!samples.empty());
    const double n = static_cast<double>(samples.size());
    double sum = 0.0;
    for (double x : samples)
        sum += x;
    const double mean = sum / n;
    double sq = 0.0;
    for (double x : samples) {
        double d = x - mean;
        sq += d * d;
    }
    // n-1 denominator for unbiased variance; guard against n=1
    const double variance = (samples.size() > 1) ? sq / (n - 1.0) : 0.0;
    return {mean, variance};
}

} // namespace detail

// ─── Strategy factory functions ──────────────────────────────────────────────
//
// Each factory captures the pricing parameters and returns a SamplerFn.
// The SamplerFn is called with (n_samples, seed) and returns exactly n_samples
// iid draws from the (possibly variance-reduced) estimator of E[payoff].

// Plain pseudo-random Monte Carlo — baseline strategy
template <StochasticProcess SP, NumericalScheme<SP> NS>
SamplerFn make_plain_mc_sampler(SP process, NS scheme, std::shared_ptr<Payoff> payoff,
                                 size_t steps, double T) {
    return [process, scheme, payoff, steps, T](size_t n, uint64_t seed) -> std::vector<double> {
        std::mt19937_64 rng{seed};
        std::vector<double> samples;
        samples.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            auto normals = detail::gen_normals(steps, rng);
            auto path = detail::generate_path(process, scheme, normals, steps, T);
            auto result = payoff->generate_payoffs({path});
            if (result && !result->empty() && !(*result)[0].empty())
                samples.push_back((*result)[0][0]);
        }
        return samples;
    };
}

// Antithetic variates — each sample is (payoff(Z) + payoff(-Z)) / 2
// Returns n samples, each consuming one pair of paths. Variance is reduced
// because the two paths are negatively correlated for monotone payoffs.
template <StochasticProcess SP, NumericalScheme<SP> NS>
SamplerFn make_antithetic_sampler(SP process, NS scheme, std::shared_ptr<Payoff> payoff,
                                   size_t steps, double T) {
    return [process, scheme, payoff, steps, T](size_t n, uint64_t seed) -> std::vector<double> {
        std::mt19937_64 rng{seed};
        std::vector<double> samples;
        samples.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            auto z = detail::gen_normals(steps, rng);
            Normals neg_z(steps);
            for (size_t j = 0; j < steps; ++j)
                neg_z[j] = -z[j];

            auto path_pos = detail::generate_path(process, scheme, z, steps, T);
            auto path_neg = detail::generate_path(process, scheme, neg_z, steps, T);

            auto r_pos = payoff->generate_payoffs({path_pos});
            auto r_neg = payoff->generate_payoffs({path_neg});

            if (r_pos && r_neg && !r_pos->empty() && !r_neg->empty() && !(*r_pos)[0].empty() &&
                !(*r_neg)[0].empty())
                samples.push_back(0.5 * ((*r_pos)[0][0] + (*r_neg)[0][0]));
        }
        return samples;
    };
}

// Control variates with in-batch OLS beta estimation.
// target_payoff: the option price we're estimating
// control_payoff: a payoff with known expectation `control_mean`
//   (e.g. EuropeanCall(0.0) gives S_T; E[S_T] = S0 * exp(mu*T))
// Returns: X_i - beta_hat * (Y_i - E[Y]) for each i
template <StochasticProcess SP, NumericalScheme<SP> NS>
SamplerFn make_control_variate_sampler(SP process, NS scheme, std::shared_ptr<Payoff> target,
                                        std::shared_ptr<Payoff> control, double control_mean,
                                        size_t steps, double T) {
    return [process, scheme, target, control, control_mean, steps,
            T](size_t n, uint64_t seed) -> std::vector<double> {
        std::mt19937_64 rng{seed};
        std::vector<double> xs, ys;
        xs.reserve(n);
        ys.reserve(n);

        for (size_t i = 0; i < n; ++i) {
            auto z = detail::gen_normals(steps, rng);
            auto path = detail::generate_path(process, scheme, z, steps, T);
            auto rx = target->generate_payoffs({path});
            auto ry = control->generate_payoffs({path});
            if (rx && ry && !rx->empty() && !ry->empty() && !(*rx)[0].empty() &&
                !(*ry)[0].empty()) {
                xs.push_back((*rx)[0][0]);
                ys.push_back((*ry)[0][0]);
            }
        }

        // OLS beta: minimises Var(X - beta*(Y - E[Y]))
        const size_t m = xs.size();
        double mean_x = 0.0, mean_y = 0.0;
        for (size_t i = 0; i < m; ++i) {
            mean_x += xs[i];
            mean_y += ys[i];
        }
        mean_x /= static_cast<double>(m);
        mean_y /= static_cast<double>(m);

        double cov_xy = 0.0, var_y = 0.0;
        for (size_t i = 0; i < m; ++i) {
            cov_xy += (xs[i] - mean_x) * (ys[i] - mean_y);
            var_y += (ys[i] - mean_y) * (ys[i] - mean_y);
        }
        const double beta = (var_y > 0.0) ? cov_xy / var_y : 0.0;

        std::vector<double> samples(m);
        for (size_t i = 0; i < m; ++i)
            samples[i] = xs[i] - beta * (ys[i] - control_mean);
        return samples;
    };
}

// Antithetic + control variates combined.
// Each raw sample is the antithetic average; OLS beta is estimated from those averages.
template <StochasticProcess SP, NumericalScheme<SP> NS>
SamplerFn make_antithetic_cv_sampler(SP process, NS scheme, std::shared_ptr<Payoff> target,
                                      std::shared_ptr<Payoff> control, double control_mean,
                                      size_t steps, double T) {
    return [process, scheme, target, control, control_mean, steps,
            T](size_t n, uint64_t seed) -> std::vector<double> {
        std::mt19937_64 rng{seed};
        std::vector<double> xs, ys;
        xs.reserve(n);
        ys.reserve(n);

        for (size_t i = 0; i < n; ++i) {
            auto z = detail::gen_normals(steps, rng);
            Normals neg_z(steps);
            for (size_t j = 0; j < steps; ++j)
                neg_z[j] = -z[j];

            auto path_p = detail::generate_path(process, scheme, z, steps, T);
            auto path_n = detail::generate_path(process, scheme, neg_z, steps, T);

            auto xp = target->generate_payoffs({path_p});
            auto xn = target->generate_payoffs({path_n});
            auto yp = control->generate_payoffs({path_p});
            auto yn = control->generate_payoffs({path_n});

            if (xp && xn && yp && yn && !xp->empty() && !xn->empty() && !yp->empty() &&
                !yn->empty() && !(*xp)[0].empty() && !(*xn)[0].empty() && !(*yp)[0].empty() &&
                !(*yn)[0].empty()) {
                xs.push_back(0.5 * ((*xp)[0][0] + (*xn)[0][0]));
                ys.push_back(0.5 * ((*yp)[0][0] + (*yn)[0][0]));
            }
        }

        const size_t m = xs.size();
        double mean_x = 0.0, mean_y = 0.0;
        for (size_t i = 0; i < m; ++i) {
            mean_x += xs[i];
            mean_y += ys[i];
        }
        mean_x /= static_cast<double>(m);
        mean_y /= static_cast<double>(m);

        double cov_xy = 0.0, var_y = 0.0;
        for (size_t i = 0; i < m; ++i) {
            cov_xy += (xs[i] - mean_x) * (ys[i] - mean_y);
            var_y += (ys[i] - mean_y) * (ys[i] - mean_y);
        }
        const double beta = (var_y > 0.0) ? cov_xy / var_y : 0.0;

        std::vector<double> samples(m);
        for (size_t i = 0; i < m; ++i)
            samples[i] = xs[i] - beta * (ys[i] - control_mean);
        return samples;
    };
}

} // namespace lfmc
