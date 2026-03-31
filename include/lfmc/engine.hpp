#pragma once

#include "lfmc/adaptive_estimator.hpp"
#include "lfmc/numerical_scheme.hpp"
#include "lfmc/payoff.hpp"
#include "lfmc/stochastic_process.hpp"
#include "lfmc/strategies.hpp"
#include "lfmc/types.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <expected>
#include <memory>
#include <string>
#include <thread>
#include <vector>

namespace lfmc {

// ─── IterativeEngine ──────────────────────────────────────────────────────────
//
// Each round:
//   1. All K strategies run on 1 "compete" thread each (they never drop out).
//   2. The current leader also gets n_exploit bonus threads.
//   3. All results are accumulated into running per-strategy statistics.
//   4. Precision weights (inverse-variance) are recomputed from ALL accumulated
//      data so far.
//   5. The strategy with the highest precision weight becomes the new leader.
//      If it changed, the bonus threads follow.
//   6. Repeat for n_rounds.
//
// Final estimate: precision-weighted mean of all strategies using their full
// accumulated statistics. Because each strategy's samples are i.i.d. and
// generated with round-unique seeds (independent across rounds and strategies),
// combining them this way is unbiased.
//
// The leader can change between rounds: a strategy that looked good early but
// had noisy variance estimates will lose the lead as more data accumulates. This
// means the engine genuinely re-evaluates every round rather than committing
// after a single pilot phase.

struct IterativeEngineConfig {
    // How many strategies compete simultaneously (1 compete thread each)
    size_t n_compete = 4;
    // Bonus threads assigned to the current leader each round
    size_t n_exploit = 8;
    // Number of competition rounds
    size_t n_rounds = 5;
    // Samples per thread per round
    size_t samples_per_thread = 2000;
};

struct RoundResult {
    size_t round;
    std::string leader;  // strategy leading after this round
    bool leader_changed; // did the leader switch this round?
    double estimate;     // precision-weighted estimate so far
    double estimated_stderr;
    std::vector<double> precision_weights; // one per strategy, after this round
};

struct IterativeEngineResult {
    double estimate;
    double estimated_stderr;
    std::string final_leader;
    std::vector<RoundResult> round_history; // diagnostics per round
    std::vector<StrategyStats> final_stats; // cumulative stats per strategy
    size_t total_samples;
};

template <StochasticProcess SP, NumericalScheme<SP> NS> class Engine {
    SP process_;
    NS scheme_;
    double control_mean_;
    size_t steps_;
    double T_;

  public:
    Engine(SP process, NS scheme, double control_mean, size_t steps, double T)
        : process_(process), scheme_(scheme), control_mean_(control_mean), steps_(steps), T_(T) {}

    std::expected<IterativeEngineResult, std::string> run(std::shared_ptr<Payoff> payoff,
                                                          IterativeEngineConfig config = {}) {
        auto all = build_all_strategies(process_, scheme_, payoff, control_mean_, steps_, T_);
        const size_t K = std::min(config.n_compete, all.size());
        if (K == 0)
            return std::unexpected("No strategies available");

        // ── Per-strategy running accumulators ─────────────────────────────────
        // We track (sum, sum_sq, n) across all rounds so precision weights can
        // be recomputed from the full history after each round.
        struct Accumulator {
            double sum = 0.0;
            double sum_sq = 0.0;
            size_t n = 0;

            void add(const std::vector<double>& samples) {
                for (double x : samples) {
                    sum += x;
                    sum_sq += x * x;
                    ++n;
                }
            }

            double mean() const {
                return n > 0 ? sum / static_cast<double>(n) : 0.0;
            }

            // Unbiased sample variance (n-1 denominator)
            double variance() const {
                if (n < 2)
                    return std::numeric_limits<double>::infinity();
                const double m = mean();
                return (sum_sq - static_cast<double>(n) * m * m) / static_cast<double>(n - 1);
            }
        };

        std::vector<Accumulator> accum(K);

        // ── Helper: compute precision weights from current accumulators ────────
        auto precision_weights = [&]() -> std::vector<double> {
            std::vector<double> prec(K);
            double total = 0.0;
            for (size_t k = 0; k < K; ++k) {
                double var = accum[k].variance();
                prec[k] = (var > 0.0 && std::isfinite(var)) ? 1.0 / var : 0.0;
                total += prec[k];
            }
            std::vector<double> w(K);
            if (total > 0.0) {
                for (size_t k = 0; k < K; ++k)
                    w[k] = prec[k] / total;
            } else {
                std::fill(w.begin(), w.end(), 1.0 / static_cast<double>(K));
            }
            return w;
        };

        size_t leader = 0; // index of current leader
        bool first_round = true;
        std::vector<RoundResult> history;
        history.reserve(config.n_rounds);

        size_t total_samples = 0;

        for (size_t round = 0; round < config.n_rounds; ++round) {
            // ── Determine thread allocation for this round ─────────────────────
            // Every strategy gets 1 compete thread.
            // The current leader also gets n_exploit bonus threads.
            // On the first round there is no established leader yet, so the
            // bonus threads are spread evenly (floor(n_exploit/K) each, with
            // remainder going to strategy 0 as a tiebreak).

            std::vector<size_t> thread_counts(K, 1); // 1 compete thread each
            if (first_round) {
                const size_t bonus_each = config.n_exploit / K;
                const size_t leftover = config.n_exploit % K;
                for (size_t k = 0; k < K; ++k)
                    thread_counts[k] += bonus_each;
                thread_counts[0] += leftover;
            } else {
                thread_counts[leader] += config.n_exploit;
            }

            // ── Run all threads for this round ─────────────────────────────────
            // Each thread generates samples_per_thread samples.
            // Seed = make_seed(strategy_idx * 10000 + thread_idx, round) ensures
            // independence across strategies, threads within a strategy, and rounds.

            // Collect samples per strategy (multiple threads concatenated)
            std::vector<std::vector<double>> round_samples(K);
            {
                // Flatten to (strategy_idx, thread_idx) pairs for easy launch
                struct Job {
                    size_t k;
                    size_t t;
                };
                std::vector<Job> jobs;
                for (size_t k = 0; k < K; ++k)
                    for (size_t t = 0; t < thread_counts[k]; ++t)
                        jobs.push_back({k, t});

                std::vector<std::vector<double>> thread_results(jobs.size());
                {
                    std::vector<std::jthread> threads;
                    threads.reserve(jobs.size());
                    for (size_t j = 0; j < jobs.size(); ++j) {
                        threads.emplace_back([&, j]() {
                            const size_t k = jobs[j].k;
                            const size_t t = jobs[j].t;
                            const uint64_t seed = detail::make_seed(k * 10000 + t, round);
                            thread_results[j] = all[k].second(config.samples_per_thread, seed);
                        });
                    }
                } // jthreads join here

                // Merge thread results into per-strategy buckets
                for (size_t j = 0; j < jobs.size(); ++j) {
                    auto& dest = round_samples[jobs[j].k];
                    auto& src = thread_results[j];
                    dest.insert(dest.end(), src.begin(), src.end());
                }
            }

            // ── Update accumulators ────────────────────────────────────────────
            for (size_t k = 0; k < K; ++k) {
                accum[k].add(round_samples[k]);
                total_samples += round_samples[k].size();
            }

            // ── Recompute precision weights and current leader ─────────────────
            auto weights = precision_weights();
            size_t new_leader = static_cast<size_t>(
                std::max_element(weights.begin(), weights.end()) - weights.begin());

            // ── Build round diagnostic ─────────────────────────────────────────
            double est = 0.0;
            for (size_t k = 0; k < K; ++k)
                est += weights[k] * accum[k].mean();

            // Var(precision-weighted sum) = sum_k w_k^2 * sigma_k^2 / n_k
            double var_est = 0.0;
            for (size_t k = 0; k < K; ++k) {
                double var_k = accum[k].variance();
                if (std::isfinite(var_k) && accum[k].n > 0)
                    var_est += weights[k] * weights[k] * var_k / static_cast<double>(accum[k].n);
            }

            history.push_back({
                .round = round + 1,
                .leader = all[new_leader].first,
                .leader_changed = (!first_round && new_leader != leader),
                .estimate = est,
                .estimated_stderr = std::sqrt(var_est),
                .precision_weights = weights,
            });

            leader = new_leader;
            first_round = false;
        }

        // ── Final estimate ────────────────────────────────────────────────────
        auto final_weights = precision_weights();
        double final_est = 0.0;
        for (size_t k = 0; k < K; ++k)
            final_est += final_weights[k] * accum[k].mean();

        double final_var = 0.0;
        for (size_t k = 0; k < K; ++k) {
            double var_k = accum[k].variance();
            if (std::isfinite(var_k) && accum[k].n > 0)
                final_var +=
                    final_weights[k] * final_weights[k] * var_k / static_cast<double>(accum[k].n);
        }

        // ── Pack cumulative StrategyStats for diagnostics ──────────────────────
        std::vector<StrategyStats> final_stats(K);
        for (size_t k = 0; k < K; ++k) {
            final_stats[k] = {
                .name = all[k].first,
                .mean = accum[k].mean(),
                .sample_variance = accum[k].variance(),
                .wall_time_ms = 0.0, // not tracked per-strategy across rounds
                .n_samples = accum[k].n,
            };
        }

        return IterativeEngineResult{
            .estimate = final_est,
            .estimated_stderr = std::sqrt(final_var),
            .final_leader = all[leader].first,
            .round_history = history,
            .final_stats = final_stats,
            .total_samples = total_samples,
        };
    }
};

} // namespace lfmc
