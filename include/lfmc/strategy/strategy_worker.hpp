#pragma once

#include "lfmc/estimator/control_variate_estimator.hpp"
#include "lfmc/estimator/monte_carlo_estimator.hpp"
#include "lfmc/pipeline/pipeline.hpp"
#include "lfmc/strategy/strategy_metrics.hpp"
#include "lfmc/strategy/types.hpp"
#include "lfmc/timing/timing.hpp"

#include <atomic>
#include <string>

namespace lfmc {

template <StochasticProcess SP, NumericalScheme<SP> NS> class StrategyWorker {
  private:
    Pipeline<SP, NS> pipeline;
    Strategy strategy;
    Timer timer;

    std::atomic<double> current_mean{0.0};
    std::atomic<double> current_variance{std::numeric_limits<double>::max()};
    std::atomic<size_t> samples_processed{0};

  public:
    StrategyWorker(Pipeline<SP, NS> pipeline, Strategy strategy)
        : pipeline(std::move(pipeline)), strategy(std::move(strategy)) {}

    // Run a single batch of simulations and update metrics I made in the other metrics file
    std::expected<void, std::string> run_batch(size_t steps, double T) {
        auto result = pipeline.run(steps, T);
        if (!result) {
            return std::unexpected(result.error());
        }

        auto* estimator = pipeline.get_estimator();

        // TODO possibly move statistics into base class for exposure so we don't need to dynamic
        // cast had to implement this, it wasn't in the base class, but need it to update the
        // metrics
        if (auto* mc_est = dynamic_cast<MonteCarloEstimator*>(estimator)) {
            current_mean.store(mc_est->mean(), std::memory_order_relaxed);
            current_variance.store(mc_est->variance(), std::memory_order_relaxed);
            samples_processed.store(mc_est->sample_count(), std::memory_order_relaxed);
        }
        // Try ControlVariateEstimator, same as above
        else if (auto* cv_est = dynamic_cast<ControlVariateEstimator*>(estimator)) {
            current_mean.store(cv_est->mean(), std::memory_order_relaxed);
            current_variance.store(cv_est->variance(), std::memory_order_relaxed);
            samples_processed.store(cv_est->sample_count(), std::memory_order_relaxed);
        }

        return {};
    }

    StrategyMetrics get_metrics() const noexcept {
        return StrategyMetrics{.strategy = strategy,
                               .mean = current_mean.load(std::memory_order_relaxed),
                               .variance = current_variance.load(std::memory_order_relaxed),
                               .std_error =
                                   std::sqrt(current_variance.load(std::memory_order_relaxed) /
                                             samples_processed.load(std::memory_order_relaxed)),
                               .samples = samples_processed.load(std::memory_order_relaxed),
                               .elapsed_ms = timer.elapsedMilliseconds()};
    }

    const Strategy& get_strategy() const noexcept {
        return strategy;
    }
};

} // namespace lfmc
