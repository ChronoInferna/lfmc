#pragma once

#include "lfmc/strategy/strategy_metrics.hpp"
#include "lfmc/strategy/strategy_runner.hpp"

#include <algorithm>
#include <memory>
#include <thread>
#include <vector>

namespace lfmc {

template <StochasticProcess SP, NumericalScheme<SP> NS> class MultiStrategyRunner {
  private:
    std::vector<std::unique_ptr<StrategyRunner<SP, NS>>> strategies;
    std::vector<std::thread> threads;

  public:
    void add_strategy(std::unique_ptr<StrategyRunner<SP, NS>> strategy) {
        strategies.push_back(std::move(strategy));
    }

    // Run all strategies for a warmup period
    std::expected<void, std::string> run_warmup(size_t steps, double T, size_t warmup_iterations) {
        threads.clear();
        threads.reserve(strategies.size());

        for (auto& strategy : strategies) {
            threads.emplace_back([&strategy, steps, T, warmup_iterations]() {
                for (size_t i = 0; i < warmup_iterations; ++i) {
                    auto result = strategy->run_batch(steps, T);
                    if (!result) {
                        // Error occurred - need this to get log
                        return;
                    }
                }
            });
        }

        for (auto& thread : threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }

        return {};
    }

    // Get current metricss
    std::vector<StrategyMetrics> get_all_metrics() const {
        std::vector<StrategyMetrics> metrics;
        metrics.reserve(strategies.size());

        for (const auto& strategy : strategies) {
            metrics.push_back(strategy->get_metrics());
        }

        return metrics;
    }

    // Find best strat
    size_t best_strategy_index() const {
        auto metrics = get_all_metrics();
        if (metrics.empty()) {
            return 0;
        }
        return std::min_element(metrics.begin(), metrics.end()) - metrics.begin();
    }

    size_t strategy_count() const noexcept {
        return strategies.size();
    }
};

} // namespace lfmc