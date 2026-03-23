#pragma once

#include "lfmc/strategy/strategy_metrics.hpp"
#include "lfmc/strategy/strategy_worker.hpp"

#include <algorithm>
#include <memory>
#include <thread>
#include <vector>

// TODO I really don't enjoy the fact that everything is past the path generator is templated and
// therefore must be in all the header files. Is that necessary? Or can we do better...

namespace lfmc {

template <StochasticProcess SP, NumericalScheme<SP> NS> class Engine {
  private:
    std::vector<std::unique_ptr<StrategyWorker<SP, NS>>> strategies;
    std::vector<std::jthread> threads;

  public:
    void add_strategy(std::unique_ptr<StrategyWorker<SP, NS>> strategy) {
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
                        return; // Just exit the thread on error, we can check metrics later to see
                                // if it failed
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

    // TODO run main simulation loop - is warmup function above necessary or is that just what this
    // is?

    // Get current metrics
    std::vector<StrategyMetrics> get_all_metrics() const {
        std::vector<StrategyMetrics> metrics;
        metrics.reserve(strategies.size());

        for (const auto& strategy : strategies) {
            metrics.push_back(strategy->get_metrics());
        }

        return metrics;
    }

    // Find best strat
    std::expected<size_t, std::string> get_best_strategy_index() const {
        auto metrics = get_all_metrics();
        if (metrics.empty()) {
            return std::unexpected("No metrics available");
        }
        return std::min_element(metrics.begin(), metrics.end()) - metrics.begin();
    }

    size_t strategy_count() const noexcept {
        return strategies.size();
    }
};

} // namespace lfmc
