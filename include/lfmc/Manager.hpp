#pragma once

#include "NumericalScheme.hpp"
#include "Simulator.hpp"
#include "StochasticProcess.hpp"
#include "VarianceReductionStrategy.hpp"

// TODO
// #include <expected>
#include <memory>
// #include <thread>
// #include <vector>

/**
 * @file manager.hpp
 * @brief Manager class for variation reduction strategies using the Strategy Pattern.
 *
 * This file defines the `Manager` class, which is responsible for managing
 * different strategies and potentially coordinating between them during runtime.
 */

namespace lfmc {

template <StochasticProcess P, NumericalScheme<P> S, Payoff PO> class Manager {
  public:
    explicit Manager(P process, S scheme, PO payoff,
                     std::unique_ptr<VarianceReductionStrategy> strategy, double initialValue,
                     double timeToMaturity, size_t stepCount) noexcept
        : process_(std::move(process)), scheme_(std::move(scheme)), payoff_(std::move(payoff)),
          currentStrategy_(std::move(strategy)), initialValue_(initialValue),
          timeToMaturity_(timeToMaturity), stepCount_(stepCount) {}

    // Temporary for minimum viable product
    double simulate(size_t numThreads = 1000) {
        simulators_.clear();
        simulators_.reserve(numThreads);

        for (size_t i{}; i < numThreads; ++i) {
            simulators_.emplace_back(process_, scheme_, payoff_, initialValue_, timeToMaturity_,
                                     stepCount_);
        }

        // TODO This should eventually be a thread pool
        double totalResult = 0.0;
        for (auto& simulator : simulators_) {
            totalResult += payoff_(simulator.generateTerminal());
        }

        return totalResult / static_cast<double>(numThreads);
    }

    // Temporary for minimum viable product
    std::pair<double, double> simulateWithError(size_t numThreads = 1000) {
        simulators_.clear();
        simulators_.reserve(numThreads);

        for (size_t i{}; i < numThreads; ++i) {
            simulators_.emplace_back(process_, scheme_, payoff_, initialValue_, timeToMaturity_,
                                     stepCount_);
        }

        // TODO This should eventually be a thread pool
        std::vector<double> results;
        results.reserve(numThreads);
        for (auto& simulator : simulators_) {
            results.push_back(payoff_(simulator.generateTerminal()));
        }

        double mean =
            std::accumulate(results.begin(), results.end(), 0.0) / static_cast<double>(numThreads);

        double variance = 0.0;
        for (const auto& r : results)
            variance += (r - mean) * (r - mean);
        variance /= static_cast<double>(numThreads - 1);
        double stdError = std::sqrt(variance / static_cast<double>(numThreads));

        return {mean, stdError};
    }

  private:
    P process_;
    S scheme_;
    PO payoff_;
    double initialValue_;
    double timeToMaturity_;
    size_t stepCount_;
    // TODO number of simulations, currently hardcoded to 1000

    std::unique_ptr<VarianceReductionStrategy> currentStrategy_;

    // std::vector<lfmc::Simulator<P, S, PO>> testingSimulators_;
    // std::vector<lfmc::Simulator<P, S, PO>> realSimulators_;

    std::vector<lfmc::Simulator<P, S, PO>> simulators_;
};

} // namespace lfmc
