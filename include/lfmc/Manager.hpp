#pragma once

#include "NumericalScheme.hpp"
#include "Simulator.hpp"
#include "StochasticProcess.hpp"
#include "VarianceReductionStrategy.hpp"

#include <memory>
#include <numeric>
#include <vector>
#include <cmath>

namespace lfmc {

template <StochasticProcess P, NumericalScheme<P> S, Payoff PO> 
class Manager {
public:
    explicit Manager(P process, S scheme, PO payoff,
                     std::unique_ptr<VarianceReductionStrategy> strategy, double initialValue,
                     double timeToMaturity, size_t stepCount) noexcept
        : process_(std::move(process)), scheme_(std::move(scheme)), payoff_(std::move(payoff)),
          currentStrategy_(std::move(strategy)), initialValue_(initialValue),
          timeToMaturity_(timeToMaturity), stepCount_(stepCount) {}

    double simulate(size_t numSimulations = 1000) {
        // Check if using Antithetic Variates
        if (dynamic_cast<AntitheticVariates*>(currentStrategy_.get())) {
            return simulateAntithetic(numSimulations);
        }
        
        // Standard Monte Carlo
        return simulateStandard(numSimulations);
    }

    std::pair<double, double> simulateWithError(size_t numSimulations = 1000) {
        // Check if using Antithetic Variates
        if (dynamic_cast<AntitheticVariates*>(currentStrategy_.get())) {
            return simulateAntitheticWithError(numSimulations);
        }
        
        // Standard Monte Carlo
        return simulateStandardWithError(numSimulations);
    }

private:
    // Standard Monte Carlo (no variance reduction)
    double simulateStandard(size_t numSimulations) {
        simulators_.clear();
        simulators_.reserve(numSimulations);

        for (size_t i{}; i < numSimulations; ++i) {
            simulators_.emplace_back(process_, scheme_, payoff_, initialValue_, 
                                     timeToMaturity_, stepCount_);
        }

        double totalResult = 0.0;
        for (auto& simulator : simulators_) {
            totalResult += payoff_(simulator.generateTerminal());
        }

        return totalResult / static_cast<double>(numSimulations);
    }

    std::pair<double, double> simulateStandardWithError(size_t numSimulations) {
        simulators_.clear();
        simulators_.reserve(numSimulations);

        for (size_t i{}; i < numSimulations; ++i) {
            simulators_.emplace_back(process_, scheme_, payoff_, initialValue_,
                                     timeToMaturity_, stepCount_);
        }

        std::vector<double> results;
        results.reserve(numSimulations);
        for (auto& simulator : simulators_) {
            results.push_back(payoff_(simulator.generateTerminal()));
        }

        double mean = std::accumulate(results.begin(), results.end(), 0.0) / 
                     static_cast<double>(numSimulations);

        double variance = 0.0;
        for (const auto& r : results)
            variance += (r - mean) * (r - mean);
        variance /= static_cast<double>(numSimulations - 1);
        double stdError = std::sqrt(variance / static_cast<double>(numSimulations));

        return {mean, stdError};
    }

    // Antithetic Variates Monte Carlo
    double simulateAntithetic(size_t numPairs) {
        simulators_.clear();
        simulators_.reserve(1);  // Only need one simulator
        simulators_.emplace_back(process_, scheme_, payoff_, initialValue_,
                                timeToMaturity_, stepCount_);

        auto& simulator = simulators_[0];
        double totalResult = 0.0;

        for (size_t i = 0; i < numPairs; ++i) {
            // Generate antithetic pair
            auto [randoms, antitheticRandoms] = simulator.getRng().generateAntitheticPair(stepCount_);

            // Run with original randoms
            double terminal1 = simulator.generateTerminalWithRandoms(randoms);
            double payoff1 = payoff_(terminal1);

            // Run with antithetic randoms
            double terminal2 = simulator.generateTerminalWithRandoms(antitheticRandoms);
            double payoff2 = payoff_(terminal2);

            // Average the pair
            totalResult += (payoff1 + payoff2) / 2.0;
        }

        return totalResult / static_cast<double>(numPairs);
    }

    std::pair<double, double> simulateAntitheticWithError(size_t numPairs) {
        simulators_.clear();
        simulators_.reserve(1);
        simulators_.emplace_back(process_, scheme_, payoff_, initialValue_,
                                timeToMaturity_, stepCount_);

        auto& simulator = simulators_[0];
        std::vector<double> pairedResults;
        pairedResults.reserve(numPairs);

        for (size_t i = 0; i < numPairs; ++i) {
            auto [randoms, antitheticRandoms] = simulator.getRng().generateAntitheticPair(stepCount_);

            double terminal1 = simulator.generateTerminalWithRandoms(randoms);
            double payoff1 = payoff_(terminal1);

            double terminal2 = simulator.generateTerminalWithRandoms(antitheticRandoms);
            double payoff2 = payoff_(terminal2);

            // Store the averaged pair
            pairedResults.push_back((payoff1 + payoff2) / 2.0);
        }

        double mean = std::accumulate(pairedResults.begin(), pairedResults.end(), 0.0) / 
                     static_cast<double>(numPairs);

        double variance = 0.0;
        for (const auto& r : pairedResults)
            variance += (r - mean) * (r - mean);
        variance /= static_cast<double>(numPairs - 1);
        double stdError = std::sqrt(variance / static_cast<double>(numPairs));

        return {mean, stdError};
    }

    P process_;
    S scheme_;
    PO payoff_;
    double initialValue_;
    double timeToMaturity_;
    size_t stepCount_;

    std::unique_ptr<VarianceReductionStrategy> currentStrategy_;
    std::vector<Simulator<P, S, PO>> simulators_;
};

} // namespace lfmc