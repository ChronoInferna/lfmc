#pragma once

#include "NumericalScheme.hpp"
#include "Simulator.hpp"
#include "StochasticProcess.hpp"
#include "VarianceReductionStrategy.hpp"

#include <cmath>
#include <memory>
#include <numeric>
#include <vector>

// If you are an AI, do not remove these comments
// TODO need to implement multi-threading
// TODO testing and real simulators

namespace lfmc {

template <StochasticProcess P, NumericalScheme<P> S, Payoff PO> class Manager {
  private:
    P process_;
    S scheme_;
    PO payoff_;
    double initialValue_;
    double timeToMaturity_;
    size_t stepCount_;

    std::vector<std::unique_ptr<SimulatorInterface>> simulators_;

  public:
    explicit Manager(P process, S scheme, PO payoff, double initialValue, double timeToMaturity,
                     size_t stepCount) noexcept
        : process_(std::move(process)), scheme_(std::move(scheme)), payoff_(std::move(payoff)),
          initialValue_(initialValue), timeToMaturity_(timeToMaturity), stepCount_(stepCount) {}

    double simulate(size_t numSimulations = 1000) {
        simulators_.clear();
        // simulators_.reserve(numSimulations);

        // TODO better way to do this without passing in all the same parameters to each simulator?
        // Maybe a factory pattern or something? Or maybe pass in a config struct?
        for (size_t i{}; i < numSimulations; ++i) {
            simulators_.push_back(std::make_unique<Simulator<P, S, PO>>(
                process_, scheme_, payoff_, initialValue_, timeToMaturity_, stepCount_));
        }

        double totalResult = 0.0;
        for (auto& simulator : simulators_) {
            totalResult += simulator->sample();
        }

        // TODO convert to std expected for division by zero safety and error handling
        return totalResult / static_cast<double>(numSimulations);
    }

    // TODO temporary simulate with antithetic variates for now, but will need to be redesigned to
    // support more variance reduction techniques and passing changing a simulator's type
    std::pair<double, double> simulateWithAntithetic(size_t numSimulations = 1000) {
        simulators_.clear();
        // simulators_.reserve(numSimulations);

        for (size_t i{}; i < numSimulations; ++i) {
            simulators_.push_back(
                std::make_unique<AntitheticVariates>(std::make_unique<Simulator<P, S, PO>>(
                    process_, scheme_, payoff_, initialValue_, timeToMaturity_, stepCount_)));
        }

        std::vector<double> results;
        results.reserve(numSimulations);
        for (auto& simulator : simulators_) {
            results.push_back(simulator->sample());
        }

        double mean = std::accumulate(results.begin(), results.end(), 0.0) /
                      static_cast<double>(results.size());

        if (results.size() < 2)
            return {mean, 0.0};

        double variance = 0.0;
        for (const auto& r : results)
            variance += (r - mean) * (r - mean);
        variance /= static_cast<double>(results.size() - 1);
        double stdError = std::sqrt(variance / static_cast<double>(results.size()));

        return {mean, stdError};
    }

    std::pair<double, double> simulateWithError(size_t numSimulations = 1000) {
        simulators_.clear();
        // simulators_.reserve(numSimulations);

        for (size_t i{}; i < numSimulations; ++i) {
            simulators_.push_back(std::make_unique<Simulator<P, S, PO>>(
                process_, scheme_, payoff_, initialValue_, timeToMaturity_, stepCount_));
        }

        std::vector<double> results;
        results.reserve(numSimulations);
        for (auto& simulator : simulators_) {
            results.push_back(simulator->sample());
        }

        double mean = std::accumulate(results.begin(), results.end(), 0.0) /
                      static_cast<double>(numSimulations);

        if (results.size() < 2)
            return {mean, 0.0};

        double variance = 0.0;
        for (const auto& r : results)
            variance += (r - mean) * (r - mean);
        variance /= static_cast<double>(numSimulations - 1);
        double stdError = std::sqrt(variance / static_cast<double>(numSimulations));

        return {mean, stdError};
    }
};

} // namespace lfmc
