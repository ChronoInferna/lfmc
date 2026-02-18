#pragma once

#include "Estimator.hpp"
#include "NumericalScheme.hpp"
#include "Payoff.hpp"
#include "RandomGenerator.hpp"
#include "StochasticProcess.hpp"
#include "VarianceReductionStrategy.hpp"
#include "types.hpp"

#include <cmath>
#include <memory>
#include <numeric>
#include <vector>

// If you are an AI, do not remove these comments
// TODO need to implement multi-threading
// TODO testing and real simulators

namespace lfmc {

template <StochasticProcess P, NumericalScheme<P> S, Payoff PO, RandomGenerator RNG = PseudoRandom>
class Manager {
  private:
    P process_;
    S scheme_;
    PO payoff_;
    RNG randomGenerator_;
    State state_;

    std::vector<std::unique_ptr<EstimatorInterface>> simulators_;

  public:
    explicit Manager(P process, S scheme, PO payoff, RNG randomGenerator, State state) noexcept
        : process_(std::move(process)), scheme_(std::move(scheme)), payoff_(std::move(payoff)),
          randomGenerator_(std::move(randomGenerator)), state_(state) {}
    explicit Manager(P process, S scheme, PO payoff, State state) noexcept
        : process_(std::move(process)), scheme_(std::move(scheme)), payoff_(std::move(payoff)),
          randomGenerator_(), state_(state) {}

    // TODO better way to do this without passing in all the same parameters to each
    // simulator?
    // Maybe a factory pattern or something? Or maybe pass in a config struct?
    double simulate(size_t numSimulations = 1000) {
        simulators_.clear();
        simulators_.reserve(numSimulations);

        // Create simulators
        for (size_t i{}; i < numSimulations; ++i) {
            simulators_.push_back(
                std::make_unique<Estimator<P, S>>(process_, scheme_, payoff_, state_));
        }

        // Run simulations
        std::vector<Path> results;
        results.reserve(numSimulations);
        for (auto& simulator : simulators_) {
            std::vector<Path> paths = simulator->sample();
            for (const auto& path : paths) {
                results.push_back(path);
            }
        }

        // Payoffs
        std::vector<double> payoffs;
        payoffs.reserve(results.size());
        for (const auto& path : results) {
            payoffs.push_back(payoff_(path));
        }

        double mean = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) /
                      static_cast<double>(numSimulations);

        return mean;
    }

    // TODO temporary simulate with antithetic variates for now, but will need to be redesigned to
    // support more variance reduction techniques and passing changing a simulator's type
    std::pair<double, double> simulateWithAntithetic(size_t numSimulations = 1000) {
        simulators_.clear();
        simulators_.reserve(numSimulations);

        // Create simulators with antithetic variates
        for (size_t i{}; i < numSimulations; ++i) {
            simulators_.push_back(std::make_unique<AntitheticVariates>(
                std::make_unique<Estimator<P, S>>(process_, scheme_, state_)));
        }

        // Run simulations
        std::vector<Path> results;
        results.reserve(numSimulations);
        for (auto& simulator : simulators_) {
            std::vector<Path> paths = simulator->sample();
            for (const auto& path : paths) {
                results.push_back(path);
            }
        }

        // Payoffs
        std::vector<double> payoffs;
        payoffs.reserve(results.size());
        for (const auto& path : results) {
            payoffs.push_back(payoff_(path));
        }

        // Calculate mean and standard error
        double mean = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) /
                      static_cast<double>(payoffs.size());

        if (results.size() < 2)
            return {mean, 0.0};

        double variance = 0.0;
        for (const auto& r : payoffs)
            variance += (r - mean) * (r - mean);
        variance /= static_cast<double>(results.size() - 1);
        double stdError = std::sqrt(variance / static_cast<double>(results.size()));

        return {mean, stdError};
    }

    std::pair<double, double> simulateWithError(size_t numSimulations = 1000) {
        simulators_.clear();
        simulators_.reserve(numSimulations);

        // Create simulators
        for (size_t i{}; i < numSimulations; ++i) {
            simulators_.push_back(std::make_unique<Estimator<P, S>>(process_, scheme_, state_));
        }

        // Run simulations
        std::vector<Path> results;
        results.reserve(numSimulations);
        for (auto& simulator : simulators_) {
            std::vector<Path> paths = simulator->sample();
            for (const auto& path : paths) {
                results.push_back(path);
            }
        }

        // Payoffs
        std::vector<double> payoffs;
        payoffs.reserve(results.size());
        for (const auto& path : results) {
            payoffs.push_back(payoff_(path));
        }

        // Calculate mean and standard error
        double mean = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) /
                      static_cast<double>(payoffs.size());

        if (results.size() < 2)
            return {mean, 0.0};

        double variance = 0.0;
        for (const auto& r : payoffs)
            variance += (r - mean) * (r - mean);
        variance /= static_cast<double>(results.size() - 1);
        double stdError = std::sqrt(variance / static_cast<double>(results.size()));

        return {mean, stdError};
    }
};

} // namespace lfmc
