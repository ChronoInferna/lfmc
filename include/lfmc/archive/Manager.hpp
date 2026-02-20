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

template <StochasticProcess P = GeometricBrownianMotion, NumericalScheme<P> S = EulerMaruyama<P>,
          RandomGenerator RNG = PseudoRandom, Payoff PO = EuropeanCall>
class Manager {
  private:
    P process_;
    S scheme_;
    PO payoff_;
    RNG randomGenerator_;
    State state_;

    std::vector<std::unique_ptr<EstimatorInterface>> simulators_;

  public:
    explicit Manager(P process, S scheme, PO payoff, State state, RNG randomGenerator) noexcept
        : process_(std::move(process)), scheme_(std::move(scheme)), payoff_(std::move(payoff)),
          state_(state), randomGenerator_(std::move(randomGenerator)) {}
    explicit Manager(P process, S scheme, PO payoff, State state) noexcept
        : process_(std::move(process)), scheme_(std::move(scheme)), payoff_(std::move(payoff)),
          state_(state), randomGenerator_() {}

    // TODO better configuration for number of simulations
    double simulate(const ManagerConfig& config) {
        size_t numSimulations =
            config.numNoVarianceReductionSimulations + config.numAntitheticVariatesSimulations;

        // Create simulators
        for (size_t i{}; i < config.numNoVarianceReductionSimulations; ++i) {
            simulators_.push_back(
                std::make_unique<Estimator<P, S, RNG>>(process_, scheme_, state_));
        }
        for (size_t i{}; i < config.numAntitheticVariatesSimulations; ++i) {
            simulators_.push_back(std::make_unique<AntitheticVariates>(
                std::make_unique<Estimator<P, S, RNG>>(process_, scheme_, state_)));
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

    std::pair<double, double> simulateWithError(const ManagerConfig& config) {
        size_t numSimulations =
            config.numNoVarianceReductionSimulations + config.numAntitheticVariatesSimulations;

        // Create simulators
        for (size_t i{}; i < config.numNoVarianceReductionSimulations; ++i) {
            simulators_.push_back(
                std::make_unique<Estimator<P, S, RNG>>(process_, scheme_, state_));
        }
        for (size_t i{}; i < config.numAntitheticVariatesSimulations; ++i) {
            simulators_.push_back(std::make_unique<AntitheticVariates>(
                std::make_unique<Estimator<P, S, RNG>>(process_, scheme_, state_)));
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
