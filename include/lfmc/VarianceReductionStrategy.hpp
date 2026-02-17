#pragma once

#include "lfmc/Simulator.hpp"

#include <memory>

// TODO: Implement derived classes for specific variance reduction techniques.
// Examples include Antithetic Variates, Control Variates, Importance Sampling, etc.
// First must decide how to design parallel infrastructure to support these techniques in Monte
// Carlo simulations - decorator design pattern?
// TODO each strategy has a window parameter for how much data you're using

/**
 * @file VarianceReductionStrategy.hpp
 * @brief Defines the VarianceReductionStrategy base class for runtime polymorphism of variance
 * reduction techniques and provides some implementations.
 */

namespace lfmc {

// Decorator for variance reduction strategies (e.g., Antithetic Variates, Control Variates, etc.)
class VarianceReductionBaseDecorator : public SimulatorInterface {
  protected:
    std::unique_ptr<SimulatorInterface> simulator_; // Pointer to the base simulator

  public:
    explicit VarianceReductionBaseDecorator(std::unique_ptr<SimulatorInterface> simulator)
        : simulator_(std::move(simulator)) {}

    double sample() override {
        // By default, just call the underlying simulator's sample method
        return simulator_->sample();
    }

    // TODO temporary now, see note in Simulator.hpp
    double sampleFromRandoms(const std::vector<double>& Z) override {
        return simulator_->sampleFromRandoms(Z);
    }
    size_t steps() const noexcept override {
        return simulator_->steps();
    }
    RandomGenerator& rng() override {
        return simulator_->rng();
    }
};

// Example implementation of Antithetic Variates strategy
class AntitheticVariates : public VarianceReductionBaseDecorator {
  public:
    using Base = VarianceReductionBaseDecorator;
    AntitheticVariates(std::unique_ptr<SimulatorInterface> simulator)
        : Base(std::move(simulator)) {}

    double sample() override {
        size_t steps = simulator_->steps();
        auto randoms = simulator_->rng().generateNormals(steps);
        std::vector<double> antitheticRandoms(steps);
        std::transform(randoms.begin(), randoms.end(), antitheticRandoms.begin(),
                       [](double r) { return -r; }); // Create antithetic randoms

        double payoff1 = simulator_->sampleFromRandoms(randoms);
        double payoff2 = simulator_->sampleFromRandoms(antitheticRandoms);

        return (payoff1 + payoff2) / 2.0;
    }
};

} // namespace lfmc
