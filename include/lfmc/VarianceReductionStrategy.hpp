#pragma once

#include "Estimator.hpp"
#include "types.hpp"

#include <functional>
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
class VarianceReductionBaseDecorator : public EstimatorInterface {
  protected:
    std::unique_ptr<EstimatorInterface> estimator_; // Pointer to the base estimator

  public:
    explicit VarianceReductionBaseDecorator(std::unique_ptr<EstimatorInterface> estimator)
        : estimator_(std::move(estimator)) {}

    std::vector<Path> sample() override {
        // By default, just call the underlying estimator's sample method
        return estimator_->sample();
    }

    State const& getState() const override {
        return estimator_->getState();
    }

    Normals generateNormals(size_t n) override {
        return estimator_->generateNormals(n);
    }

    Path generatePath(std::span<const double> randomNormals) override {
        return estimator_->generatePath(randomNormals);
    }
};

// Example implementation of Antithetic Variates strategy
class AntitheticVariates : public VarianceReductionBaseDecorator {
  public:
    using Base = VarianceReductionBaseDecorator;
    AntitheticVariates(std::unique_ptr<EstimatorInterface> estimator)
        : Base(std::move(estimator)) {}

    std::vector<Path> sample() override {
        const State& state = estimator_->getState();

        Normals normals = estimator_->generateNormals(state.steps);
        Normals antitheticNormals(normals);
        std::transform(normals.begin(), normals.end(), antitheticNormals.begin(), std::negate());

        Path path = estimator_->generatePath(normals);
        Path antitheticPath = estimator_->generatePath(antitheticNormals);

        return {path, antitheticPath};
    }
};

} // namespace lfmc
