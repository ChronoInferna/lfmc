#pragma once

#include <concepts>

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

class VarianceReductionStrategy {
  public:
    virtual ~VarianceReductionStrategy() = 0;

    virtual double apply(double data) noexcept = 0;
};

inline VarianceReductionStrategy::~VarianceReductionStrategy() = default;

template <typename T>
concept VRStrategy = std::derived_from<T, VarianceReductionStrategy>;

class NoVarianceReduction : public VarianceReductionStrategy {
  public:
    NoVarianceReduction() noexcept = default;
    ~NoVarianceReduction() noexcept override = default;

    double apply(double data) noexcept override {
        return data;
    }
};

class AntitheticVariates : public VarianceReductionStrategy {
  public:
    AntitheticVariates() noexcept = default;
    ~AntitheticVariates() noexcept override = default;

    double apply(double data) noexcept override {
        // Placeholder implementation - in practice, this would need to be integrated with the
        // random number generation and path simulation to create antithetic pairs.
        return data;
        // Negate the data as a simple example of an antithetic transformation
    }
};

} // namespace lfmc
