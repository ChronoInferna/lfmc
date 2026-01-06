#pragma once

// TODO: Implement derived classes for specific variance reduction techniques.
// Examples include Antithetic Variates, Control Variates, Importance Sampling, etc.
// First must decide how to design parallel infrastructure to support these techniques in Monte
// Carlo simulations.

/**
 * @file VarianceReductionStrategy.hpp
 * @brief Defines the VarianceReductionStrategy base class for runtime polymorphism of variance
 * reduction techniques and provides some implementations.
 */

namespace lfmc {

/**
 * @brief Base class for variance reduction strategies.
 *
 * This abstract class defines the interface for different variance reduction techniques.
 * Derived classes must implement the `apply` method to modify the input data accordingly.
 */
class VarianceReductionStrategy {
  public:
    /// @brief Pure virtual destructor to ensure proper cleanup of derived classes.
    virtual ~VarianceReductionStrategy() = 0;

    /**
     * @brief Apply the variance reduction technique to the input data.
     * @param data The input data to be modified.
     * @return The modified data after applying the variance reduction technique.
     */
    virtual double apply(double data) const noexcept = 0;
};

/// @brief Definition of the pure virtual destructor.
inline VarianceReductionStrategy::~VarianceReductionStrategy() = default;

class NoVarianceReduction : public VarianceReductionStrategy {
  public:
    NoVarianceReduction() noexcept = default;
    ~NoVarianceReduction() noexcept override = default;

    double apply(double data) const noexcept override {
        return data; // No modification
    }
};

} // namespace lfmc
