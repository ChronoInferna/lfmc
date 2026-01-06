#pragma once
#include "NumericalScheme.hpp"
#include "StochasticProcess.hpp"
#include "VarianceReductionStrategy.hpp"

// #include <expected> TODO
#include <memory>

/**
 * @file manager.hpp
 * @brief Manager class for variation reduction strategies using the Strategy Pattern.
 *
 * This file defines the `Manager` class, which is responsible for managing
 * different strategies and potentially coordinating between them during runtime.
 */

namespace lfmc {

/**
 * @brief Manager class for handling stochastic processes, numerical schemes,
 *        and variance reduction strategies.
 *
 * The `Manager` class utilizes the Strategy Pattern to manage different stochastic processes and
 * numerical schemes at compile-time, and variance reduction techniques at runtime. It holds
 * instances of a stochastic process, a numerical scheme, and a variance reduction strategy.
 *
 * @tparam P The stochastic process type.
 * @tparam S The numerical scheme type.
 */
template <StochasticProcess P, NumericalScheme<P> S> class Manager {
  public:
    /**
     * @brief Construct a Manager with the given stochastic process, numerical scheme,
     *        and optional variance reduction strategy.
     *
     * @param process The stochastic process instance.
     * @param scheme The numerical scheme instance.
     * @param strategy Optional variance reduction strategy instance.
     */
    explicit Manager(P process, S scheme,
                     std::unique_ptr<VarianceReductionStrategy> strategy) noexcept
        : process_(std::move(process)), scheme_(std::move(scheme)), strategy_(std::move(strategy)) {
    }

    /**
     * @brief Set the variance reduction strategy at runtime.
     * @param strategy The new variance reduction strategy to be used.
     */
    void setStrategy(std::unique_ptr<VarianceReductionStrategy>&& strategy) noexcept {
        strategy_ = std::move(strategy);
    }

    /**
     * @brief Simulate a single step of the stochastic process using the numerical scheme
     *        and apply the variance reduction strategy.
     *
     * @param x The current state variable.
     * @param dt The time step size.
     * @param dW The Wiener increment.
     * @return The next state variable after applying the variance reduction strategy,
     *         or an error message if the strategy is not set.
     */
    double simulateStep(double x, double dt, double dW) const noexcept {
        double nextX = scheme_.step(process_, x, dt, dW);
        double reducedX = strategy_->apply(nextX);
        return reducedX;
    }

    // TODO big 5? do we even need it?

  private:
    P process_;
    S scheme_;
    std::unique_ptr<VarianceReductionStrategy> strategy_;
};

} // namespace lfmc
