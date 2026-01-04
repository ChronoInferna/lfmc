#pragma once
#include "StochasticProcess.hpp"

#include <concepts>

/**
 * @file NumericalScheme.hpp
 * @brief Defines the NumericalScheme concept for numerical methods solving SDEs and
 * provides some implementations.
 */

namespace lfmc {

/**
 * @brief Concept for numerical schemes used to solve SDEs.
 *
 * A type `S` satisfies the `NumericalScheme` concept for a stochastic process `P` if it provides
 * the following member function:
 * - `double step(P const& process, double x, double dt, double dW) const`: Computes the next state
 *   given the current state `x`, time step `dt`, and Wiener increment `dW`.
 *
 * @tparam S The numerical scheme type to be checked against the concept.
 * @tparam P The stochastic process type used in the SDE.
 * @param process The stochastic process defining the SDE.
 * @param x The current state variable.
 * @param dt The time step size.
 * @param dW The Wiener increment.
 */
template <class S, class P>
concept NumericalScheme =
    StochasticProcess<P> && requires(S const& s, P const& p, double x, double dt, double dW) {
        { s.step(p, x, dt, dW) } -> std::same_as<double>;
    };

/**
 * @brief Euler-Maruyama numerical scheme for solving SDEs.
 *
 * @tparam P The stochastic process type used in the SDE.
 *
 * The Euler-Maruyama method is a simple and widely used numerical scheme for approximating
 * solutions to stochastic differential equations. It is defined by the update rule:
 * X_{n+1} = X_n + drift(X_n) * dt + diffusion(X_n) * dW
 */
template <StochasticProcess P> struct EulerMaruyama {
    /**
     * @brief Compute the next state using the Euler-Maruyama method.
     * @param process The stochastic process defining the SDE.
     * @param x The current state variable.
     * @param dt The time step size.
     * @param dW The Wiener increment.
     * @return The next state variable.
     */
    double step(P const& process, double x, double dt, double dW) const noexcept {
        return x + process.drift(x) * dt + process.diffusion(x) * dW;
    }
};

} // namespace lfmc
