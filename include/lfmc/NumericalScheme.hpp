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
 * A type `S` satisfies the `NumericalScheme` concept if it provides the following member
 * functions:
 * - `double step(double x, double dt, double dW) const`: Computes the next state given the current
 *   state `x`, time step `dt`, and Wiener increment `dW`.
 *
 * @tparam S The type to be checked against the concept.
 * @tparam P The stochastic process type used in the SDE.
 */
template <class S, class P>
concept NumericalScheme = requires(S const& s, P const& p, double x, double dt, double dW) {
    { s.step(x, dt, dW, p) } -> std::same_as<double>;
};

/**
 * @brief Euler-Maruyama numerical scheme for solving SDEs.
 *
 * The Euler-Maruyama method is a simple and widely used numerical scheme for approximating
 * solutions to stochastic differential equations. It is defined by the update rule:
 * X_{n+1} = X_n + drift(X_n) * dt + diffusion(X_n) * dW
 */
struct EulerMaruyama {
    /**
     * @brief Compute the next state using the Euler-Maruyama method.
     * @tparam P The stochastic process type used in the SDE.
     * @param process The stochastic process defining the SDE.
     * @param x The current state variable.
     * @param dt The time step size.
     * @param dW The Wiener increment.
     * @return The next state variable.
     */
    template <StochasticProcess P>
    double step(P const& process, double x, double dt, double dW) const noexcept {
        return x + process.drift(x) * dt + process.diffusion(x) * dW;
    }
};

} // namespace lfmc
