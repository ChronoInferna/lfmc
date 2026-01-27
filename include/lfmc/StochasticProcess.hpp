#pragma once

#include <concepts>

/**
 * @file StochasticProcess.hpp
 * @brief Defines the StochasticProcess concept for stochastic differential equations (SDEs) and
 * provides some implementations.
 */

namespace lfmc {

/**
 * @brief Concept for stochastic processes used in SDEs.
 *
 * A type `P` satisfies the `StochasticProcess` concept if it provides the following member
 * functions:
 * - `double drift(double x) const`: Computes the drift term at state `x`.
 * - `double diffusion(double x) const`: Computes the diffusion term at state `x`.
 *
 * @tparam P The type to be checked against the concept.
 * @param x The current state variable.
 */
template <class P>
concept StochasticProcess = requires(P const& p, double x) {
    { p.drift(x) } -> std::same_as<double>;
    { p.diffusion(x) } -> std::same_as<double>;
};

/**
 * @brief Geometric Brownian Motion (GBM) stochastic process.
 *
 * The GBM process is defined by the stochastic differential equation:
 * dX_t = mu * X_t * dt + sigma * X_t * dW_t
 * where `mu` is the drift coefficient and `sigma` is the diffusion coefficient.
 *
 * @param mu The drift coefficient.
 * @param sigma The diffusion coefficient.
 */
struct GeometricBrownianMotion {
    double mu;
    double sigma;

    /**
     * @brief Compute drift term at state x.
     * @param x The current state variable.
     * @return The drift term mu * x.
     */
    double drift(double x) const noexcept {
        return mu * x;
    }

    /**
     * @brief Compute diffusion term at state x.
     * @param x The current state variable.
     * @return The diffusion term sigma * x.
     */

    double diffusion(double x) const noexcept {
        return sigma * x;
    }
};

} // namespace lfmc
