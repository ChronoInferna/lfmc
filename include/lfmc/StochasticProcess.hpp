#pragma once
#include <stdexcept>
#include <concepts>

/**
 * @file StochasticProcess.hpp
 * @brief Defines the StochasticProcess concept for stochastic differential equations (SDEs) and
 * provides some implementations.
 */

namespace lfmc {

template <class P>
concept StochasticProcess = requires(P const& p, double x) {
    { p.drift(x) } -> std::same_as<double>;
    { p.diffusion(x) } -> std::same_as<double>;
};

struct GeometricBrownianMotion {
    double mu;
    double sigma;

    /**
     * @brief Constructor to initialize GBM parameters.
     * @param drift_coef The drift coefficient (mu).
     * @param diffusion_coef The diffusion coefficient (sigma/volatility).
     */
    GeometricBrownianMotion(double drift_coef, double diffusion_coef)
        : mu(drift_coef), sigma(diffusion_coef) {
        if (sigma < 0.0) {
            throw std::invalid_argument("Diffusion coefficient (sigma) must be non-negative");
        }
    }

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