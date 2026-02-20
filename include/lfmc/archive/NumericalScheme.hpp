#pragma once

#include "StochasticProcess.hpp"

#include <cmath>
#include <concepts>

/**
 * @file NumericalScheme.hpp
 * @brief Defines the NumericalScheme concept for numerical methods solving SDEs and
 * provides some implementations.
 */

namespace lfmc {

/**
 * @brief Concept for numerical schemes solving SDEs.
 *
 * A NumericalScheme must implement a step function that computes the next state
 * given the current state, time step, and a standard normal random variable.
 *
 * @tparam S Numerical scheme type.
 * @tparam P Stochastic process type.
 *
 * Requires:
 * - S must have a method:
 *   double step(const P& process, double x_current, double dt, double z) const noexcept;
 *   where:
 *     - process: Stochastic process defining drift and diffusion.
 *     - x_current: Current state.
 *     - dt: Time step size.
 *     - z: Standard normal random variable N(0,1).
 *   The method returns the next state X_{t+dt}.
 */
template <class S, class P>
concept NumericalScheme =
    StochasticProcess<P> && requires(S const& s, P const& p, double x, double dt, double z) {
        { s.step(p, x, dt, z) } -> std::same_as<double>;
    };

template <StochasticProcess P> struct EulerMaruyama {
    /**
     * @brief Compute the next state using Euler-Maruyama.
     * @param Stochastic process defining drift and diffusion.
     * @param x_current Current state.
     * @param dt Time step size.
     * @param z Standard normal random variable N(0,1).
     * @return Next state X_{t+dt}.
     */
    double step(const P& process, double x_current, double dt, double z) const noexcept {
        double drift = process.drift(x_current);
        double diffusion = process.diffusion(x_current);
        return x_current + drift * dt + diffusion * std::sqrt(dt) * z;
    }
};

/**
 * @brief Exact simulation for Geometric Brownian Motion.
 *
 * Uses the closed-form solution:
 * X_T = X_0 * exp((mu - 0.5*sigma^2)*T + sigma*sqrt(T)*Z)
 *
 * This is faster and more accurate than Euler-Maruyama for GBM.
 */
// class GBMExact {
//   public:
//     explicit GBMExact(GeometricBrownianMotion gbm) : gbm_(gbm) {}
//
//     /**
//      * @brief Simulate terminal value using exact solution.
//      * @param x0 Initial value.
//      * @param T Time to maturity.
//      * @param z Standard normal random variable.
//      * @return Terminal value X_T.
//      */
//     double simulate_terminal(double x0, double T, double z) const noexcept {
//         double drift_adjusted = (gbm_.mu - 0.5 * gbm_.sigma * gbm_.sigma) * T;
//         double diffusion_term = gbm_.sigma * std::sqrt(T) * z;
//         return x0 * std::exp(drift_adjusted + diffusion_term);
//     }
//
//   private:
//     GeometricBrownianMotion gbm_;
// };

} // namespace lfmc
