#pragma once
#include "StochasticProcess.hpp"

#include <cmath>
#include <concepts>
#include <cstddef>
#include <utility>
#include <vector>

/**
 * @file NumericalScheme.hpp
 * @brief Defines the NumericalScheme concept for numerical methods solving SDEs and
 * provides some implementations.
 */

namespace lfmc {

template <typename S, typename P>
concept NumericalScheme = requires(S const& s, P const& p, double x, double dt, double dW) {
    { s.step(p, x, dt, dW) } -> std::same_as<double>;
};

template <StochasticProcess P> class EulerMaruyama {
  public:
    explicit EulerMaruyama(P process) : process_(std::move(process)) {}

    /**
     * @brief Compute the next state using Euler-Maruyama.
     * @param x_current Current state.
     * @param dt Time step size.
     * @param z Standard normal random variable N(0,1).
     * @return Next state X_{t+dt}.
     */
    double step(double x_current, double dt, double z) const noexcept {
        double drift = process_.drift(x_current);
        double diffusion = process_.diffusion(x_current);
        return x_current + drift * dt + diffusion * std::sqrt(dt) * z;
    }

    /**
     * @brief Simulate an entire path.
     * @param x0 Initial value.
     * @param T Total time.
     * @param n_steps Number of time steps.
     * @param random_normals Vector of N(0,1) random variables (size = n_steps).
     * @return Vector of path values [X_0, X_1, ..., X_T].
     */
    std::vector<double> simulate_path(double x0, double T, size_t n_steps,
                                      std::vector<double> const& random_normals) const {
        double dt = T / static_cast<float>(n_steps);
        std::vector<double> path;
        path.reserve(n_steps + 1);
        path.push_back(x0);

        double x = x0;
        for (size_t i = 0; i < n_steps; ++i) {
            x = step(x, dt, random_normals[i]);
            path.push_back(x);
        }
        return path;
    }

    /**
     * @brief Simulate only the terminal value (no full path).
     * @param x0 Initial value vector of N(0,1) random variables.
     * @return Terminal value X_T.
     */
    double simulate_terminal(double x0, double T, size_t n_steps,
                             std::vector<double> const& random_normals) const {
        double dt = T / static_cast<float>(n_steps);
        double x = x0;
        for (size_t i = 0; i < n_steps; ++i) {
            x = step(x, dt, random_normals[i]);
        }
        return x;
    }

  private:
    P process_;
};

/**
 * @brief Exact simulation for Geometric Brownian Motion.
 *
 * Uses the closed-form solution:
 * X_T = X_0 * exp((mu - 0.5*sigma^2)*T + sigma*sqrt(T)*Z)
 *
 * This is faster and more accurate than Euler-Maruyama for GBM.
 */
class GBMExact {
  public:
    explicit GBMExact(GeometricBrownianMotion gbm) : gbm_(gbm) {}

    /**
     * @brief Simulate terminal value using exact solution.
     * @param x0 Initial value.
     * @param T Time to maturity.
     * @param z Standard normal random variable.
     * @return Terminal value X_T.
     */
    double simulate_terminal(double x0, double T, double z) const noexcept {
        double drift_adjusted = (gbm_.mu - 0.5 * gbm_.sigma * gbm_.sigma) * T;
        double diffusion_term = gbm_.sigma * std::sqrt(T) * z;
        return x0 * std::exp(drift_adjusted + diffusion_term);
    }

  private:
    GeometricBrownianMotion gbm_;
};

} // namespace lfmc
