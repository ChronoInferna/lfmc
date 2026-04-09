#pragma once

#include "lfmc/stochastic_process/stochastic_process.hpp"

#include <cmath>

namespace lfmc {

// Euler-Maruyama discretisation of a generic SDE:
//   X_{t+dt} = X_t + drift(t, X_t)*dt + diffusion(t, X_t)*sqrt(dt)*Z
//
// DISCRETISATION BIAS FOR GBM
// For Geometric Brownian Motion (drift=mu*X, diffusion=sigma*X) this is an
// ARITHMETIC approximation of the exact LOG-NORMAL step:
//   Exact: X_{t+dt} = X_t * exp((mu - 0.5*sigma^2)*dt + sigma*sqrt(dt)*Z)
//   Euler: X_{t+dt} = X_t * (1 + mu*dt + sigma*sqrt(dt)*Z)
//
// The two distributions differ by O(dt) terms. Consequences:
//   - Terminal price distribution is not exactly log-normal.
//   - Option prices carry a systematic O(dt) discretisation bias.
//   - Measured bias for ATM European call (sigma=0.2, T=1, steps=1): ~-0.22 (~2%).
//   - The bias shrinks proportionally to 1/steps as steps increases.
//
// For paper benchmarks, use at least 52 weekly steps; 252 daily steps gives
// a bias well below statistical noise at typical sample counts.
// The GBMExact class (commented out below) eliminates this bias entirely.
template <StochasticProcess P> class EulerMaruyama {
  public:
    double step(P const& process, double t, double x, double dt, double z) const noexcept {
        double drift = process.drift(t, x);
        double diffusion = process.diffusion(t, x);
        return x + drift * dt + diffusion * std::sqrt(dt) * z;
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
