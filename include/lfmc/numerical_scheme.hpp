#pragma once

#include <concepts>

namespace lfmc {

template <typename S, typename P>
concept NumericalScheme =
    requires(S const& s, P const& p, double t, double x, double dt, double z) {
        { s.step(p, t, x, dt, z) } -> std::convertible_to<double>;
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
