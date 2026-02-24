#pragma once

#include <concepts>

/**
 * @file StochasticProcess.hpp
 * @brief Defines the StochasticProcess concept for stochastic differential equations (SDEs) and
 * provides some implementations.
 */

namespace lfmc {

template <typename P>
concept StochasticProcess = requires(P const& p, double t, double x) {
    { p.initial() } -> std::convertible_to<double>;
    { p.drift(t, x) } -> std::convertible_to<double>;
    { p.diffusion(t, x) } -> std::convertible_to<double>;
};

class GeometricBrownianMotion {
  private:
    double mu_;
    double sigma_;
    double x0_;

  public:
    GeometricBrownianMotion(double mu, double sigma, double x0);

    double initial() const noexcept;
    double drift(double, double x) const noexcept;
    double diffusion(double, double x) const noexcept;
};

} // namespace lfmc
