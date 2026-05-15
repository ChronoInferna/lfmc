#pragma once

/**
 * @file StochasticProcess.hpp
 * @brief Defines the StochasticProcess concept for stochastic differential equations (SDEs) and
 * provides some implementations.
 */

namespace lfmc {

class GeometricBrownianMotion {
  private:
    double mu_;
    double sigma_;
    double x0_;

  public:
    GeometricBrownianMotion(double mu, double sigma, double x0);

    double initial() const noexcept;
    double drift(double x, double) const noexcept;
    double diffusion(double x, double) const noexcept;

    double mu() const noexcept;
    double sigma() const noexcept;
};

} // namespace lfmc
