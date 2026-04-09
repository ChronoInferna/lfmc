#pragma once

namespace lfmc {

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

    double mu() const noexcept;
    double sigma() const noexcept;
};

} // namespace lfmc
