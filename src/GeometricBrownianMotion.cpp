namespace lfmc {

class GeometricBrownianMotion {
  private:
    double mu_;
    double sigma_;
    double x0_;

  public:
    GeometricBrownianMotion(double mu, double sigma, double x0) : mu_(mu), sigma_(sigma), x0_(x0) {}

    double initial() const noexcept {
        return x0_;
    }

    double drift(double, double x) const noexcept {
        return mu_ * x;
    }

    double diffusion(double, double x) const noexcept {
        return sigma_ * x;
    }
};

} // namespace lfmc
