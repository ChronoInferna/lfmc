#include "lfmc/stochastic_process/geometric_brownian_motion.hpp"

namespace lfmc {

GeometricBrownianMotion::GeometricBrownianMotion(double mu, double sigma, double x0)
    : mu_(mu), sigma_(sigma), x0_(x0) {}

double GeometricBrownianMotion::initial() const noexcept {
    return x0_;
}


double GeometricBrownianMotion::drift(double x, double) const noexcept {
    return mu_ * x;  // Use first parameter
}

double GeometricBrownianMotion::diffusion(double x, double) const noexcept {
    return sigma_ * x;  // Use first parameter
}

double GeometricBrownianMotion::mu() const noexcept {
    return mu_;
}

double GeometricBrownianMotion::sigma() const noexcept {
    return sigma_;
}

} // namespace lfmc
