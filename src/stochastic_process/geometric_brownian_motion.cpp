#include "lfmc/stochastic_process.hpp"

namespace lfmc {

GeometricBrownianMotion::GeometricBrownianMotion(double mu, double sigma, double x0)
    : mu_(mu), sigma_(sigma), x0_(x0) {}

double GeometricBrownianMotion::initial() const noexcept {
    return x0_;
}

double GeometricBrownianMotion::drift(double, double x) const noexcept {
    return mu_ * x;
}

double GeometricBrownianMotion::diffusion(double, double x) const noexcept {
    return sigma_ * x;
}

double GeometricBrownianMotion::mu() const noexcept {
    return mu_;
}

double GeometricBrownianMotion::sigma() const noexcept {
    return sigma_;
}

} // namespace lfmc
