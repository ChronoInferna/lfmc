#include "lfmc/stochastic_process.hpp"

#include <cmath>

namespace lfmc {

template <StochasticProcess P = GeometricBrownianMotion> class EulerMaruyama {
  public:
    double step(P const& process, double t, double x, double dt, double z) const noexcept {
        double drift = process.drift(t, x);
        double diffusion = process.diffusion(t, x);
        return x + drift * dt + diffusion * std::sqrt(dt) * z;
    }
};

} // namespace lfmc
