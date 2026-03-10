#pragma once

#include "lfmc/stochastic_process/stochastic_process.hpp"

#include <cmath>

namespace lfmc {

template <StochasticProcess P> class EulerMaruyama {
  public:
    double step(P const& process, double x, double t, double dt, double z) const noexcept {
        double drift = process.drift(x, t);
        double diffusion = process.diffusion(x, t);
        return x + drift * dt + diffusion * std::sqrt(dt) * z;
    }
};
} // namespace lfmc
