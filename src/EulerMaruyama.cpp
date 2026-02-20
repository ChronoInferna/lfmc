#include "lfmc/StochasticProcess.hpp"

#include <cmath>

namespace lfmc {

class EulerMaruyama {
  public:
    template <StochasticProcess P>
    double step(P const& process, double t, double x, double dt, double z) const noexcept {
        double drift = process.drift(t, x);
        double diffusion = process.diffusion(t, x);
        return x + drift * dt + diffusion * std::sqrt(dt) * z;
    }
};

} // namespace lfmc
