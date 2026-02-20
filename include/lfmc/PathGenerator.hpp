#pragma once

#include "lfmc/NumericalScheme.hpp"
#include "lfmc/StochasticProcess.hpp"
#include "lfmc/types.hpp"

namespace lfmc {

template <StochasticProcess Process, typename Scheme>
    requires NumericalScheme<Scheme, Process>
class PathGenerator {
  private:
    Process process_;
    Scheme scheme_;
    double T_;
    size_t steps_;

  public:
    PathGenerator(Process process, Scheme scheme, double T, size_t steps)
        : process_(std::move(process)), scheme_(std::move(scheme)), T_(T), steps_(steps) {}

    Path generate(const Normals& normals) const {
        const double dt = T_ / steps_;

        Path path(steps_ + 1);

        double t = 0.0;
        double x = process_.initial();
        path.push_back(x);

        for (size_t i = 0; i < steps_; ++i) {
            x = scheme_.step(process_, t, x, dt, normals[i]);
            path.push_back(x);
            t += dt;
        }

        return path;
    }
};

} // namespace lfmc
