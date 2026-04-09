#pragma once

#include "lfmc/core/types.hpp"
#include "lfmc/numerical_scheme/numerical_scheme.hpp"
#include "lfmc/stochastic_process/stochastic_process.hpp"

namespace lfmc {

template <StochasticProcess Process, typename Scheme>
    requires NumericalScheme<Scheme, Process>
class PathGenerator {
  private:
    Process process_;
    Scheme scheme_;

  public:
    PathGenerator(Process process, Scheme scheme)
        : process_(std::move(process)), scheme_(std::move(scheme)) {}

    // TODO move to cpp file
    std::vector<Path> generate_paths(const std::vector<Normals>& normals, size_t steps,
                                     double T) const {
        const double dt = T / static_cast<double>(steps);

        std::vector<Path> paths;
        for (const auto& norm : normals) {
            Path path;
            path.reserve(steps + 1);

            double t = 0.0;
            double x = process_.initial();
            path.push_back(x);

            for (size_t i = 0; i < steps; ++i) {
                x = scheme_.step(process_, t, x, dt, norm[i]);
                path.push_back(x);
                t += dt;
            }
            paths.push_back(std::move(path));
        }

        return paths;
    }
};

} // namespace lfmc
