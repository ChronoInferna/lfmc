#pragma once

#include "NumericalScheme.hpp"
#include "StochasticProcess.hpp"
#include "types.hpp"

#include <span>

namespace lfmc {

template <StochasticProcess P, NumericalScheme<P> S> class PathGenerator {
  private:
    P process_;
    S scheme_;
    State state_;
    double dt_;

  public:
    PathGenerator(P process, S scheme, State state)
        : process_(std::move(process)), scheme_(std::move(scheme)), state_(state),
          dt_(state.timeToMaturity / static_cast<double>(state.steps)) {}

    Path generate(std::span<const double> randomNormals) {
        Path path(state_.steps + 1);
        path[0] = state_.initialValue;

        for (size_t i = 0; i < state_.steps; ++i) {
            path[i + 1] = scheme_.step(process_, path[i], dt_, randomNormals[i]);
        }

        return path;
    }
};

} // namespace lfmc
