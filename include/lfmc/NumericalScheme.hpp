#pragma once
#include "StochasticProcess.hpp"

#include <concepts>

/**
 * @file NumericalScheme.hpp
 * @brief Defines the NumericalScheme concept for numerical methods solving SDEs and
 * provides some implementations.
 */

namespace lfmc {

template <typename S, typename P>
concept NumericalScheme = requires(S const& s, P const& p, double x, double dt, double dW) {
    { s.step(p, x, dt, dW) } -> std::same_as<double>;
};

template <StochasticProcess P> struct EulerMaruyama {
    double step(P const& process, double x, double dt, double dW) const noexcept {
        return x + process.drift(x) * dt + process.diffusion(x) * dW;
    }
};

} // namespace lfmc
