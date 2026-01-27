#pragma once

#include <concepts>

/**
 * @file StochasticProcess.hpp
 * @brief Defines the StochasticProcess concept for stochastic differential equations (SDEs) and
 * provides some implementations.
 */

namespace lfmc {

template <class P>
concept StochasticProcess = requires(P const& p, double x) {
    { p.drift(x) } -> std::same_as<double>;
    { p.diffusion(x) } -> std::same_as<double>;
};

struct GeometricBrownianMotion {
    double mu;
    double sigma;

    double drift(double x) const noexcept {
        return mu * x;
    }

    double diffusion(double x) const noexcept {
        return sigma * x;
    }
};

} // namespace lfmc
