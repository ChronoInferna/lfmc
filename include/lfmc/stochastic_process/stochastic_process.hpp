#pragma once

#include <concepts>

/**
 * @file StochasticProcess.hpp
 * @brief Defines the StochasticProcess concept for stochastic differential equations (SDEs) and
 * provides some implementations.
 */

namespace lfmc {

template <typename P>
concept StochasticProcess = requires(P const& p, double t, double x) {
    { p.initial() } -> std::convertible_to<double>;
    { p.drift(x, t) } -> std::convertible_to<double>;
    { p.diffusion(x, t) } -> std::convertible_to<double>;
};

} // namespace lfmc
