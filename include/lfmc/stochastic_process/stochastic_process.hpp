#pragma once

#include <concepts>

/**
 * @file stochastic_process.hpp
 * @brief Defines the StochasticProcess concept for stochastic differential equations (SDEs).
 */

namespace lfmc {

template <typename P>
concept StochasticProcess = requires(P const& p, double t, double x) {
    { p.initial() } -> std::convertible_to<double>;
    { p.drift(t, x) } -> std::convertible_to<double>;
    { p.diffusion(t, x) } -> std::convertible_to<double>;
};

} // namespace lfmc
