#pragma once

#include <concepts>

namespace lfmc {

template <typename S, typename P>
concept NumericalScheme =
    requires(S const& s, P const& p, double t, double x, double dt, double z) {
        { s.step(p, t, x, dt, z) } -> std::convertible_to<double>;
    };

} // namespace lfmc
