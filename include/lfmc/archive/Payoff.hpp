#pragma once

#include "types.hpp"

#include <algorithm>
#include <cmath>

namespace lfmc {

template <class P>
concept Payoff = requires(P const& p, const Path& x) {
    { p(x) } -> std::same_as<double>;
};

/**
 * @brief European Call option payoff.
 * Payoff = max(S_T - K, 0)
 */
struct EuropeanCall {
    double strike;

    double operator()(const Path& path) const noexcept {
        return std::max(path[path.size() - 1] - strike, 0.0);
    }
};

/**
 * @brief European Put option payoff.
 * Payoff = max(K - S_T, 0)
 */
struct EuropeanPut {
    double strike;

    double operator()(const Path& path) const noexcept {
        return std::max(strike - path[path.size() - 1], 0.0);
    }
};

} // namespace lfmc
