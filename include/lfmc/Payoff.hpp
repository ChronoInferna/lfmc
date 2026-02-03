#pragma once

#include <algorithm>
#include <cmath>

namespace lfmc {

template <class P>
concept Payoff = requires(P const& p, double x) {
    { p(x) } -> std::same_as<double>;
};

/**
 * @brief European Call option payoff.
 * Payoff = max(S_T - K, 0)
 */
struct EuropeanCall {
    double strike;

    double operator()(double terminal_value) const noexcept {
        return std::max(terminal_value - strike, 0.0);
    }
};

/**
 * @brief European Put option payoff.
 * Payoff = max(K - S_T, 0)
 */
class EuropeanPut {
    double strike;

    double operator()(double terminal_value) const noexcept {
        return std::max(strike - terminal_value, 0.0);
    }
};

} // namespace lfmc
