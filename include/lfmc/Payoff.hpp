#pragma once
#include <algorithm>
#include <cmath>
#include <vector>

namespace lfmc {

/**
 * @brief European Call option payoff.
 * Payoff = max(S_T - K, 0)
 */
class EuropeanCall {
public:
    explicit EuropeanCall(double strike) : strike_(strike) {}
    
    double operator()(double terminal_value) const noexcept {
        return std::max(terminal_value - strike_, 0.0);
    }
    
    double strike() const noexcept { return strike_; }

private:
    double strike_;
};

/**
 * @brief European Put option payoff.
 * Payoff = max(K - S_T, 0)
 */
class EuropeanPut {
public:
    explicit EuropeanPut(double strike) : strike_(strike) {}
    
    double operator()(double terminal_value) const noexcept {
        return std::max(strike_ - terminal_value, 0.0);
    }
    
    double strike() const noexcept { return strike_; }

private:
    double strike_;
};

} // namespace lfmc