#include "lfmc/payoff.hpp"

#include <cmath>

namespace lfmc {

EuropeanCall::EuropeanCall(double strike) : strike_(strike) {}

Payoffs EuropeanCall::generate_payoffs(const std::vector<Path>& paths) const {
    Payoffs payoffs;
    for (const auto& path : paths) {
        double final_price = path.back();
        payoffs.push_back(std::max(final_price - strike_, 0.0));
    }
    return payoffs;
}

EuropeanPut::EuropeanPut(double strike) : strike_(strike) {}

Payoffs EuropeanPut::generate_payoffs(const std::vector<Path>& paths) const {
    Payoffs payoffs;
    for (const auto& path : paths) {
        double final_price = path.back();
        payoffs.push_back(std::max(strike_ - final_price, 0.0));
    }
    return payoffs;
}

} // namespace lfmc
