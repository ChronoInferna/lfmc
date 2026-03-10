#include "lfmc/payoff/european_payoffs.hpp"

#include <cmath>
#include <expected>
#include <vector>

namespace lfmc {

EuropeanCall::EuropeanCall(double strike) : strike_(strike) {}

std::expected<std::vector<Payoffs>, std::string>
EuropeanCall::generate_payoffs(const std::vector<Path>& paths) const {
    Payoffs payoffs;
    for (const auto& path : paths) {
        double final_price = path.back();
        payoffs.push_back(std::max(final_price - strike_, 0.0));
    }
    return std::vector<Payoffs>{payoffs};
}

EuropeanPut::EuropeanPut(double strike) : strike_(strike) {}

std::expected<std::vector<Payoffs>, std::string>
EuropeanPut::generate_payoffs(const std::vector<Path>& paths) const {
    Payoffs payoffs;
    for (const auto& path : paths) {
        double final_price = path.back();
        payoffs.push_back(std::max(strike_ - final_price, 0.0));
    }
    return std::vector<Payoffs>{payoffs};
}

} // namespace lfmc
