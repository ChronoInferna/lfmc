#include "lfmc/payoff/lookback_payoffs.hpp"

#include "lfmc/core/types.hpp"

#include <algorithm>
#include <expected>
#include <string>
#include <vector>

namespace lfmc {

std::expected<std::vector<Payoffs>, std::string>
LookbackCall::generate_payoffs(const std::vector<Path>& paths) const {
    Payoffs payoffs;
    payoffs.reserve(paths.size());

    for (const auto& path : paths) {
        if (path.empty())
            return std::unexpected("Empty path encountered in LookbackCall");

        double min_price = *std::min_element(path.begin(), path.end());
        payoffs.push_back(path.back() - min_price);
    }

    return std::vector<Payoffs>{payoffs};
}

std::expected<std::vector<Payoffs>, std::string>
LookbackPut::generate_payoffs(const std::vector<Path>& paths) const {
    Payoffs payoffs;
    payoffs.reserve(paths.size());

    for (const auto& path : paths) {
        if (path.empty())
            return std::unexpected("Empty path encountered in LookbackPut");

        double max_price = *std::max_element(path.begin(), path.end());
        payoffs.push_back(max_price - path.back());
    }

    return std::vector<Payoffs>{payoffs};
}

} // namespace lfmc
