#include "lfmc/payoff.hpp"
#include "lfmc/types.hpp"

#include <algorithm>
#include <expected>
#include <numeric>
#include <string>
#include <vector>

namespace lfmc {

AsianCall::AsianCall(double strike) : strike_(strike) {}

std::expected<std::vector<Payoffs>, std::string>
AsianCall::generate_payoffs(const std::vector<Path>& paths) const {
    Payoffs payoffs;
    payoffs.reserve(paths.size());

    for (const auto& path : paths) {
        if (path.empty())
            return std::unexpected("Empty path encountered in AsianCall");

        double mean = std::reduce(path.begin(), path.end(), 0.0) / static_cast<double>(path.size());
        payoffs.push_back(std::max(mean - strike_, 0.0));
    }

    return std::vector<Payoffs>{payoffs};
}

AsianPut::AsianPut(double strike) : strike_(strike) {}

std::expected<std::vector<Payoffs>, std::string>
AsianPut::generate_payoffs(const std::vector<Path>& paths) const {
    Payoffs payoffs;
    payoffs.reserve(paths.size());

    for (const auto& path : paths) {
        if (path.empty())
            return std::unexpected("Empty path encountered in AsianPut");

        double mean = std::reduce(path.begin(), path.end(), 0.0) / static_cast<double>(path.size());
        payoffs.push_back(std::max(strike_ - mean, 0.0));
    }

    return std::vector<Payoffs>{payoffs};
}

} // namespace lfmc
