#include "lfmc/payoff/barrier_payoffs.hpp"

#include "lfmc/core/types.hpp"

#include <algorithm>
#include <expected>
#include <string>
#include <vector>

namespace lfmc {

UpAndOutCall::UpAndOutCall(double strike, double barrier) : strike_(strike), barrier_(barrier) {}

std::expected<std::vector<Payoffs>, std::string>
UpAndOutCall::generate_payoffs(const std::vector<Path>& paths) const {
    Payoffs payoffs;
    payoffs.reserve(paths.size());

    for (const auto& path : paths) {
        if (path.empty())
            return std::unexpected("Empty path encountered in UpAndOutCall");

        bool knocked_out =
            std::any_of(path.begin(), path.end(), [this](double s) { return s >= barrier_; });

        payoffs.push_back(knocked_out ? 0.0 : std::max(path.back() - strike_, 0.0));
    }

    return std::vector<Payoffs>{payoffs};
}

DownAndInPut::DownAndInPut(double strike, double barrier) : strike_(strike), barrier_(barrier) {}

std::expected<std::vector<Payoffs>, std::string>
DownAndInPut::generate_payoffs(const std::vector<Path>& paths) const {
    Payoffs payoffs;
    payoffs.reserve(paths.size());

    for (const auto& path : paths) {
        if (path.empty())
            return std::unexpected("Empty path encountered in DownAndInPut");

        bool knocked_in =
            std::any_of(path.begin(), path.end(), [this](double s) { return s <= barrier_; });

        payoffs.push_back(knocked_in ? std::max(strike_ - path.back(), 0.0) : 0.0);
    }

    return std::vector<Payoffs>{payoffs};
}

} // namespace lfmc
