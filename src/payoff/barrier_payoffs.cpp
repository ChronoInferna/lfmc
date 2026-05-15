#include "lfmc/payoff/barrier_payoffs.hpp"

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

        if (knocked_out) {
            payoffs.push_back(0.0);
        } else {
            payoffs.push_back(std::max(path.back() - strike_, 0.0));
        }
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

        if (knocked_in) {
            payoffs.push_back(std::max(strike_ - path.back(), 0.0));

        } else {
            payoffs.push_back(0.0);
        }
    }

    return std::vector<Payoffs>{payoffs};
}

} // namespace lfmc
