#include "lfmc/estimator.hpp"
#include "lfmc/types.hpp"

namespace lfmc {

std::expected<void, std::string> MonteCarloEstimator::add_payoffs(const Payoffs& payoffs) {
    if (payoffs.empty()) {
        return std::unexpected("No payoffs provided to add to the estimator.");
    }

    for (double payoff : payoffs) {
        sum += payoff;
        ++count;
    }

    return {};
}

bool MonteCarloEstimator::converged() const {
    return count >= 10000; // Placeholder convergence criterion
}

std::expected<double, std::string> MonteCarloEstimator::result() const {
    if (count <= 0) {
        return std::unexpected("No samples added to the estimator.");
    }
    // if (!self.converged()) {
    //     return std::unexpected("Estimator has not yet converged.");
    // }

    return {sum / static_cast<double>(count)};
}

// void MonteCarloEstimator::merge(Estimator const& other) {
//     auto const& mcOther = dynamic_cast<MonteCarloEstimator const&>(other);
//     sum += mcOther.sum;
//     count += mcOther.count;
// }

} // namespace lfmc
