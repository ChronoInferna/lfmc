#include "lfmc/estimator/monte_carlo_estimator.hpp"

#include <vector>

namespace lfmc {

std::expected<void, std::string>
MonteCarloEstimator::add_payoffs(const std::vector<Payoffs>& payoffs) {
    if (payoffs.empty()) {
        return std::unexpected("No payoffs provided to add to the estimator.");
    }

    for (const auto& payoff : payoffs) {
        sum += payoff[0];
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
    if (!converged()) {
        return std::unexpected("Estimator has not yet converged.");
    }

    return {sum / static_cast<double>(count)};
}

std::expected<void, std::string> MonteCarloEstimator::merge(Estimator const& other) {
    const auto* other_estimator = dynamic_cast<const MonteCarloEstimator*>(&other);
    if (!other_estimator) {
        return std::unexpected("Incompatible estimator type for merging.");
    }

    sum += other_estimator->sum;
    count += other_estimator->count;

    return {};
}

} // namespace lfmc
