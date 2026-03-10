#include "lfmc/estimator/monte_carlo_estimator.hpp"

#include <vector>
#include <cmath>    // ← ADD THIS for std::sqrt
#include <limits>   // ← ADD THIS for std::numeric_limits
namespace lfmc {

std::expected<void, std::string>
MonteCarloEstimator::add_payoffs(const std::vector<Payoffs>& payoffs) {
    if (payoffs.empty()) {
        return std::unexpected("No payoffs provided to add to the estimator.");
    }

    for (const auto& payoff : payoffs) {
        double value = payoff[0];
        sum += value;
        sum_sq += value * value;
        ++count;
    }

    return {};
}


double MonteCarloEstimator::mean() const noexcept {
    return count > 0 ? sum / static_cast<double>(count) : 0.0;
}

double MonteCarloEstimator::variance() const noexcept {
    if (count < 2) return std::numeric_limits<double>::max();
    double m = mean();
    return (sum_sq / static_cast<double>(count)) - (m * m);  // ← Add explicit cast
}

double MonteCarloEstimator::std_error() const noexcept {
    if (count < 2) return std::numeric_limits<double>::max();
    return std::sqrt(variance() / static_cast<double>(count));
}

std::size_t MonteCarloEstimator::sample_count() const noexcept {
    return count;
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
    sum_sq += other_estimator->sum_sq;
    count += other_estimator->count;

    return {};
}

} // namespace lfmc
