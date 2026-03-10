#include "lfmc/estimator/control_variate_estimator.hpp"

namespace lfmc {

std::expected<void, std::string>
ControlVariateEstimator::add_payoffs(const std::vector<Payoffs>& payoffs) {
    if (payoffs.empty()) {
        return std::unexpected("No payoffs provided to add to the estimator.");
    }

    for (const auto& row : payoffs) {
        double x = row[0];
        double y = row[1];

        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_yy += y * y;
        sum_xx += x * x;
        ++count;
    }

    return {};
}

bool ControlVariateEstimator::converged() const {
    return count >= 10000; // Placeholder convergence criterion
}

std::expected<double, std::string> ControlVariateEstimator::result() const {
    if (count <= 0) {
        return std::unexpected("No samples added to the estimator.");
    }
    if (!converged()) {
        return std::unexpected("Estimator has not yet converged.");
    }

    double n = static_cast<double>(count);

    double mean_x = sum_x / n;
    double mean_y = sum_y / n;

    double cov_xy = (sum_xy / n) - mean_x * mean_y;
    double var_y = (sum_yy / n) - mean_y * mean_y;

    if (var_y == 0.0)
        return std::unexpected("Zero variance in control variate");

    double beta = cov_xy / var_y;

    return mean_x - beta * (mean_y - control_expectation_);
}

std::expected<void, std::string> ControlVariateEstimator::merge(Estimator const& other) {
    const auto* other_estimator = dynamic_cast<const ControlVariateEstimator*>(&other);
    if (!other_estimator) {
        return std::unexpected("Incompatible estimator type for merging");
    }

    count += other_estimator->count;
    sum_x += other_estimator->sum_x;
    sum_y += other_estimator->sum_y;
    sum_xy += other_estimator->sum_xy;
    sum_yy += other_estimator->sum_yy;
    sum_xx += other_estimator->sum_xx;

    return {};
}

double ControlVariateEstimator::mean() const noexcept {
    if (count == 0) return 0.0;
    
    double n = static_cast<double>(count);
    double mean_x = sum_x / n;
    double mean_y = sum_y / n;
    
    double cov_xy = (sum_xy / n) - mean_x * mean_y;
    double var_y = (sum_yy / n) - mean_y * mean_y;
    
    if (var_y == 0.0) return mean_x;
    
    double beta = cov_xy / var_y;
    return mean_x - beta * (mean_y - control_expectation_);
}

double ControlVariateEstimator::variance() const noexcept {
    if (count < 2) return std::numeric_limits<double>::max();
    
    double n = static_cast<double>(count);
    double mean_x = sum_x / n;
    double mean_y = sum_y / n;
    
    double var_x = (sum_xx / n) - mean_x * mean_x;
    double var_y = (sum_yy / n) - mean_y * mean_y;
    double cov_xy = (sum_xy / n) - mean_x * mean_y;
    
    if (var_y == 0.0) return var_x;
    
    double beta = cov_xy / var_y;
    
    // Var(X - β(Y - E[Y])) = Var(X) + β²Var(Y) - 2βCov(X,Y)
    return var_x + beta * beta * var_y - 2.0 * beta * cov_xy;
}

double ControlVariateEstimator::std_error() const noexcept {
    if (count < 2) return std::numeric_limits<double>::max();
    return std::sqrt(variance() / static_cast<double>(count));
}

std::size_t ControlVariateEstimator::sample_count() const noexcept {
    return count;
}

} // namespace lfmc
