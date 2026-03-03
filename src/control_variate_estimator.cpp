#include "lfmc/estimator.hpp"

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

} // namespace lfmc
