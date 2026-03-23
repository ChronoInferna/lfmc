#pragma once

#include "lfmc/core/types.hpp"
#include "lfmc/estimator/estimator.hpp"

#include <cmath>
#include <expected>
#include <string>
#include <vector>

namespace lfmc {

// TODO if control is analytically known, can put it at payoff level
class ControlVariateEstimator : public Estimator {
  private:
    std::size_t count = 0;

    double sum_x = 0.0;  // Sum of original payoffs
    double sum_y = 0.0;  // Sum of control variate values
    double sum_xy = 0.0; // Sum of products of payoffs and control variate
    double sum_yy = 0.0; // Sum of squares of control variate values
    double sum_xx = 0.0;

    double control_expectation_ = 0.0; // Expected value of control variate (known analytically)

  public:
    explicit ControlVariateEstimator(double control_expectation)
        : control_expectation_(control_expectation) {}

    std::expected<void, std::string> add_payoffs(const std::vector<Payoffs>& payoffs) override;
    bool converged() const override;
    std::expected<double, std::string> result() const override;
    std::expected<void, std::string> merge(Estimator const& other) override;

    double mean() const noexcept;
    double variance() const noexcept;
    double std_error() const noexcept;
    std::size_t sample_count() const noexcept;
};

} // namespace lfmc
