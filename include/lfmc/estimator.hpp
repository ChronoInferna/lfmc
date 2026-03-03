#pragma once

#include "lfmc/types.hpp"

#include <expected>
#include <string>
#include <vector>

namespace lfmc {

class Estimator {
  public:
    virtual std::expected<void, std::string> add_payoffs(const std::vector<Payoffs>& payoffs) = 0;
    virtual bool converged() const = 0;
    virtual std::expected<double, std::string> result() const = 0;
    // virtual void merge(Estimator const& other) = 0;
    virtual ~Estimator() = default;
};

class MonteCarloEstimator : public Estimator {
  private:
    double sum = 0.0;
    std::size_t count = 0;

  public:
    std::expected<void, std::string> add_payoffs(const std::vector<Payoffs>& payoffs) override;
    bool converged() const override;
    std::expected<double, std::string> result() const override;
    // void merge(Estimator const& other) override;
};

// TODO if control is analytically known, can put it at payoff level
class ControlVariateEstimator : public Estimator {
  private:
    std::size_t count = 0;

    double sum_x = 0.0;  // Sum of original payoffs
    double sum_y = 0.0;  // Sum of control variate values
    double sum_xy = 0.0; // Sum of products of payoffs and control variate
    double sum_yy = 0.0; // Sum of squares of control variate values

    double control_expectation_ = 0.0; // Expected value of control variate (known analytically)

  public:
    explicit ControlVariateEstimator(double control_expectation)
        : control_expectation_(control_expectation) {}

    std::expected<void, std::string> add_payoffs(const std::vector<Payoffs>& payoffs) override;
    bool converged() const override;
    std::expected<double, std::string> result() const override;
    // void merge(Estimator const& other) override;
};

} // namespace lfmc
