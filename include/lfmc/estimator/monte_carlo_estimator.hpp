#pragma once

#include "lfmc/estimator/estimator.hpp"

#include <expected>
#include <string>
#include <vector>

namespace lfmc {

class MonteCarloEstimator : public Estimator {
  private:
    double sum = 0.0;
    std::size_t count = 0;

  public:
    std::expected<void, std::string> add_payoffs(const std::vector<Payoffs>& payoffs) override;
    bool converged() const override;
    std::expected<double, std::string> result() const override;
};

} // namespace lfmc
