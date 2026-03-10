#pragma once

#include "lfmc/core/types.hpp"
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
    std::expected<void, std::string> merge(Estimator const& other) override;
};

} // namespace lfmc
