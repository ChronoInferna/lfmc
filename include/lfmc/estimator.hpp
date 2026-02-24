#pragma once

#include "lfmc/types.hpp"

#include <expected>
#include <string>

namespace lfmc {

class Estimator {
  public:
    virtual std::expected<void, std::string> add_payoffs(const Payoffs& payoffs) = 0;
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
    std::expected<void, std::string> add_payoffs(const Payoffs& payoffs) override;
    bool converged() const override;
    std::expected<double, std::string> result() const override;
    // void merge(Estimator const& other) override;
};

} // namespace lfmc
