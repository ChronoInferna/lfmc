#pragma once

#include "lfmc/core/types.hpp"

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

} // namespace lfmc
