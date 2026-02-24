#pragma once

#include <expected>
#include <string>

namespace lfmc {

class Estimator {
  public:
    virtual void add(double x) = 0;
    // virtual bool converged() const = 0;
    virtual std::expected<double, std::string> result() const = 0;
    // virtual void merge(Estimator const& other) = 0;
    virtual ~Estimator() = default;
};

} // namespace lfmc
