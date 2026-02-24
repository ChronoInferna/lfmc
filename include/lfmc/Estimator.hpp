#pragma once

namespace lfmc {

class Estimator {
  public:
    virtual void add(double x) = 0;
    // virtual bool converged() const = 0;
    // virtual Result result() const = 0;
    // virtual void merge(Estimator const& other) = 0;
    virtual ~Estimator() = default;
};

} // namespace lfmc
