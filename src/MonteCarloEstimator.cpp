#include "lfmc/Estimator.hpp"

#include <cstddef>

namespace lfmc {

class MonteCarloEstimator : public Estimator {
  private:
    double sum = 0.0;
    std::size_t count = 0;

  public:
    void add(double x) override {
        sum += x;
        ++count;
    }

    // bool converged() const override {
    //     return count >= 10000; // Placeholder convergence criterion
    // }

    std::expected<double, std::string> result() const override {
        if (count == 0) {
            return std::unexpected{"No samples added"};
        }
        return {sum / static_cast<double>(count)};
    }

    // void merge(Estimator const& other) override {
    //     auto const& mcOther = dynamic_cast<MonteCarloEstimator const&>(other);
    //     sum += mcOther.sum;
    //     count += mcOther.count;
    // }
};

} // namespace lfmc
