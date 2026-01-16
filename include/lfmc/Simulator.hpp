#include "NumericalScheme.hpp"
#include "StochasticProcess.hpp"
#include "VarianceReductionStrategy.hpp"

#include <memory>
#include <thread>

namespace lfmc {

template <StochasticProcess P, NumericalScheme S>
    requires std::same_as<typename S::process_type, P>
class Simulator {
  public:
    explicit Simulator(std::unique_ptr<VarianceReductionStrategy> strategy) noexcept
        : currentStrategy_(std::move(strategy)) {}

    // TODO change return type?
    void setStrategy(std::unique_ptr<VarianceReductionStrategy> strategy) {
        currentStrategy_ = std::move(strategy);
    }

  private:
    std::jthread thread;
    std::unique_ptr<VarianceReductionStrategy> currentStrategy_;
    int window;
};
} // namespace lfmc
