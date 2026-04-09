#pragma once

#include "lfmc/payoff/payoff.hpp"

#include <expected>
#include <memory>
#include <vector>

namespace lfmc {

// TODO subpar in terms of architecture because we have both control variate payoff and estimator,
// but for now it's simpler to implement it this way - can enforce with a factory. Can refactor
// later if needed.
class ControlVariatePayoff : public Payoff {
  private:
    std::unique_ptr<Payoff> target_payoff_;
    std::unique_ptr<Payoff> control_payoff_;

  public:
    explicit ControlVariatePayoff(std::unique_ptr<Payoff> target_payoff,
                                  std::unique_ptr<Payoff> control_payoff);

    std::expected<std::vector<Payoffs>, std::string>
    generate_payoffs(const std::vector<Path>& paths) const override;
};

} // namespace lfmc
