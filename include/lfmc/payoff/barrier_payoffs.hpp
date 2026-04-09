#pragma once

#include "lfmc/core/types.hpp"
#include "lfmc/payoff/payoff.hpp"

#include <expected>
#include <vector>

namespace lfmc {

class UpAndOutCall : public Payoff {
  private:
    double strike_;
    double barrier_;

  public:
    UpAndOutCall(double strike, double barrier);

    std::expected<std::vector<Payoffs>, std::string>
    generate_payoffs(const std::vector<Path>& paths) const override;
};

class DownAndInPut : public Payoff {
  private:
    double strike_;
    double barrier_;

  public:
    DownAndInPut(double strike, double barrier);

    std::expected<std::vector<Payoffs>, std::string>
    generate_payoffs(const std::vector<Path>& paths) const override;
};

} // namespace lfmc
