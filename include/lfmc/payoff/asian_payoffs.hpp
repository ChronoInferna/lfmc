#pragma once

#include "lfmc/core/types.hpp"
#include "lfmc/payoff/payoff.hpp"

#include <expected>
#include <vector>

namespace lfmc {

class AsianCall : public Payoff {
  private:
    double strike_;

  public:
    explicit AsianCall(double strike);

    std::expected<std::vector<Payoffs>, std::string>
    generate_payoffs(const std::vector<Path>& paths) const override;
};

class AsianPut : public Payoff {
  private:
    double strike_;

  public:
    explicit AsianPut(double strike);

    std::expected<std::vector<Payoffs>, std::string>
    generate_payoffs(const std::vector<Path>& paths) const override;
};

} // namespace lfmc
