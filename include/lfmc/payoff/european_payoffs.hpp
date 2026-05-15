#pragma once

#include "lfmc/core/types.hpp"
#include "lfmc/payoff/payoff.hpp"

#include <expected>
#include <vector>

namespace lfmc {

class EuropeanCall : public Payoff {
  private:
    double strike_;

  public:
    explicit EuropeanCall(double strike);

    std::expected<std::vector<Payoffs>, std::string>
    generate_payoffs(const std::vector<Path>& paths) const override;
};

class EuropeanPut : public Payoff {
  private:
    double strike_;

  public:
    explicit EuropeanPut(double strike);

    std::expected<std::vector<Payoffs>, std::string>
    generate_payoffs(const std::vector<Path>& paths) const override;
};

} // namespace lfmc
