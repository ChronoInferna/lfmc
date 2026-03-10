#pragma once

#include "lfmc/payoff/payoff.hpp"

namespace lfmc {

class LookbackCall : public Payoff {
  public:
    std::expected<std::vector<Payoffs>, std::string>

    generate_payoffs(const std::vector<Path>& paths) const override;
};

class LookbackPut : public Payoff {
  public:
    std::expected<std::vector<Payoffs>, std::string>
    generate_payoffs(const std::vector<Path>& paths) const override;
};

} // namespace lfmc
