#pragma once

#include "lfmc/types.hpp"

namespace lfmc {

class Payoff {
  public:
    virtual Payoffs generate_payoffs(const std::vector<Path>& paths) const = 0;
    virtual ~Payoff() = default;
};

class EuropeanCall : public Payoff {
  private:
    double strike_;

  public:
    explicit EuropeanCall(double strike);

    Payoffs generate_payoffs(const std::vector<Path>& paths) const override;
};

class EuropeanPut : public Payoff {
  private:
    double strike_;

  public:
    explicit EuropeanPut(double strike);

    Payoffs generate_payoffs(const std::vector<Path>& paths) const override;
};

} // namespace lfmc
