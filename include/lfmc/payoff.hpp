#pragma once

#include "lfmc/types.hpp"

#include <expected>
#include <memory>
#include <vector>

namespace lfmc {

class Payoff {
  public:
    virtual std::expected<std::vector<Payoffs>, std::string>
    generate_payoffs(const std::vector<Path>& paths) const = 0;
    virtual ~Payoff() = default;
};

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
