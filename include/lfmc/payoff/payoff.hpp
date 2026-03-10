#pragma once

#include "lfmc/core/types.hpp"

#include <expected>
#include <string>
#include <vector>

namespace lfmc {

class Payoff {
  public:
    virtual std::expected<std::vector<Payoffs>, std::string>
    generate_payoffs(const std::vector<Path>& paths) const = 0;
    virtual ~Payoff() = default;
};

} // namespace lfmc
