#include "lfmc/Payoff.hpp"

namespace lfmc {

struct EuropeanCall {
    double strike;

    double operator()(double terminal) const {
        return std::max(terminal - strike, 0.0);
    }
};

struct EuropeanPut {
    double strike;

    double operator()(double terminal) const {
        return std::max(strike - terminal, 0.0);
    }
};

} // namespace lfmc
