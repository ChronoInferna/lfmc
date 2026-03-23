#pragma once

#include "lfmc/strategy/types.hpp"

#include <compare> // For operator<=>
#include <cstddef>
#include <ostream>

namespace lfmc {

struct StrategyMetrics {
    Strategy strategy;
    double mean;
    double variance;
    double std_error;
    std::size_t samples;
    long long elapsed_ms;

    // Define ordering based on variance
    auto operator<=>(const StrategyMetrics& other) const noexcept {
        return variance <=> other.variance;
    }
};

inline std::ostream& operator<<(std::ostream& os, const Strategy& s) {
    return os << to_string(s); // or custom formatting
}

} // namespace lfmc
