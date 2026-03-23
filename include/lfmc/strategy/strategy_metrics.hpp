#pragma once

#include "lfmc/strategy/types.hpp"

#include <compare> // For operator<=>
#include <cstddef>

namespace lfmc {

struct StrategyMetrics {
    Strategy strategy;
    double mean;
    double variance;
    double std_error;
    size_t samples;
    long long elapsed_ms;

    // Define ordering based on variance
    auto operator<=>(const StrategyMetrics& other) const noexcept {
        return variance <=> other.variance;
    }
};

} // namespace lfmc
