#pragma once

#include <cstddef>
#include <string>

namespace lfmc {

struct StrategyMetrics {
    std::string name;
    double mean;
    double variance;
    double std_error;
    std::size_t samples;
    long long elapsed_ms;

    bool operator<(const StrategyMetrics& other) const noexcept {
        return variance < other.variance;
    }
};

} // namespace lfmc