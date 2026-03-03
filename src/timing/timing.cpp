#include "lfmc/timing.hpp"

namespace lfmc {

void Timer::reset() noexcept {
    start_time_ = clock::now();
}

long long Timer::elapsedMilliseconds() const noexcept {
    return std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - start_time_)
        .count();
}

ScopedTimer::~ScopedTimer() noexcept {
    *out_ =
        std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - start_time_).count();
}

} // namespace lfmc
