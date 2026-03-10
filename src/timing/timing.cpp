#include "lfmc/timing/timing.hpp"

#include <chrono>

namespace lfmc {

Timer::Timer() noexcept : start_time_(clock::now()) {}

void Timer::reset() noexcept {
    start_time_ = clock::now();
}

long long Timer::elapsedMilliseconds() const noexcept {
    return std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - start_time_)
        .count();
}

ScopedTimer::ScopedTimer(long long& out) noexcept : out_(&out), start_time_(clock::now()) {}

ScopedTimer::~ScopedTimer() noexcept {
    *out_ =
        std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - start_time_).count();
}

} // namespace lfmc
