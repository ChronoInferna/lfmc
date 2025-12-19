#include "lfmc/timer.hpp"

namespace lfmc {

void Timer::reset() {
    start_time_ = clock::now();
}

long long Timer::elapsedMilliseconds() const {
    return std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - start_time_)
        .count();
}

ScopedTimer::~ScopedTimer() noexcept {
    *out_ =
        std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - start_time_).count();
}

} // namespace lfmc
