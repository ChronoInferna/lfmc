#include "lfmc/timer.hpp"

namespace lfmc {

void Timer::reset() {
    start_time_ = clock::now();
}

long long Timer::elapsedMilliseconds() const {
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - start_time_).count();
    return duration;
}

ScopedTimer::~ScopedTimer() noexcept {
    *out_ =
        std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - start_time_).count();
}

} // namespace lfmc
