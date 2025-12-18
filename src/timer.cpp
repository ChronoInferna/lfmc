#include "lfmc/timer.hpp"

namespace lfmc {

void Timer::reset() {
    start_time_ = std::chrono::high_resolution_clock::now();
}

long long Timer::elapsedMilliseconds() const {
    auto current_time = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time_).count();
    return duration;
}

ScopedTimer::~ScopedTimer() {
    // auto end_time = std::chrono::high_resolution_clock::now();
    // auto duration =
    //     std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time_).count();
}

} // namespace lfmc
