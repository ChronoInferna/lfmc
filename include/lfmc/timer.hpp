#pragma once

#include <chrono>

namespace lfmc {

class Timer {
    using clock = std::chrono::high_resolution_clock;

  private:
    clock::time_point start_time_;

  public:
    Timer() : start_time_(clock::now()) {}
    ~Timer() = default;

    void reset();
    long long elapsedMilliseconds() const;
};

// RAII timer that measures the time duration of a scope
class ScopedTimer {
    using clock = std::chrono::high_resolution_clock;

  private:
    long long* out_;
    clock::time_point start_time_;

  public:
    explicit ScopedTimer(long long& out) noexcept : out_(&out), start_time_(clock::now()) {}
    ~ScopedTimer() noexcept;
};

} // namespace lfmc
