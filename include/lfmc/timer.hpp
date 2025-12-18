#pragma once

#include <chrono>

namespace lfmc {

class Timer {
  private:
    std::chrono::high_resolution_clock::time_point start_time_;

  public:
    Timer() : start_time_(std::chrono::high_resolution_clock::now()) {}
    ~Timer() = default;

    void reset();
    long long elapsedMilliseconds() const;
};

// RAII timer that measures the time duration of a scope
// TODO unsure of how we want to implement this
class ScopedTimer {
  private:
    std::chrono::high_resolution_clock::time_point start_time_;

  public:
    ScopedTimer() : start_time_(std::chrono::high_resolution_clock::now()) {}
    // TODO use destructor? Then how to return the time?
    ~ScopedTimer();
};

} // namespace lfmc
