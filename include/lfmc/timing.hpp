#pragma once
#include <chrono>

/**
 * @file timing.hpp
 * @brief High-resolution Timer and Scoped Timer classes for measuring elapsed time.
 *
 * This file defines two classes: `Timer` and `ScopedTimer`.
 * - `Timer`: A high-resolution timer that can be started, stopped, and reset to measure elapsed
 * time in milliseconds.
 * - `ScopedTimer`: A RAII-style timer that measures the duration of a scope and outputs the elapsed
 * time to a provided variable.
 */

namespace lfmc {

/**
 * @class Timer
 * @brief A high-resolution timer for measuring elapsed time in milliseconds.
 *
 * The `Timer` class provides functionality to start, stop, and reset a timer,
 * allowing users to measure elapsed time intervals with high precision.
 */
class Timer {
    using clock = std::chrono::high_resolution_clock;

  public:
    /// Constructor
    Timer() noexcept : start_time_(clock::now()) {}

    /// Resets the timer to the current time.
    void reset() noexcept;
    /** @brief Returns the elapsed time in milliseconds since the timer was started or last reset.
     *
     * @return Elapsed time in milliseconds.
     */
    long long elapsedMilliseconds() const noexcept;

  private:
    clock::time_point start_time_;
};

/**
 * @class ScopedTimer
 * @brief A RAII-style timer that measures the duration of a scope.
 *
 * The `ScopedTimer` class starts timing upon construction and stops timing
 * upon destruction, outputting the elapsed time in milliseconds to a provided variable.
 */
class ScopedTimer {
    using clock = std::chrono::high_resolution_clock;

  public:
    /** @brief Constructs a ScopedTimer that outputs elapsed time to the provided variable.
     *
     * @param out Reference to a long long variable where the elapsed time in milliseconds will be
     * stored upon destruction.
     */
    explicit ScopedTimer(long long& out) noexcept : out_(&out), start_time_(clock::now()) {}
    /** @brief Destructor calculates and stores the elapsed time in milliseconds.
     *
     * Upon destruction, the elapsed time since construction is calculated and stored
     * in the variable provided during construction.
     */
    ~ScopedTimer() noexcept;

  private:
    long long* out_;
    clock::time_point start_time_;
};

} // namespace lfmc
