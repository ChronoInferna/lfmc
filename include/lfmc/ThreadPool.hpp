#pragma once

#include <thread>

/**
 * @file ThreadPool.hpp
 * @brief Defines a ThreadPool for managing a pool of threads.
 *
 * This file provides an implementation of a thread pool specific to Monte Carlo simulations and
 * variance reduction techniques.
 */

namespace lfmc {

/**
 * @brief ThreadPool class for managing a pool of threads.
 *
 * The ThreadPool class provides a simple interface for creating and managing a pool of threads
 * that can be used to execute tasks concurrently. This is particularly useful in Monte Carlo
 * simulations where multiple independent simulations can be run in parallel.
 */
class ThreadPool {
  public:
    /**
     * @brief Construct a ThreadPool with the specified number of threads.
     * @param numThreads The number of threads in the pool.
     */
    explicit ThreadPool(std::size_t numThreads) noexcept : numThreads_(numThreads) {
        // Implementation for initializing the thread pool can be added here.
    }

    /**
     * @brief Get the number of threads in the pool.
     * @return The number of threads.
     */
    std::size_t getNumThreads() const noexcept {
        return numThreads_;
    }

    /** Additional methods for submitting tasks and managing threads can be added here. */

  private:
    std::size_t numThreads_;
};

} // namespace lfmc
