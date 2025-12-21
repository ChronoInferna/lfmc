#pragma once

/** @file process.hpp
 * @brief Process class representing a stochastic process.
 *
 * This file defines the `Process` class and implementations of different stochastic processes.
 */

namespace lfmc {

class Process {
  public:
    virtual ~Process() = default;
    // Define common interface methods for processes here
    // e.g., virtual void simulateStep() = 0;
    // TODO flesh out based on requirements
};

} // namespace lfmc
