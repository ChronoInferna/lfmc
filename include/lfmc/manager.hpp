#pragma once

#include "lfmc/process.hpp"
#include "lfmc/scheme.hpp"

/**
 * @file manager.hpp
 * @brief Manager class for handling processes and schemes using the Strategy Pattern.
 *
 * This file defines the `Manager` class, which is responsible for managing
 * multiple processes and schemes within the LFMC framework using the Strategy Pattern.
 *
 * @see Process
 * @see Scheme
 */

namespace lfmc {

/**
 * @class Manager
 * @brief Manages processes and schemes using the Strategy Pattern.
 *
 * The `Manager` class is responsible for handling multiple processes and schemes,
 * allowing for dynamic selection and management of different strategies at runtime.
 */
class Manager {
  public:
    /// Constructor
    Manager() = default;
    // Additional methods for managing processes and schemes will be added here

  private:
    // Internal data structures for managing processes and schemes will be added here
    // TODO decide on the appropriate strategy implementation - refer to notes and design patterns -
    // runtime vs. compile-time determines more
};

} // namespace lfmc
