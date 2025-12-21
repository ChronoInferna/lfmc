#pragma once
#include "lfmc/process.hpp"
#include "lfmc/scheme.hpp"

#include <memory>

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
    // TODO as a starting point, we use runtime polymorphism to manage strategies
  public:
    explicit Manager(std::unique_ptr<Process>&& process = {}, std::unique_ptr<Scheme>&& scheme = {})
        : scheme_(std::move(scheme)), process_(std::move(process)) {}

    // TODO move to implementation file once actual design is fleshed out
    void setScheme(std::unique_ptr<Scheme>&& scheme) {
        scheme_ = std::move(scheme);
    }
    void setProcess(std::unique_ptr<Process>&& process) {
        process_ = std::move(process);
    }

  private:
    // TODO Internal data structures for managing processes and schemes will be added here
    // TODO decide on the appropriate strategy implementation - refer to notes and design patterns -
    // runtime vs. compile-time

    std::unique_ptr<Scheme> scheme_;
    std::unique_ptr<Process> process_;
};

} // namespace lfmc
