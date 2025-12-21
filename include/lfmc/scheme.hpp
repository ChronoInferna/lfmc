#pragma once

/* @file scheme.hpp
 * @brief Scheme class representing numerical schemes for stochastic processes.
 *
 * This file defines the `Scheme` class and implementations of different numerical schemes
 * used for simulating stochastic processes.
 */

namespace lfmc {

class Scheme {
  public:
    virtual ~Scheme() = default;
    // Define common interface methods for schemes here
};

} // namespace lfmc
