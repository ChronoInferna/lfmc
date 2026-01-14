#pragma once
#include "NumericalScheme.hpp"
#include "StochasticProcess.hpp"
#include "VarianceReductionStrategy.hpp"

// #include <expected> TODO
#include <memory>
#include <thread>
#include <variant>
#include <vector>

/**
 * @file manager.hpp
 * @brief Manager class for variation reduction strategies using the Strategy Pattern.
 *
 * This file defines the `Manager` class, which is responsible for managing
 * different strategies and potentially coordinating between them during runtime.
 */

namespace lfmc {

/**
 * @brief Concept for variant types used in variance reduction strategies.
 *
 * A type `T` satisfies the `Variant` concept if it is a `std::variant`
 * with at least one alternative type.
 *
 * @tparam T The type to be checked against the concept.
 */
template <typename T>
concept Variant = requires { std::variant_size_v<T>; } && std::variant_size_v<T> > 0;

/**
 * @brief Manager class for handling variance reduction strategies.
 *
 * The `Manager` class utilizes the Strategy Pattern to manage different
 * variance reduction strategies during Monte Carlo simulations. It holds
 * instances of a stochastic process, a numerical scheme, and a selected
 * variance reduction strategy.
 *
 * @tparam P The type of the stochastic process, satisfying the `StochasticProcess` concept.
 * @tparam S The type of the numerical scheme, satisfying the `NumericalScheme` concept, which
 * involves having the correct process_type P.
 * @tparam VRVariant A variant type containing different variance reduction strategies,
 *                   satisfying the `Variant` concept.
 */
template <StochasticProcess P, NumericalScheme S, Variant VRVariant> class Manager {
  public:
    explicit Manager(P process, S scheme,
                     std::unique_ptr<VarianceReductionStrategy> strategy =
                         std::make_unique<lfmc::NoVarianceReduction>()) noexcept
        : process_(std::move(process)), scheme_(std::move(scheme)),
          currentStrategy_(std::move(strategy)) {}

    // TODO Randomness here? Or in GBM?

  private:
    P process_;
    S scheme_;

    std::unique_ptr<VarianceReductionStrategy> currentStrategy_;

    static constexpr std::size_t numStrategies = std::variant_size_v<VRVariant>;
    std::array<VRVariant, numStrategies> testingThreads;

    std::vector<std::thread> realThreads;
};

} // namespace lfmc
