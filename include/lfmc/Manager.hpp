#pragma once
#include "NumericalScheme.hpp"
#include "Simulator.hpp"
#include "StochasticProcess.hpp"
#include "VarianceReductionStrategy.hpp"

// TODO
// #include <expected>
#include <memory>
// TODO
// #include <thread>
#include <vector>

/**
 * @file manager.hpp
 * @brief Manager class for variation reduction strategies using the Strategy Pattern.
 *
 * This file defines the `Manager` class, which is responsible for managing
 * different strategies and potentially coordinating between them during runtime.
 */

namespace lfmc {

template <StochasticProcess P, NumericalScheme S, VRStrategy... VRStrategies>
    requires std::same_as<typename S::process_type, P> && (sizeof...(VRStrategies) > 0)
class Manager {
  public:
    explicit Manager(P process, S scheme,
                     std::unique_ptr<VarianceReductionStrategy> strategy) noexcept
        : process_(std::move(process)), scheme_(std::move(scheme)),
          currentStrategy_(std::move(strategy)) {}
    explicit Manager(P process, S scheme)
        : Manager(std::move(process), std::move(scheme),
                  std::make_unique<lfmc::NoVarianceReduction>()) {}

  private:
    // TODO do we need these?
    P process_;
    S scheme_;

    std::unique_ptr<VarianceReductionStrategy> currentStrategy_;

    // std::tuple<VRStrategies...> strategies_;
    std::array<lfmc::Simulator<P, S>, sizeof...(VRStrategies)> testingThreads;

    // std::vector<std::thread> realThreads;
    std::vector<lfmc::Simulator<P, S>> realThreads;

    // TODO Trying to think about how we represent each strategy within each thread and if managing
    // threads through a vector of objects is viable
};

} // namespace lfmc
