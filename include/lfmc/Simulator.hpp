#pragma once

#include "NumericalScheme.hpp"
#include "Payoff.hpp"
#include "RandomGenerator.hpp"
#include "StochasticProcess.hpp"
// #include "VarianceReductionStrategy.hpp"

// #include <thread>
#include <vector>

// TODO
// #include <barrier>
// #include <memory>

namespace lfmc {

template <StochasticProcess P, NumericalScheme<P> S, Payoff PO> class Simulator {
  public:
    // NOTE VRS is const pointer since we do not need ownership - alternatively use shared_ptr but
    // probably overkill? More importantly, will threads outlive the manager that owns the strategy?
    Simulator(P process, S scheme, PO payoff, double initialValue, double timeToMaturity,
              size_t stepCount)
        : process_(std::move(process)), scheme_(std::move(scheme)), payoff_(std::move(payoff)),
          initialValue_(initialValue), timeToMaturity_(timeToMaturity), stepCount_(stepCount),
          rng_() {}

    std::vector<double> generatePath() {
        double dt = timeToMaturity_ / static_cast<float>(stepCount_);
        std::vector<double> path;
        path.reserve(stepCount_ + 1);
        path.push_back(initialValue_);
        std::vector<double> randomNormals = rng_.generateNormals(stepCount_);

        double x = initialValue_;
        for (size_t i{}; i < stepCount_; ++i) {
            x = scheme_.step(process_, x, dt, randomNormals[i]);
            path.push_back(x);
        }
        return path;
    }

    double generateTerminal() {
        double dt = timeToMaturity_ / static_cast<float>(stepCount_);
        std::vector<double> randomNormals = rng_.generateNormals(stepCount_);

        double x = initialValue_;
        for (size_t i{}; i < stepCount_; ++i) {
            x = scheme_.step(process_, x, dt, randomNormals[i]);
        }
        return x;
    }

    double generateTerminalWithRandoms(const std::vector<double>& randomNormals) {
        double dt = timeToMaturity_ / static_cast<float>(stepCount_);
        double x = initialValue_;

        for (size_t i{}; i < stepCount_; ++i) {
            x = scheme_.step(process_, x, dt, randomNormals[i]);
        }
        return x;
    }

    RandomGenerator& getRng() {
        return rng_;
    }

    // TODO change return type?
    // void setStrategy() noexcept {
    // Restart thread?
    // this->restartThread();

    // TODO update atomic? some way for the test threads to update some shared state to indicate
    // which is the best

    // void restartThread();

  private:
    P process_;
    S scheme_;
    PO payoff_;
    RandomGenerator rng_;
    double initialValue_;   // Initial value
    double timeToMaturity_; // Time to maturity
    size_t stepCount_;      // Number of time steps
};

} // namespace lfmc
