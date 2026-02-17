#pragma once

#include "NumericalScheme.hpp"
#include "Payoff.hpp"
#include "RandomGenerator.hpp"
#include "StochasticProcess.hpp"
// #include "VarianceReductionStrategy.hpp"

// #include <thread>
#include <vector>

// #include <barrier>
#include <memory>

namespace lfmc {

// TODO break this up into multiple engines for various parts like path generation, aggregation,
// etc. since otherwise this will cause a lot of blowup - i.e. every new variation technique
// requires overriding every single method here
class SimulatorInterface {
  public:
    virtual ~SimulatorInterface() = default;

    // This is ideally the only thing that the manager should call, and the rest of the methods need
    // to be moved and are for variance reduction techniques to use when they override this method
    // virtual double sample() = 0;

    virtual double generateTerminal() = 0;
    virtual double generatePayoff() = 0;
    virtual double generateTerminalFromRandoms(const std::vector<double>& randomNormals) = 0;
    virtual double generatePayoffFromRandoms(const std::vector<double>& randomNormals) = 0;
    // Unneeded as of now (since we only need the terminals for the estimator), but may be useful
    // later for some variance reduction techniques virtual std::vector<double> generatePath() = 0;

    virtual RandomGenerator& getRng() = 0;
    virtual size_t getStepCount() const noexcept = 0;
};

template <StochasticProcess P, NumericalScheme<P> S, Payoff PO>
class Simulator : public SimulatorInterface {
  private:
    P process_;
    S scheme_;
    PO payoff_;
    std::unique_ptr<RandomGenerator> rng_; // Random number generator for the simulation
    double initialValue_;                  // Initial value
    double timeToMaturity_;                // Time to maturity
    size_t stepCount_;                     // Number of time steps

  public:
    Simulator(P process, S scheme, PO payoff, double initialValue, double timeToMaturity,
              size_t stepCount)
        : process_(std::move(process)), scheme_(std::move(scheme)), payoff_(std::move(payoff)),
          initialValue_(initialValue), timeToMaturity_(timeToMaturity), stepCount_(stepCount),
          rng_(std::make_unique<RandomGenerator>()) {}

    // std::vector<double> generatePath() {
    //     double dt = timeToMaturity_ / static_cast<float>(stepCount_);
    //     std::vector<double> path;
    //     path.reserve(stepCount_ + 1);
    //     path.push_back(initialValue_);
    //     std::vector<double> randomNormals = rng_->generateNormals(stepCount_);
    //
    //     double x = initialValue_;
    //     for (size_t i{}; i < stepCount_; ++i) {
    //         x = scheme_.step(process_, x, dt, randomNormals[i]);
    //         path.push_back(x);
    //     }
    //     return path;
    // }

    double generateTerminal() override {
        double dt = timeToMaturity_ / static_cast<float>(stepCount_);
        std::vector<double> randomNormals = rng_->generateNormals(stepCount_);

        double x = initialValue_;
        for (size_t i{}; i < stepCount_; ++i) {
            x = scheme_.step(process_, x, dt, randomNormals[i]);
        }
        return x;
    }

    double generatePayoff() override {
        double terminal = this->generateTerminal();
        return payoff_(terminal);
    }

    double generateTerminalFromRandoms(const std::vector<double>& randomNormals) override {
        double dt = timeToMaturity_ / static_cast<float>(stepCount_);
        double x = initialValue_;

        for (size_t i{}; i < stepCount_; ++i) {
            x = scheme_.step(process_, x, dt, randomNormals[i]);
        }
        return x;
    }

    double generatePayoffFromRandoms(const std::vector<double>& randomNormals) override {
        double terminal = this->generateTerminalFromRandoms(randomNormals);
        return payoff_(terminal);
    }

    RandomGenerator& getRng() override {
        return *rng_;
    }

    size_t getStepCount() const noexcept override {
        return stepCount_;
    }

    // TODO change return type?
    // Restart thread?
    // this->restartThread();

    // TODO update atomic? some way for the test threads to update some shared state to indicate
    // which is the best

    // void restartThread();
};

} // namespace lfmc
