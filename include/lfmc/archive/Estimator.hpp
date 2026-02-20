#pragma once

#include "NumericalScheme.hpp"
#include "PathGenerator.hpp"
#include "RandomGenerator.hpp"
#include "StochasticProcess.hpp"
#include "types.hpp"

// #include <thread>
// #include <vector>

// #include <barrier>
// #include <memory>

namespace lfmc {

class EstimatorInterface {
  public:
    virtual ~EstimatorInterface() = default;
    virtual std::vector<Path> sample() = 0;

    // Expose for decorators to use - maybe not the best design, but allows for more flexible
    // variance reduction techniques
    virtual State const& getState() const = 0;
    virtual Normals generateNormals(size_t n) = 0;
    virtual Path generatePath(std::span<const double> randomNormals) = 0;
};

template <StochasticProcess P = GeometricBrownianMotion, NumericalScheme<P> S = EulerMaruyama<P>,
          RandomGenerator RNG = PseudoRandom>
class Estimator : public EstimatorInterface {
  private:
    P process_;
    S scheme_;
    RNG randomGenerator_;
    State state_;
    PathGenerator<P, S> pathGenerator_;

  public:
    explicit Estimator(P process, S scheme, State state, RNG randomGenerator)
        : process_(std::move(process)), scheme_(std::move(scheme)), state_(state),
          pathGenerator_(process_, scheme_, state_), randomGenerator_(std::move(randomGenerator)) {}
    explicit Estimator(P process, S scheme, State state)
        : process_(std::move(process)), scheme_(std::move(scheme)), state_(state),
          randomGenerator_(), pathGenerator_(process_, scheme_, state_) {}

    std::vector<Path> sample() override {
        Normals normals = randomGenerator_.generate(state_.steps);
        Path path = pathGenerator_.generate(normals);
        return {path};
    }

    State const& getState() const override {
        return state_;
    }

    Normals generateNormals(size_t n) override {
        return randomGenerator_.generate(n);
    }

    Path generatePath(std::span<const double> randomNormals) override {
        return pathGenerator_.generate(randomNormals);
    }

    // TODO change return type?
    // Restart thread?
    // this->restartThread();

    // TODO update atomic? some way for the test threads to update some shared state to indicate
    // which is the best

    // void restartThread();
};

} // namespace lfmc
