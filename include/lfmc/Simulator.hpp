#include "NumericalScheme.hpp"
#include "StochasticProcess.hpp"
#include "VarianceReductionStrategy.hpp"

// TODO
// #include <barrier>
// #include <memory>
#include <thread>

namespace lfmc {

/* @brief Simulator worker class that runs Monte Carlo simulations using a specified stochastic
 * process, numerical scheme, and variance reduction strategy.
 *
 * This class is responsible for executing simulations in a separate thread, applying the chosen
 * variance reduction technique to improve the efficiency of the Monte Carlo method.
 *
 * @tparam P The type of the stochastic process, must satisfy the StochasticProcess concept.
 * @tparam S The type of the numerical scheme, must satisfy the NumericalScheme concept and be
 * associated with the stochastic process P.
 */
template <StochasticProcess P, NumericalScheme<P> S> class Simulator {
  public:
    // NOTE VRS is const pointer since we do not need ownership - alternatively use shared_ptr but
    // probably overkill? More importantly, will threads outlive the manager that owns the strategy?
    explicit Simulator(const P& process, const S& scheme,
                       const VarianceReductionStrategy& strategy) noexcept
        : process_(process), scheme_(scheme), thread_(&Simulator::simulate(), this) {}
    // TODO how to actually initialize thread?
    // TODO destructor

    // TODO change return type?
    void setStrategy() noexcept {
        // Restart thread?
        // this->restartThread();
    }

    // TODO update atomic? some way for the test threads to update some shared state to indicate
    // which is the best

    // void restartThread();

  private:
    std::jthread thread_;
    // std::barrier state_;

    // TODO function that actually simulates
    void simulate() {}

    P process_;
    S scheme_;
    // TODO
    // int window;
};

} // namespace lfmc
