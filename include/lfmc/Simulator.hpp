#pragma once
#include "NumericalScheme.hpp"
#include "Payoff.hpp"
#include "RandomGenerator.hpp"
#include "StochasticProcess.hpp"
#include "VarianceReductionStrategy.hpp"

#include <numeric>
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
    Simulator(P process, S scheme, PO payoff, double x0, double T, size_t n_steps,
              size_t n_simulations)
        : process_(std::move(process)), scheme_(std::move(scheme)), payoff_(std::move(payoff)),
          x0_(x0), T_(T), n_steps_(n_steps), n_simulations_(n_simulations) {}
    // TODO initialize other members

    // Used Claude to develope this sim. didn't feel like putting in the work since it will be
    // replaced later...
    double run_simulations() {
        double sum = 0.0;

        for (size_t i = 0; i < n_simulations_; ++i) {
            // Generate random normals for this path
            auto randoms = rng_.generate_normals(n_steps_);

            // Simulate terminal value
            double terminal = scheme_.simulate_terminal(x0_, T_, n_steps_, randoms);

            // Evaluate payoff
            double payoff_value = payoff_(terminal);

            sum += payoff_value;
        }

        // Return average (Monte Carlo estimate)
        return sum / n_simulations_;
    }

    /**
     * @brief Run simulations and return both mean and standard error.
     * @return Pair of (mean, standard_error).
     */
    std::pair<double, double> run_simulations_with_error() {
        std::vector<double> payoffs(n_simulations_);

        for (size_t i = 0; i < n_simulations_; ++i) {
            auto randoms = rng_.generate_normals(n_steps_);
            double terminal = scheme_.simulate_terminal(x0_, T_, n_steps_, randoms);
            payoffs[i] = payoff_(terminal);
        }

        // Calculate mean
        double mean = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / n_simulations_;

        // Calculate standard deviation
        double variance = 0.0;
        for (double payoff : payoffs) {
            double diff = payoff - mean;
            variance += diff * diff;
        }
        variance /= (n_simulations_ - 1);

        // Standard error = std_dev / sqrt(n)
        double std_error = std::sqrt(variance / n_simulations_);

        return {mean, std_error};
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
    double x0_;            // Initial value
    double T_;             // Time to maturity
    size_t n_steps_;       // Number of time steps
    size_t n_simulations_; // Number of Monte Carlo paths
    RandomGenerator rng_;
};

} // namespace lfmc
