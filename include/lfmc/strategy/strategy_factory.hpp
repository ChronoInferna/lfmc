#pragma once

#include "lfmc/estimator/control_variate_estimator.hpp"
#include "lfmc/estimator/monte_carlo_estimator.hpp"
#include "lfmc/numerical_scheme/euler_maruyama.hpp"
#include "lfmc/path_generator/path_generator.hpp"
#include "lfmc/payoff/control_variate_payoffs.hpp"
#include "lfmc/payoff/payoff.hpp"
#include "lfmc/pipeline/pipeline.hpp"
#include "lfmc/random_source/antithetic_random_source.hpp"
#include "lfmc/random_source/pseudo_random_source.hpp"
#include "lfmc/stochastic_process/stochastic_process.hpp"
#include "lfmc/strategy/strategy_runner.hpp"

#include <memory>
#include <random>
#include <string>

namespace lfmc {

template <StochasticProcess SP, NumericalScheme<SP> NS> class StrategyFactory {
  public:
    static std::unique_ptr<StrategyRunner<SP, NS>>
    create_pseudo_random_strategy(SP process, NS scheme, std::unique_ptr<Payoff> payoff,
                                  unsigned seed = std::random_device{}()) {
        auto rs = std::make_unique<PseudoRandomSource>(seed);
        auto pg = std::make_unique<PathGenerator<SP, NS>>(process, scheme);
        auto est = std::make_unique<MonteCarloEstimator>();

        Pipeline<SP, NS> pipeline(std::move(rs), std::move(pg), std::move(payoff), std::move(est));

        return std::make_unique<StrategyRunner<SP, NS>>(std::move(pipeline), "PseudoRandom");
    }

    static std::unique_ptr<StrategyRunner<SP, NS>>
    create_antithetic_strategy(SP process, NS scheme, std::unique_ptr<Payoff> payoff,
                               unsigned seed = std::random_device{}()) {
        auto rs = std::make_unique<AntitheticRandomSource>(seed);
        auto pg = std::make_unique<PathGenerator<SP, NS>>(process, scheme);
        auto est = std::make_unique<MonteCarloEstimator>();

        Pipeline<SP, NS> pipeline(std::move(rs), std::move(pg), std::move(payoff), std::move(est));

        return std::make_unique<StrategyRunner<SP, NS>>(std::move(pipeline), "Antithetic");
    }

    static std::unique_ptr<StrategyRunner<SP, NS>>
    create_control_variate_strategy(SP process, NS scheme, std::unique_ptr<Payoff> target_payoff,
                                    std::unique_ptr<Payoff> control_payoff,
                                    double control_expectation,
                                    unsigned seed = std::random_device{}()) {
        auto rs = std::make_unique<PseudoRandomSource>(seed);
        auto pg = std::make_unique<PathGenerator<SP, NS>>(process, scheme);
        auto po = std::make_unique<ControlVariatePayoff>(std::move(target_payoff),
                                                         std::move(control_payoff));
        auto est = std::make_unique<ControlVariateEstimator>(control_expectation);

        Pipeline<SP, NS> pipeline(std::move(rs), std::move(pg), std::move(po), std::move(est));

        return std::make_unique<StrategyRunner<SP, NS>>(std::move(pipeline), "ControlVariate");
    }

    // Add more factory methods as we make/implement more VR strategies:
    // - create_quasi_monte_carlo_strategy()
    // - create_importance_sampling_strategy()
    // etc.
};

} // namespace lfmc