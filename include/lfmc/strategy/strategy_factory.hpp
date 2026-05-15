#pragma once

#include "lfmc/estimator/control_variate_estimator.hpp"
#include "lfmc/path_generator/path_generator.hpp"
#include "lfmc/payoff/control_variate_payoffs.hpp"
#include "lfmc/payoff/payoff.hpp"
#include "lfmc/pipeline/pipeline_builder.hpp"
#include "lfmc/random_source/antithetic_random_source.hpp"
#include "lfmc/random_source/pseudo_random_source.hpp"
#include "lfmc/stochastic_process/stochastic_process.hpp"
#include "lfmc/strategy/strategy_worker.hpp"
#include "lfmc/strategy/types.hpp"

#include <expected>
#include <memory>
#include <random>

namespace lfmc {

template <StochasticProcess SP, NumericalScheme<SP> NS> class StrategyFactory {
  public:
    static std::expected<std::unique_ptr<StrategyWorker<SP, NS>>, std::string>
    create_pseudo_random_strategy(SP process, NS scheme, std::unique_ptr<Payoff> payoff,
                                  unsigned seed = std::random_device{}()) {
        auto rs = std::make_unique<PseudoRandomSource>(seed);
        auto pg = std::make_unique<PathGenerator<SP, NS>>(process, scheme);

        auto pipeline = PipelineBuilder<SP, NS>()
                            .random_source(std::move(rs))
                            .path_generator(std::move(pg))
                            .payoff(std::move(payoff))
                            .estimator()
                            .build();

        if (!pipeline) {
            return std::unexpected("Failed to build pipeline: " + pipeline.error());
        }

        return std::make_unique<StrategyWorker<SP, NS>>(std::move(pipeline.value()),
                                                        Strategy::PseudoRandom);
    }

    static std::expected<std::unique_ptr<StrategyWorker<SP, NS>>, std::string>
    create_antithetic_strategy(SP process, NS scheme, std::unique_ptr<Payoff> payoff,
                               unsigned seed = std::random_device{}()) {
        auto rs = std::make_unique<AntitheticRandomSource>(seed);
        auto pg = std::make_unique<PathGenerator<SP, NS>>(process, scheme);

        auto pipeline = PipelineBuilder<SP, NS>()
                            .random_source(std::move(rs))
                            .path_generator(std::move(pg))
                            .payoff(std::move(payoff))
                            .estimator()
                            .build();

        if (!pipeline) {
            return std::unexpected("Failed to build pipeline: " + pipeline.error());
        }

        return std::make_unique<StrategyWorker<SP, NS>>(std::move(pipeline.value()),
                                                        Strategy::Antithetic);
    }

    static std::expected<std::unique_ptr<StrategyWorker<SP, NS>>, std::string>
    create_control_variate_strategy(SP process, NS scheme, std::unique_ptr<Payoff> target_payoff,
                                    std::unique_ptr<Payoff> control_payoff,
                                    double control_expectation,
                                    unsigned seed = std::random_device{}()) {
        auto rs = std::make_unique<PseudoRandomSource>(seed);
        auto pg = std::make_unique<PathGenerator<SP, NS>>(process, scheme);
        auto po = std::make_unique<ControlVariatePayoff>(std::move(target_payoff),
                                                         std::move(control_payoff));
        auto est = std::make_unique<ControlVariateEstimator>(control_expectation);

        auto pipeline = PipelineBuilder<SP, NS>()
                            .random_source(std::move(rs))
                            .path_generator(std::move(pg))
                            .payoff(std::move(po))
                            .estimator(std::move(est))
                            .build();

        if (!pipeline) {
            return std::unexpected("Failed to build pipeline: " + pipeline.error());
        }

        return std::make_unique<StrategyWorker<SP, NS>>(std::move(pipeline.value()),
                                                        Strategy::ControlVariate);
    }

    // Add more factory methods as we make/implement more VR strategies:
    // - create_quasi_monte_carlo_strategy()
    // - create_importance_sampling_strategy()
    // etc., and add it to the enum
};

} // namespace lfmc
