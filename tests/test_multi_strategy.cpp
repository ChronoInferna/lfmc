#include "lfmc/numerical_scheme/euler_maruyama.hpp"
#include "lfmc/payoff/european_payoffs.hpp"
#include "lfmc/stochastic_process/geometric_brownian_motion.hpp"
#include "lfmc/strategy/multi_strategy_runner.hpp"
#include "lfmc/strategy/strategy_factory.hpp"

#include <catch2/catch_test_macros.hpp>
#include <memory>

using namespace lfmc;

// TODO add more test cases for different option types, more strategies, error handling, etc.

TEST_CASE("Multi-strategy runner compares variance reduction techniques") {
    GeometricBrownianMotion gbm(0.05, 0.20, 100.0);
    EulerMaruyama<GeometricBrownianMotion> euler;

    constexpr size_t steps = 252;
    constexpr double T = 1.0;
    constexpr size_t warmup_iterations = 100;

    MultiStrategyRunner<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>> runner;

    // Add pseudo-random strategy
    runner.add_strategy(
        StrategyFactory<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>>::
            create_pseudo_random_strategy(gbm, euler, std::make_unique<EuropeanCall>(100.0), 42u));

    // Add antithetic variates strategy
    runner.add_strategy(
        StrategyFactory<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>>::
            create_antithetic_strategy(gbm, euler, std::make_unique<EuropeanCall>(100.0), 42u));

    auto result = runner.run_warmup(steps, T, warmup_iterations);

    REQUIRE(result.has_value());
    REQUIRE(runner.strategy_count() == 2);

    auto metrics = runner.get_all_metrics();

    REQUIRE(metrics.size() == 2);
    REQUIRE(metrics[0].samples > 0);
    REQUIRE(metrics[1].samples > 0);
    REQUIRE(metrics[0].variance > 0.0);
    REQUIRE(metrics[1].variance > 0.0);

    size_t best_idx = runner.best_strategy_index();
    INFO("Best strategy: " << metrics[best_idx].name
                           << " with variance: " << metrics[best_idx].variance);

    REQUIRE(best_idx < metrics.size());
}

TEST_CASE("Multi-strategy runner with control variates") {
    // TODO implement once control variate bugs are fixed
}

TEST_CASE("Multi-strategy runner handles empty strategy list") {
    // TODO test error handling for no strategies added
}

TEST_CASE("Multi-strategy runner can run multiple warmup cycles") {
    // TODO test running warmup multiple times and verify metrics update
}