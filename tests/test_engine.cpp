#include "lfmc/engine/engine.hpp"
#include "lfmc/numerical_scheme/euler_maruyama.hpp"
#include "lfmc/payoff/european_payoffs.hpp"
#include "lfmc/stochastic_process/geometric_brownian_motion.hpp"
#include "lfmc/strategy/strategy_factory.hpp"

#include <catch2/catch_test_macros.hpp>
#include <memory>

using namespace lfmc;

// TODO add more test cases for different option types, more strategies, error handling, etc.

TEST_CASE("Engine compares variance reduction techniques", "[engine]") {
    GeometricBrownianMotion gbm(0.05, 0.20, 100.0);
    EulerMaruyama<GeometricBrownianMotion> euler;

    constexpr size_t steps = 252;
    constexpr double T = 1.0;
    constexpr size_t warmup_iterations = 100;

    Engine<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>> engine;

    // Add pseudo-random strategy
    auto pseudo_random_strategy_result =
        StrategyFactory<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>>::
            create_pseudo_random_strategy(gbm, euler, std::make_unique<EuropeanCall>(100.0), 42u);
    if (!pseudo_random_strategy_result) {
        FAIL("Failed to create pseudo-random strategy: " << pseudo_random_strategy_result.error());
    }
    engine.add_strategy(std::move(pseudo_random_strategy_result.value()));

    // Add antithetic variates strategy
    auto antithetic_strategy_result =
        StrategyFactory<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>>::
            create_antithetic_strategy(gbm, euler, std::make_unique<EuropeanCall>(100.0), 42u);
    if (!antithetic_strategy_result) {
        FAIL("Failed to create antithetic strategy: " << antithetic_strategy_result.error());
    }
    engine.add_strategy(std::move(antithetic_strategy_result.value()));

    auto result = engine.run_warmup(steps, T, warmup_iterations);

    REQUIRE(result.has_value());
    REQUIRE(engine.strategy_count() == 2);

    auto metrics = engine.get_all_metrics();

    REQUIRE(metrics.size() == 2);
    REQUIRE(metrics[0].samples > 0);
    REQUIRE(metrics[1].samples > 0);
    REQUIRE(metrics[0].variance > 0.0);
    REQUIRE(metrics[1].variance > 0.0);

    auto best_idx = engine.get_best_strategy_index();
    if (!best_idx) {
        FAIL("Failed to get best strategy index: " << best_idx.error());
    }
    size_t idx = best_idx.value();

    REQUIRE(idx < metrics.size());
}

TEST_CASE("Multi-strategy engine with control variates") {
    // TODO implement once control variate bugs are fixed
}

TEST_CASE("Multi-strategy engine handles empty strategy list") {
    // TODO test error handling for no strategies added
}

TEST_CASE("Multi-strategy engine can run multiple warmup cycles") {
    // TODO test running warmup multiple times and verify metrics update
}
