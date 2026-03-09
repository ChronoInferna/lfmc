#include "lfmc/pipeline.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace lfmc;

// TODO add more test cases, e.g. for convergence criteria, error handling, etc. and dummy types

TEST_CASE("Pipeline runs until estimator converges") {
    GeometricBrownianMotion gbm(0.05, 0.2, 100.0);
    EulerMaruyama<GeometricBrownianMotion> euler;

    auto rs = std::make_unique<PseudoRandomSource>();
    auto pg = std::make_unique<
        PathGenerator<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>>>(gbm, euler);
    auto po = std::make_unique<EuropeanCall>(100.0);
    auto est = std::make_unique<MonteCarloEstimator>();

    Pipeline<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>> pipeline(
        std::move(rs), std::move(pg), std::move(po), std::move(est));

    auto result = pipeline.run(10, 1.0);

    INFO("Pipeline result: " << (result.has_value() ? std::to_string(result.value())
                                                    : result.error()));
    REQUIRE(result.has_value());
    REQUIRE(result.value() > 0.0); // Price should be positive
}

TEST_CASE("Pipeline with antithetic variates") {
    GeometricBrownianMotion gbm(0.05, 0.2, 100.0);
    EulerMaruyama<GeometricBrownianMotion> euler;

    auto rs = std::make_unique<AntitheticRandomSource>();
    auto pg = std::make_unique<
        PathGenerator<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>>>(gbm, euler);
    auto po = std::make_unique<EuropeanCall>(100.0);
    auto est = std::make_unique<MonteCarloEstimator>();

    Pipeline<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>> pipeline(
        std::move(rs), std::move(pg), std::move(po), std::move(est));

    auto result = pipeline.run(10, 1.0);

    INFO("Pipeline result: " << (result.has_value() ? std::to_string(result.value())
                                                    : result.error()));
    REQUIRE(result.has_value());
    REQUIRE(result.value() > 0.0); // Price should be positive
}

TEST_CASE("Pipeline with control variates") {
    GeometricBrownianMotion gbm(0.05, 0.2, 100.0);
    EulerMaruyama<GeometricBrownianMotion> euler;

    auto rs = std::make_unique<PseudoRandomSource>();
    auto pg = std::make_unique<
        PathGenerator<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>>>(gbm, euler);
    auto po1 = std::make_unique<EuropeanCall>(100.0);
    auto po2 = std::make_unique<EuropeanCall>(
        100.0); // Control variate payoff (same as target for simplicity)
    auto po = std::make_unique<ControlVariatePayoff>(std::move(po1), std::move(po2));
    auto est = std::make_unique<ControlVariateEstimator>(
        gbm.initial() * std::exp(gbm.mu() * 1.0)); // Control variate expectation

    Pipeline<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>> pipeline(
        std::move(rs), std::move(pg), std::move(po), std::move(est));

    auto result = pipeline.run(10, 1.0);

    INFO("Pipeline result: " << (result.has_value() ? std::to_string(result.value())
                                                    : result.error()));
    REQUIRE(result.has_value());
    REQUIRE(result.value() > 0.0); // Price should be positive
                                   // TODO Currently around 100, but should be around 10 - need to
                                   // debug control variate implementation
}

TEST_CASE("Pipeline stops early if estimator add_payoffs fails") {
    // TODO
}

TEST_CASE("Pipeline performs correct number of iterations") {
    // TODO
}
