#include "lfmc/stochastic_process.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

TEST_CASE("GeometricBrownianMotion computes drift and diffusion correctly",
          "[StochasticProcess][GeometricBrownianMotion]") {
    lfmc::GeometricBrownianMotion gbm(0.1, 0.2, 100.0);

    double x = 100.0;
    double expectedDrift = 0.1 * x;
    double expectedDiffusion = 0.2 * x;

    REQUIRE_THAT(gbm.drift(x, x), WithinAbs(expectedDrift, 1e-10));
    REQUIRE_THAT(gbm.diffusion(x, x), WithinAbs(expectedDiffusion, 1e-10));
}
