#include "lfmc/stochastic_process/geometric_brownian_motion.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

TEST_CASE("StochasticProcess concept works correctly", "[StochasticProcess]") {
    // TODO
    // struct ValidProcess {
    //     double initial() const noexcept {
    //         return 1.0;
    //     }
    //     double drift(double x) const noexcept {
    //         return x;
    //     }
    //     double diffusion(double x) const noexcept {
    //         return x;
    //     }
    // };
    //
    // struct InvalidProcessNoDrift {
    //     double diffusion(double x) const noexcept {
    //         return x;
    //     }
    // };
    //
    // struct InvalidProcessWrongDiffusion {
    //     double drift(double x) const noexcept {
    //         return x;
    //     }
    //     int diffusion(double x) const noexcept {
    //         return static_cast<int>(x);
    //     }
    // };
    //
    // REQUIRE(lfmc::StochasticProcess<ValidProcess>);
    // REQUIRE_FALSE(lfmc::StochasticProcess<InvalidProcessNoDrift>);
    // REQUIRE_FALSE(lfmc::StochasticProcess<InvalidProcessWrongDiffusion>);
}

TEST_CASE("GeometricBrownianMotion computes drift and diffusion correctly",
          "[StochasticProcess][GeometricBrownianMotion]") {
    lfmc::GeometricBrownianMotion gbm(0.1, 0.2, 100.0);

    double x = 100.0;
    double expectedDrift = 0.1 * x;
    double expectedDiffusion = 0.2 * x;

    REQUIRE_THAT(gbm.drift(x, x), WithinAbs(expectedDrift, 1e-10));
    REQUIRE_THAT(gbm.diffusion(x, x), WithinAbs(expectedDiffusion, 1e-10));
}
