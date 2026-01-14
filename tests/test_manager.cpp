#include "lfmc/Manager.hpp"

#include <catch2/catch_test_macros.hpp>
#include <memory>

TEST_CASE("Manager initializes", "[Manager]") {
    // Define a variant type for variance reduction strategies
    using VRs = std::variant<lfmc::NoVarianceReduction>;
    // using VRs = std::variant<>; // Used for testing concept

    SECTION("Initialization with temporaries") {
        REQUIRE_NOTHROW(lfmc::Manager<lfmc::GeometricBrownianMotion,
                                      lfmc::EulerMaruyama<lfmc::GeometricBrownianMotion>, VRs>(
            lfmc::GeometricBrownianMotion{.mu = 0.05, .sigma = 0.2},
            lfmc::EulerMaruyama<lfmc::GeometricBrownianMotion>(),
            std::make_unique<lfmc::NoVarianceReduction>()));
    };

    SECTION("Initialization with named objects") {
        lfmc::GeometricBrownianMotion gbm{.mu = 0.05, .sigma = 0.2};
        lfmc::EulerMaruyama<lfmc::GeometricBrownianMotion> scheme;
        std::unique_ptr<lfmc::VarianceReductionStrategy> vrStrategy =
            std::make_unique<lfmc::NoVarianceReduction>();

        REQUIRE_NOTHROW(lfmc::Manager<lfmc::GeometricBrownianMotion,
                                      lfmc::EulerMaruyama<lfmc::GeometricBrownianMotion>, VRs>(
            gbm, scheme, std::move(vrStrategy)));
    };
}
