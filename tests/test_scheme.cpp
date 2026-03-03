#include "lfmc/numerical_scheme.hpp"
#include "lfmc/stochastic_process.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

namespace test1 {

struct ValidProcess {
    double initial() const noexcept {
        return 1.0;
    }
    double drift(double x) const noexcept {
        return x;
    }
    double diffusion(double x) const noexcept {
        return x;
    }
};

template <lfmc::StochasticProcess P> struct ValidScheme {
    double step(P const& process, double x, double dt, double dW) const noexcept {
        return x + process.drift(x) * dt + process.diffusion(x) * dW * std::sqrt(dt);
    }
};

struct InvalidSchemeNoStep {
    // Missing step function
};

template <lfmc::StochasticProcess P> struct InvalidSchemeWrongStep {
    double step(P const& process, double x, double dt, double dW) const noexcept {
        return x + process.drift(x) * dt; // Missing diffusion term
    }
};

} // namespace test1

TEST_CASE("NumericalScheme concept works correctly", "[NumericalScheme]") {
    using namespace test1;

    // TODO
    // REQUIRE(lfmc::NumericalScheme<ValidProcess, ValidScheme<ValidProcess>>);
    // REQUIRE_FALSE(lfmc::NumericalScheme<InvalidSchemeNoStep, ValidProcess>);
    // REQUIRE_FALSE(lfmc::NumericalScheme<InvalidSchemeWrongStep<ValidProcess>, ValidProcess>);
};

TEST_CASE("EulerMaruyama computes next step correctly", "[NumericalScheme][EulerMaruyama]") {
    double x = 100.0;
    double dt = 0.01;
    double dW = 0.05;

    lfmc::GeometricBrownianMotion gbm(0.1, 0.2, 100.0);
    lfmc::EulerMaruyama<lfmc::GeometricBrownianMotion> scheme;

    double expectedNextX = x + gbm.drift(x, x) * dt + gbm.diffusion(x, x) * dW * std::sqrt(dt);
    double computedNextX = scheme.step(gbm, x, x, dt, dW);

    REQUIRE_THAT(computedNextX, WithinAbs(expectedNextX, 1e-10));
}

TEST_CASE("EulerMaruyama works with different stochastic processes",
          "[NumericalScheme][EulerMaruyama]") {
    // TODO
    // struct CustomProcess {
    //     double initial() const noexcept {
    //         return 50.0;
    //     }
    //     double drift(double x) const noexcept {
    //         return 2.0 * x;
    //     }
    //     double diffusion(double x) const noexcept {
    //         return 0.5 * x;
    //     }
    // };
    //
    // CustomProcess process;
    // lfmc::EulerMaruyama<CustomProcess> scheme;
    //
    // double x = 50.0;
    // double dt = 0.02;
    // double dW = 0.03;
    //
    // double expectedNextX = x + process.drift(x) * dt + process.diffusion(x) * dW * std::sqrt(dt);
    // double computedNextX = scheme.step(process, x, dt, dW);
    //
    // REQUIRE_THAT(computedNextX, WithinAbs(expectedNextX, 1e-10));
}
