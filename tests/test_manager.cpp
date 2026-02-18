#include "lfmc/Manager.hpp"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("Full Monte Carlo tests ", "[Manager]") {
    auto blackScholesCall = [](double S, double K, double r, double sigma, double T) {
        double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
        double d2 = d1 - sigma * std::sqrt(T);

        // NormCDF approx
        auto normCdf = [](double x) { return 0.5 * std::erfc(-x * M_SQRT1_2); };

        return S * normCdf(d1) - K * std::exp(-r * T) * normCdf(d2);
    };

    SECTION("No variance reduction") {
        // Parameters
        double S0 = 100.0;   // Initial stock price
        double K = 100.0;    // Strike price
        double r = 0.05;     // Risk-free rate (5%)
        double sigma = 0.20; // Volatility (20%)
        double T = 1.0;      // Time to maturity (1 year)
        size_t nSteps = 252; // Daily steps

        lfmc::GeometricBrownianMotion gbm(r, sigma);
        lfmc::EulerMaruyama<lfmc::GeometricBrownianMotion> euler;
        lfmc::EuropeanCall call(K);
        lfmc::State state{S0, T, nSteps};
        lfmc::ManagerConfig config{1000, 0};

        // Create manager
        lfmc::Manager<> manager(gbm, euler, call, state);

        // Run sims
        auto [mcPrice, stdError] = manager.simulateWithError(config);

        // Discount to present value
        double discountedPrice = mcPrice * std::exp(-r * T);
        double discountedError = stdError * std::exp(-r * T);

        // Calculate Black-Scholes for comparison
        double bsPrice = blackScholesCall(S0, K, r, sigma, T);

        // Assertions on results
        REQUIRE(std::abs(discountedPrice - bsPrice) < 0.05); // Check absolute error is small
        REQUIRE((std::abs(discountedPrice - bsPrice) / bsPrice * 100) <
                10.0); // Check relative error

        // 95% confidence interval
        double ciLower = discountedPrice - 1.96 * discountedError;
        double ciUpper = discountedPrice + 1.96 * discountedError;

        // Check if Black-Scholes price is within confidence interval
        REQUIRE(bsPrice >= ciLower);
        REQUIRE(bsPrice <= ciUpper);
    }

    SECTION("Antithetic variates") {
        // Parameters
        double S0 = 100.0;   // Initial stock price
        double K = 100.0;    // Strike price
        double r = 0.05;     // Risk-free rate (5%)
        double sigma = 0.20; // Volatility (20%)
        double T = 1.0;      // Time to maturity (1 year)
        size_t nSteps = 252; // Daily steps

        lfmc::GeometricBrownianMotion gbm(r, sigma);
        lfmc::EulerMaruyama<lfmc::GeometricBrownianMotion> euler;
        lfmc::EuropeanCall call(K);
        lfmc::State state{S0, T, nSteps};
        lfmc::ManagerConfig config{0, 1000};

        // Create manager with antithetic variates
        lfmc::Manager<> manager(gbm, euler, call, state);

        // Run sims with antithetic variates
        auto [mcPrice, stdError] = manager.simulateWithError(config);

        // Discount to present value
        double discountedPrice = mcPrice * std::exp(-r * T);
        double discountedError = stdError * std::exp(-r * T);

        // Calculate Black-Scholes for comparison
        double bsPrice = blackScholesCall(S0, K, r, sigma, T);

        // Assertions on results
        REQUIRE(std::abs(discountedPrice - bsPrice) < 0.05); // Check absolute error is small
        REQUIRE((std::abs(discountedPrice - bsPrice) / bsPrice * 100) <
                10.0); // Check relative error

        // 95% confidence interval
        double ciLower = discountedPrice - 1.96 * discountedError;
        double ciUpper = discountedPrice + 1.96 * discountedError;

        // Check if Black-Scholes price is within confidence interval
        REQUIRE(bsPrice >= ciLower);
        REQUIRE(bsPrice <= ciUpper);
    }
}
