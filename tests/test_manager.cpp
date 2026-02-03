#include "lfmc/Manager.hpp"

#include <catch2/catch_test_macros.hpp>
#include <memory>

TEST_CASE("Full Monte Carlo test", "[Manager]") {
    SECTION("Black-Scholes Function") {
        auto black_scholes_call = [](double S, double K, double r, double sigma, double T) {
            double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
            double d2 = d1 - sigma * std::sqrt(T);

            // NormCDF approx
            auto norm_cdf = [](double x) { return 0.5 * std::erfc(-x * M_SQRT1_2); };

            return S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
        };

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

        // Create manager
        lfmc::Manager<lfmc::GeometricBrownianMotion,
                      lfmc::EulerMaruyama<lfmc::GeometricBrownianMotion>, lfmc::EuropeanCall>
            manager(gbm, euler, call, std::make_unique<lfmc::NoVarianceReduction>(), S0, T, nSteps);

        // Run sims
        auto [mc_price, std_error] = manager.simulateWithError();

        // Discount to present value
        double discounted_price = mc_price * std::exp(-r * T);
        double discounted_error = std_error * std::exp(-r * T);

        // Calculate Black-Scholes for comparison
        double bs_price = black_scholes_call(S0, K, r, sigma, T);

        // Assertions on results
        REQUIRE(std::abs(discounted_price - bs_price) < 0.05); // Check absolute error is small
        REQUIRE((std::abs(discounted_price - bs_price) / bs_price * 100) <
                10.0); // Check relative error

        // 95% confidence interval
        double ci_lower = discounted_price - 1.96 * discounted_error;
        double ci_upper = discounted_price + 1.96 * discounted_error;

        // Check if Black-Scholes price is within confidence interval
        REQUIRE(bs_price >= ci_lower);
        REQUIRE(bs_price <= ci_upper);
    }
}
