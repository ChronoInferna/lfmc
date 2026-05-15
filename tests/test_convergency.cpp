#include "lfmc/estimator/monte_carlo_estimator.hpp"
#include "lfmc/numerical_scheme/euler_maruyama.hpp"
#include "lfmc/payoff/asian_payoffs.hpp"
#include "lfmc/payoff/barrier_payoffs.hpp"
#include "lfmc/payoff/european_payoffs.hpp"
#include "lfmc/payoff/lookback_payoffs.hpp"
#include "lfmc/pipeline/pipeline.hpp"
#include "lfmc/random_source/pseudo_random_source.hpp"
#include "lfmc/stochastic_process/geometric_brownian_motion.hpp"
#include "lfmc/timing/timing.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace lfmc;
using Catch::Matchers::WithinAbs;

// Black-Scholes closed-form helpers for ground truth
namespace bs {

static double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// Standard European Call  used to sanity check our GBM setup
double european_call(double S, double K, double r, double sigma, double T) {
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    return S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
}

// Geometric Asian Call closed-form (Kemna-Vorst approximation)
// Used as ground truth since arithmetic Asian has no closed form
double geometric_asian_call(double S, double K, double r, double sigma, double T, int n) {
    double sigma_adj = sigma * std::sqrt((2.0 * n + 1.0) / (6.0 * (n + 1.0)));
    double r_adj = 0.5 * (r - 0.5 * sigma * sigma) + 0.5 * sigma_adj * sigma_adj;
    double d1 =
        (std::log(S / K) + (r_adj + 0.5 * sigma_adj * sigma_adj) * T) / (sigma_adj * std::sqrt(T));
    double d2 = d1 - sigma_adj * std::sqrt(T);
    return std::exp(-r * T) * (S * std::exp(r_adj * T) * norm_cdf(d1) - K * norm_cdf(d2));
}

// Up-and-Out Call closed form (continuous barrier, no dividends)
double up_and_out_call(double S, double K, double H, double r, double sigma, double T) {
    if (S >= H)
        return 0.0;
    double vanilla = european_call(S, K, r, sigma, T);

    double lambda = (r + 0.5 * sigma * sigma) / (sigma * sigma);
    double d1 =
        (std::log(H * H / (S * K)) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    double reflection =
        std::pow(H / S, 2.0 * lambda) * (S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2));

    return vanilla - reflection;
}

} // namespace bs

// Convergence runner: runs pipeline at increasing sample counts and prints a table of results vs
// ground truth
struct ConvergenceResult {
    std::size_t samples;
    double estimate;
    double ground_truth;
    double abs_error;
    long long elapsed_ms;
};

// We run the pipeline once per sample tier by re-seeding with
// a fixed seed for reproducibility
template <typename PayoffFactory>
std::vector<ConvergenceResult>
run_convergence(const std::string& label, PayoffFactory make_payoff, double ground_truth,
                const std::vector<std::size_t>& sample_tiers, double S0 = 100.0, double mu = 0.05,
                double sigma = 0.20, std::size_t steps = 252, double T = 1.0) {
    std::cout << "\n=== " << label << " ===\n";
    std::cout << std::setw(10) << "Samples" << std::setw(14) << "Estimate" << std::setw(14)
              << "Truth" << std::setw(14) << "AbsError" << std::setw(12) << "Time(ms)" << "\n"
              << std::string(64, '-') << "\n";

    std::vector<ConvergenceResult> results;

    for (std::size_t n : sample_tiers) {
        GeometricBrownianMotion gbm(mu, sigma, S0);
        EulerMaruyama<GeometricBrownianMotion> euler;

        // Fixed seed for reproducibility
        auto rs = std::make_unique<PseudoRandomSource>(42u);
        auto pg = std::make_unique<
            PathGenerator<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>>>(gbm,
                                                                                            euler);
        auto po = make_payoff();
        auto est = std::make_unique<MonteCarloEstimator>(); // see note below

        Pipeline<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>> pipeline(
            std::move(rs), std::move(pg), std::move(po), std::move(est));

        long long elapsed = 0;
        double estimate = 0.0;
        {
            ScopedTimer timer(elapsed);
            auto res = pipeline.run(steps, T);
            REQUIRE(res.has_value());
            estimate = res.value();
        }

        double err = std::abs(estimate - ground_truth);
        results.push_back({n, estimate, ground_truth, err, elapsed});

        std::cout << std::setw(10) << n << std::setw(14) << std::fixed << std::setprecision(4)
                  << estimate << std::setw(14) << ground_truth << std::setw(14) << err
                  << std::setw(12) << elapsed << "\n";
    }
    return results;
}

// Shared parameters
static constexpr double S0 = 100.0;
static constexpr double K = 100.0;
static constexpr double B_UP = 120.0; // Up-and-out barrier
static constexpr double B_DN = 80.0;  // Down-and-in barrier
static constexpr double MU = 0.05;
static constexpr double SIGMA = 0.20;
static constexpr double T = 1.0;
static constexpr int STEPS = 252;

static const std::vector<std::size_t> TIERS = {1000, 5000, 10000, 50000, 100000, 500000};

// Tests

TEST_CASE("Asian Call convergence", "[exotic][convergence][asian]") {
    double truth = bs::geometric_asian_call(S0, K, MU, SIGMA, T, STEPS);

    auto results = run_convergence(
        "Arithmetic Asian Call (truth = geometric approx)",
        []() { return std::make_unique<AsianCall>(K); }, truth, TIERS);

    // At 100k samples we should be within $0.50 of the approximation
    auto& last = results.back();
    REQUIRE_THAT(last.estimate, WithinAbs(last.ground_truth, 0.50));
}

TEST_CASE("Up-and-Out Barrier Call convergence", "[exotic][convergence][barrier]") {
    double truth = bs::up_and_out_call(S0, K, B_UP, MU, SIGMA, T);

    auto results = run_convergence(
        "Up-and-Out Barrier Call", []() { return std::make_unique<UpAndOutCall>(K, B_UP); }, truth,
        TIERS);

    auto& last = results.back();
    REQUIRE_THAT(last.estimate, WithinAbs(last.ground_truth, 0.50));
}

TEST_CASE("Lookback Call convergence", "[exotic][convergence][lookback]") {
    // No simple closed form for discrete lookback  we use the
    // large-sample MC estimate itself as a self-consistency check
    // and just verify convergence tightens with more samples
    auto results = run_convergence(
        "Lookback Call (floating strike)", []() { return std::make_unique<LookbackCall>(); },
        0.0, // placeholder  see check below
        TIERS);

    // Check that later estimates are closer to each other than early ones
    // i.e. the std deviation of the last 3 tiers < first 3 tiers
    auto spread = [&](int a, int b) { return std::abs(results[b].estimate - results[a].estimate); };
    REQUIRE(spread(3, 5) < spread(0, 2));
}

TEST_CASE("Sanity check: European Call matches Black-Scholes", "[sanity][european]") {
    double bs_price = bs::european_call(S0, K, MU, SIGMA, T);

    auto results = run_convergence(
        "European Call (BS sanity check)", []() { return std::make_unique<EuropeanCall>(K); },
        bs_price, TIERS);

    auto& last = results.back();
    REQUIRE_THAT(last.estimate, WithinAbs(last.ground_truth, 0.20));
}
