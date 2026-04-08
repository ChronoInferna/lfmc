#include "lfmc/adaptive_estimator.hpp"
#include "lfmc/numerical_scheme.hpp"
#include "lfmc/payoff.hpp"
#include "lfmc/stochastic_process.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numeric>

using namespace lfmc;

// --- Black-Scholes closed-form price ----------------------------------------
// Used as the ground truth to verify unbiasedness.

static double standard_normal_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

static double bs_call(double S, double K, double r, double sigma, double T) {
    const double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    const double d2 = d1 - sigma * std::sqrt(T);
    return S * standard_normal_cdf(d1) - K * std::exp(-r * T) * standard_normal_cdf(d2);
}

// --- Shared test parameters -------------------------------------------------
// ATM European call on GBM: S0=100, K=100, mu=0.05, sigma=0.2, T=1
//
// The library computes the UNDISCOUNTED expected payoff E[max(S_T-K,0)].
// With mu = risk-free rate r, the physical and risk-neutral measures coincide,
// so the correct reference is:
//
//   E[max(S_T-K,0)] = BS_call(S0,K,r,sigma,T) * exp(r*T)
//                   ≈ 10.45 * exp(0.05) ≈ 10.99
//
// The Black-Scholes discounted price is NOT the right comparison target here.

static constexpr double S0    = 100.0;
static constexpr double MU    = 0.05;
static constexpr double SIGMA = 0.2;
static constexpr double K     = 100.0;
static constexpr double T     = 1.0;
static constexpr size_t STEPS = 52; // weekly steps
// control: E[S_T] = S0 * exp(mu*T)
static const double CONTROL_MEAN = S0 * std::exp(MU * T);

// Build the standard four-strategy suite
static std::vector<NamedStrategy> build_strategies() {
    GeometricBrownianMotion gbm{MU, SIGMA, S0};
    EulerMaruyama<GeometricBrownianMotion> euler;

    return {
        {"plain_mc",
         make_plain_mc_sampler(gbm, euler, std::make_shared<EuropeanCall>(K), STEPS, T)},

        {"antithetic",
         make_antithetic_sampler(gbm, euler, std::make_shared<EuropeanCall>(K), STEPS, T)},

        // Control variate: target = call payoff, control = S_T (EuropeanCall(0) = max(S_T,0) = S_T
        // for GBM), E[S_T] = S0 * exp(mu*T)
        {"control_variate",
         make_control_variate_sampler(gbm, euler, std::make_shared<EuropeanCall>(K),
                                       std::make_shared<EuropeanCall>(0.0), CONTROL_MEAN, STEPS,
                                       T)},

        {"antithetic_cv",
         make_antithetic_cv_sampler(gbm, euler, std::make_shared<EuropeanCall>(K),
                                     std::make_shared<EuropeanCall>(0.0), CONTROL_MEAN, STEPS, T)},
    };
}

// --- Tests ------------------------------------------------------------------

TEST_CASE("ASVR: estimate is close to undiscounted expected payoff") {
    // ASVR with 20 000 total samples.
    // True target: E[max(S_T-K,0)] = BS_call * exp(r*T) ≈ 10.99
    // We use a wide tolerance (±1.0) since the estimator has no discounting
    // and discretization error from Euler-Maruyama adds a small systematic bias.

    const double undiscounted_target = bs_call(S0, K, MU, SIGMA, T) * std::exp(MU * T);

    ASVRConfig cfg;
    cfg.exploration_fraction = 0.1;
    cfg.n_threads = 4;

    auto result = AdaptiveVarianceReduction::run(build_strategies(), 20'000, cfg);

    INFO("Undiscounted target: " << undiscounted_target);
    INFO("ASVR estimate: " << result.estimate << " ± " << result.estimated_stderr);
    REQUIRE(result.estimate > 0.0);
    // Within ±1.0 of the true undiscounted payoff (very conservative bound)
    REQUIRE(std::abs(result.estimate - undiscounted_target) < 1.0);
}

TEST_CASE("ASVR: variance reduction ratio > 1 (beats plain MC)") {
    // The core claim: same total sample budget, lower estimated variance.
    // For ATM European calls on GBM, antithetic and CV each achieve >50% variance
    // reduction. Even accounting for the 10% exploration overhead, ASVR should
    // easily achieve VR ratio > 1.

    ASVRConfig cfg;
    cfg.exploration_fraction = 0.1;
    cfg.n_threads = 4;

    auto result = AdaptiveVarianceReduction::run(build_strategies(), 20'000, cfg);

    INFO("VR ratio: " << result.variance_reduction_ratio);
    INFO("ASVR variance: " << result.estimated_variance);
    INFO("Plain MC variance: " << result.plain_mc_variance_estimate);
    REQUIRE(result.variance_reduction_ratio > 1.0);
}

TEST_CASE("ASVR: precision weights sum to 1") {
    ASVRConfig cfg;
    cfg.exploration_fraction = 0.1;
    cfg.n_threads = 4;

    auto result = AdaptiveVarianceReduction::run(build_strategies(), 10'000, cfg);

    double weight_sum = 0.0;
    for (double w : result.precision_weights)
        weight_sum += w;

    REQUIRE(weight_sum == Catch::Approx(1.0).epsilon(1e-10));
}

TEST_CASE("ASVR: exploration identifies a VR strategy as best") {
    // The best strategy (lowest exploration variance) should NOT be plain_mc.
    // For a European call, antithetic/CV should always outperform plain MC.

    ASVRConfig cfg;
    cfg.exploration_fraction = 0.15;
    cfg.n_threads = 4;

    auto result = AdaptiveVarianceReduction::run(build_strategies(), 20'000, cfg);

    // Find which strategy got the highest weight
    const auto& weights = result.precision_weights;
    const size_t best_k = static_cast<size_t>(
        std::max_element(weights.begin(), weights.end()) - weights.begin());

    INFO("Best strategy: " << result.exploration_stats[best_k].name);
    INFO("Weights: ");
    for (size_t k = 0; k < weights.size(); ++k)
        INFO("  " << result.exploration_stats[k].name << ": " << weights[k]);

    // plain_mc is index 0; it should not have the highest weight
    REQUIRE(best_k != 0);
}

TEST_CASE("ASVR: sample budget splits correctly") {
    const size_t N = 10'000;
    const size_t K = 4;

    ASVRConfig cfg;
    cfg.exploration_fraction = 0.1;
    cfg.n_threads = 4;
    cfg.min_exploration_per_strategy = 0; // allow the fraction to fully control the split

    auto result = AdaptiveVarianceReduction::run(build_strategies(), N, cfg);

    // Exploration: K * floor(0.1 * N / K) = K * floor(250) = 1000
    // Exploitation: 9000
    REQUIRE(result.n_exploration + result.n_exploitation == N);
    REQUIRE(result.total_samples == N);

    // Exploitation counts sum to n_exploitation
    size_t exploit_sum = 0;
    for (size_t c : result.exploitation_counts)
        exploit_sum += c;
    REQUIRE(exploit_sum == result.n_exploitation);
}

// --- Empirical variance comparison ------------------------------------------
//
// This test is the paper-quality benchmark: run both ASVR and plain MC many
// times, compute empirical MSE for each, and verify ASVR MSE < plain MC MSE.
// It is slower (~seconds) but produces the core empirical result.

TEST_CASE("ASVR: empirical MSE lower than plain MC with same sample budget", "[.slow]") {
    // Use a fixed reference: the mean of a large plain MC run as the "true" value.
    // This avoids the discounting ambiguity and tests variance reduction in isolation.
    // We compute a reference from 200k plain MC samples (essentially noiseless).
    GeometricBrownianMotion gbm{MU, SIGMA, S0};
    EulerMaruyama<GeometricBrownianMotion> euler;
    SamplerFn plain =
        make_plain_mc_sampler(gbm, euler, std::make_shared<EuropeanCall>(K), STEPS, T);

    // Reference: 200k-sample plain MC mean (variance ≈ 0 relative to RUNS-run experiment)
    auto ref_samples = plain(200'000, detail::make_seed(0, 42));
    double ref_mean = 0.0;
    for (double x : ref_samples)
        ref_mean += x;
    ref_mean /= static_cast<double>(ref_samples.size());

    const size_t N = 10'000; // samples per run
    const size_t RUNS = 100; // independent replications

    ASVRConfig cfg;
    cfg.exploration_fraction = 0.1;
    cfg.n_threads = 4;

    double mse_asvr = 0.0;
    double mse_plain = 0.0;

    for (size_t r = 0; r < RUNS; ++r) {
        // ASVR run - N total samples
        auto asvr_res = AdaptiveVarianceReduction::run(build_strategies(), N, cfg);
        double err_asvr = asvr_res.estimate - ref_mean;
        mse_asvr += err_asvr * err_asvr;

        // Plain MC run - same N samples
        auto plain_samples = plain(N, detail::make_seed(r + 10000, 99));
        double plain_mean = 0.0;
        for (double x : plain_samples)
            plain_mean += x;
        plain_mean /= static_cast<double>(plain_samples.size());
        double err_plain = plain_mean - ref_mean;
        mse_plain += err_plain * err_plain;
    }

    mse_asvr /= static_cast<double>(RUNS);
    mse_plain /= static_cast<double>(RUNS);

    INFO("Reference value:        " << ref_mean);
    INFO("Empirical MSE ASVR:     " << mse_asvr);
    INFO("Empirical MSE plain MC: " << mse_plain);
    INFO("Empirical VR ratio:     " << mse_plain / mse_asvr);

    REQUIRE(mse_asvr < mse_plain);
}
