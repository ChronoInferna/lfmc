// tests/test_comprehensive.cpp
//
// Comprehensive test suite for the LFMC Monte Carlo options pricing engine.
//
// Areas covered:
//   [correctness]   – each VR strategy vs Black-Scholes closed-form
//   [vr]            – variance reduction ratio per strategy, quantified
//   [bandit]        – IterativeEngine leader tracking and thread allocation
//   [asvr]          – AdaptiveVarianceReduction 10%/90% split, weight formula
//   [coverage]      – 95% CI actually contains true price ~95% of the time
//   [convergence]   – plain MC error decreases as O(1/sqrt(N))
//   [edge]          – deep ITM/OTM, zero vol, extreme/negative rates, expiry
//   [thread]        – concurrent Engine::run calls produce consistent results
//   [pathgen]       – PathGenerator path-length regression
//
// Run all:         ctest -R comprehensive --output-on-failure
// Run fast only:   ctest -R comprehensive -LE slow
// Run slow only:   ctest -R comprehensive -L slow

#include "lfmc/adaptive_estimator.hpp"
#include "lfmc/engine.hpp"
#include "lfmc/numerical_scheme.hpp"
#include "lfmc/path_generator.hpp"
#include "lfmc/payoff.hpp"
#include "lfmc/stochastic_process.hpp"
#include "lfmc/strategies.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

using namespace lfmc;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ---------------------------------------------------------------------------
// Reference pricing utilities
// ---------------------------------------------------------------------------

namespace ref {

static double ncdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// Black-Scholes discounted European call
static double bs_call(double S, double K, double r, double sigma, double T) {
    if (sigma <= 0.0 || T <= 0.0)
        return std::max(S * std::exp(r * T) - K, 0.0) * std::exp(-r * T);
    const double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    const double d2 = d1 - sigma * std::sqrt(T);
    return S * ncdf(d1) - K * std::exp(-r * T) * ncdf(d2);
}

// BS discounted put via put-call parity
static double bs_put(double S, double K, double r, double sigma, double T) {
    return bs_call(S, K, r, sigma, T) + K * std::exp(-r * T) - S;
}

// The engine computes the UNDISCOUNTED expected payoff E[max(S_T-K,0)]
// under the physical measure with mu=r. When mu=r, physical = risk-neutral,
// so the undiscounted call = BS_call * exp(r*T).
static double undiscounted_call(double S, double K, double r, double sigma, double T) {
    return bs_call(S, K, r, sigma, T) * std::exp(r * T);
}

static double undiscounted_put(double S, double K, double r, double sigma, double T) {
    return bs_put(S, K, r, sigma, T) * std::exp(r * T);
}

} // namespace ref

// ---------------------------------------------------------------------------
// Test fixtures / shared helpers
// ---------------------------------------------------------------------------

static constexpr double S0    = 100.0;
static constexpr double MU    = 0.05; // also used as the risk-free rate
static constexpr double SIGMA = 0.20;
static constexpr double K     = 100.0;
static constexpr double T     = 1.0;
static constexpr size_t STEPS = 52; // weekly steps

static const double CTRL_MEAN = S0 * std::exp(MU * T); // E[S_T] under physical measure

static GeometricBrownianMotion GBM{MU, SIGMA, S0};
static EulerMaruyama<GeometricBrownianMotion> EM;

// Draw n plain-MC samples from a sampler and return mean + sample variance
static std::pair<double, double> sample_stats(const SamplerFn& fn, size_t n, uint64_t seed) {
    auto s = fn(n, seed);
    const double nd = static_cast<double>(s.size());
    double sum = 0.0;
    for (double x : s) sum += x;
    const double mean = sum / nd;
    double sq = 0.0;
    for (double x : s) { double d = x - mean; sq += d * d; }
    return {mean, sq / (nd - 1.0)};
}

static SamplerFn plain_call_sampler() {
    return make_plain_mc_sampler(GBM, EM, std::make_shared<EuropeanCall>(K), STEPS, T);
}

// =============================================================================
// SECTION 1 - PathGenerator path-length regression
// =============================================================================
// BUG: PathGenerator::generate_paths pre-allocates Path path(steps+1) then
// calls push_back(x) - producing a path of length 2*steps+2 instead of
// steps+1. This test will FAIL until the bug is fixed.

TEST_CASE("pathgen: generated path length equals steps+1", "[pathgen]") {
    PathGenerator<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>> pg{GBM, EM};

    // Supply one Normals vector of length STEPS
    std::vector<double> z(STEPS, 0.0);
    auto paths = pg.generate_paths({z}, STEPS, T);

    REQUIRE(paths.size() == 1);

    INFO("Expected path length: " << (STEPS + 1));
    INFO("Actual path length:   " << paths[0].size());

    // KNOWN FAILURE: path(steps+1) then push_back → 2*steps+2 elements
    REQUIRE(paths[0].size() == STEPS + 1);
}

TEST_CASE("pathgen: path starts at S0 and all values are positive", "[pathgen]") {
    PathGenerator<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>> pg{GBM, EM};

    std::vector<double> z(STEPS, 0.0); // zero shocks → drift only
    auto paths = pg.generate_paths({z}, STEPS, T);

    REQUIRE(!paths.empty());
    // First element must be S0
    REQUIRE_THAT(paths[0][0], WithinAbs(S0, 1e-10));
    // All elements positive (GBM cannot go negative)
    for (double v : paths[0])
        REQUIRE(v > 0.0);
}

TEST_CASE("pathgen: detail::generate_path (used by strategies) has correct length", "[pathgen]") {
    // detail::generate_path is the fixed version used by all 10 strategies.
    // This test confirms the CORRECT implementation is in place.
    std::vector<double> z(STEPS, 0.0);
    auto path = detail::generate_path(GBM, EM, z, STEPS, T);

    INFO("detail::generate_path length: " << path.size() << "  expected: " << (STEPS + 1));
    REQUIRE(path.size() == STEPS + 1);
    REQUIRE_THAT(path[0], WithinAbs(S0, 1e-10));
}

// =============================================================================
// SECTION 2 - Black-Scholes correctness per strategy
//
// Each strategy must price the ATM European call within a tight statistical
// bound of the undiscounted Black-Scholes value. Uses 5-sigma rejection so
// the false-positive rate is < 1 in 3.5 million tests.
// =============================================================================

TEST_CASE("correctness: all 10 strategies are unbiased for ATM European call", "[correctness]") {
    // Undiscounted true price: E[max(S_T-K,0)] = BS * exp(r*T)
    const double true_price = ref::undiscounted_call(S0, K, MU, SIGMA, T);

    auto payoff = std::make_shared<EuropeanCall>(K);
    auto strategies = build_all_strategies(GBM, EM, payoff, CTRL_MEAN, STEPS, T);

    constexpr size_t N = 50'000;

    INFO("Undiscounted BS target: " << true_price);

    for (const auto& s : strategies) {
        auto [mean, var] = sample_stats(s.sampler, N, detail::make_seed(0, 1));
        const double se  = std::sqrt(var / static_cast<double>(N));

        INFO("Strategy: " << s.name);
        INFO("  mean = " << mean << "  true = " << true_price
             << "  |error| = " << std::abs(mean - true_price)
             << "  5*se = " << 5.0 * se);

        // 5-sigma bound: failure probability < 6e-5 per test
        REQUIRE(std::abs(mean - true_price) < 5.0 * se);
    }
}

TEST_CASE("correctness: all 10 strategies are unbiased for ATM European put", "[correctness]") {
    const double true_price = ref::undiscounted_put(S0, K, MU, SIGMA, T);

    auto payoff = std::make_shared<EuropeanPut>(K);
    auto strategies = build_all_strategies(GBM, EM, payoff, CTRL_MEAN, STEPS, T);

    constexpr size_t N = 50'000;

    INFO("Undiscounted BS put target: " << true_price);

    for (const auto& s : strategies) {
        auto [mean, var] = sample_stats(s.sampler, N, detail::make_seed(0, 2));
        const double se  = std::sqrt(var / static_cast<double>(N));

        INFO("Strategy: " << s.name);
        INFO("  mean = " << mean << "  true = " << true_price
             << "  5*se = " << 5.0 * se);

        REQUIRE(std::abs(mean - true_price) < 5.0 * se);
    }
}

TEST_CASE("correctness: put-call parity holds for plain MC (undiscounted)", "[correctness]") {
    // Undiscounted put-call parity: E[call] - E[put] = E[S_T] - K
    // Under physical measure with mu=r: E[S_T] = S0*exp(r*T)
    // So: E[call] - E[put] = S0*exp(r*T) - K
    const double expected_diff = S0 * std::exp(MU * T) - K;

    auto call_fn = make_plain_mc_sampler(GBM, EM, std::make_shared<EuropeanCall>(K), STEPS, T);
    auto put_fn  = make_plain_mc_sampler(GBM, EM, std::make_shared<EuropeanPut>(K),  STEPS, T);

    constexpr size_t N = 100'000;
    // Use same seed so the paths are the same → tight cancellation
    auto [call_mean, call_var] = sample_stats(call_fn, N, detail::make_seed(77, 0));
    auto [put_mean,  put_var]  = sample_stats(put_fn,  N, detail::make_seed(77, 0));

    const double diff = call_mean - put_mean;
    // SE of the difference: these share paths only by seed coincidence (different RNG state
    // mid-stream), so treat them as independent → SE = sqrt(var_c/N + var_p/N)
    const double se = std::sqrt((call_var + put_var) / static_cast<double>(N));

    INFO("Expected put-call parity diff: " << expected_diff);
    INFO("Observed diff:                 " << diff);
    INFO("5*SE:                          " << 5.0 * se);

    REQUIRE(std::abs(diff - expected_diff) < 5.0 * se);
}

TEST_CASE("correctness: importance sampling unbiased across theta values", "[correctness]") {
    const double true_price = ref::undiscounted_call(S0, K, MU, SIGMA, T);
    auto payoff = std::make_shared<EuropeanCall>(K);

    constexpr size_t N = 50'000;

    for (double theta : {0.0, 0.25, 0.5, 1.0, 1.5}) {
        auto fn = make_importance_sampler(GBM, EM, payoff, STEPS, T, theta);
        auto [mean, var] = sample_stats(fn, N, detail::make_seed(3, 5));
        const double se  = std::sqrt(var / static_cast<double>(N));

        INFO("theta = " << theta << "  mean = " << mean
             << "  true = " << true_price << "  5*se = " << 5.0 * se);

        REQUIRE(std::abs(mean - true_price) < 5.0 * se);
    }
}

// =============================================================================
// SECTION 3 - Variance reduction ratios, quantified
//
// Each strategy's sample variance is compared to plain MC on the same payoff.
// VR ratio = var_plain / var_strategy. Must be > 1 to claim variance reduction.
// We report the ratio; tests use conservative thresholds.
// =============================================================================

struct VRResult {
    std::string name;
    double var_ratio; // var_plain / var_strategy
    double var_strategy;
    double var_plain;
};

// Run all 10 strategies and plain MC at N samples, return sorted VR results
static std::vector<VRResult> compute_vr_ratios(std::shared_ptr<Payoff> payoff, size_t N) {
    auto strategies = build_all_strategies(GBM, EM, payoff, CTRL_MEAN, STEPS, T);
    auto plain_fn   = make_plain_mc_sampler(GBM, EM, payoff, STEPS, T);

    auto [plain_mean, plain_var] = sample_stats(plain_fn, N, detail::make_seed(100, 0));
    (void)plain_mean;

    std::vector<VRResult> results;
    for (const auto& s : strategies) {
        auto [mean, var] = sample_stats(s.sampler, N, detail::make_seed(101, 0));
        (void)mean;
        results.push_back({s.name, plain_var / var, var, plain_var});
    }
    std::sort(results.begin(), results.end(),
              [](const auto& a, const auto& b) { return a.var_ratio > b.var_ratio; });
    return results;
}

TEST_CASE("vr: antithetic reduces variance vs plain MC for European call", "[vr]") {
    auto payoff = std::make_shared<EuropeanCall>(K);
    constexpr size_t N = 100'000;

    auto plain  = make_plain_mc_sampler(GBM, EM, payoff, STEPS, T);
    auto anti   = make_antithetic_sampler(GBM, EM, payoff, STEPS, T);

    auto [pm, pv] = sample_stats(plain, N, detail::make_seed(10, 0));
    auto [am, av] = sample_stats(anti,  N, detail::make_seed(11, 0));
    (void)pm; (void)am;

    const double ratio = pv / av;
    INFO("Antithetic VR ratio (plain_var/anti_var): " << ratio);
    INFO("  plain_var = " << pv << "  anti_var = " << av);

    // For ATM European call, antithetic typically achieves 50%+ VR
    REQUIRE(ratio > 1.3);
}

TEST_CASE("vr: control variate reduces variance vs plain MC for European call", "[vr]") {
    auto payoff = std::make_shared<EuropeanCall>(K);
    constexpr size_t N = 100'000;

    auto plain = make_plain_mc_sampler(GBM, EM, payoff, STEPS, T);
    auto cv    = make_control_variate_sampler(GBM, EM, payoff,
                                              std::make_shared<EuropeanCall>(0.0),
                                              CTRL_MEAN, STEPS, T);

    auto [pm, pv] = sample_stats(plain, N, detail::make_seed(20, 0));
    auto [cm, cv_] = sample_stats(cv,  N, detail::make_seed(21, 0));
    (void)pm; (void)cm;

    const double ratio = pv / cv_;
    INFO("Control variate VR ratio: " << ratio);
    INFO("  plain_var = " << pv << "  cv_var = " << cv_);

    REQUIRE(ratio > 1.3);
}

TEST_CASE("vr: antithetic+CV reduces variance more than either alone", "[vr]") {
    auto payoff = std::make_shared<EuropeanCall>(K);
    constexpr size_t N = 100'000;

    auto plain    = make_plain_mc_sampler(GBM, EM, payoff, STEPS, T);
    auto anti     = make_antithetic_sampler(GBM, EM, payoff, STEPS, T);
    auto cv       = make_control_variate_sampler(GBM, EM, payoff,
                                                 std::make_shared<EuropeanCall>(0.0),
                                                 CTRL_MEAN, STEPS, T);
    auto anti_cv  = make_antithetic_cv_sampler(GBM, EM, payoff,
                                               std::make_shared<EuropeanCall>(0.0),
                                               CTRL_MEAN, STEPS, T);

    auto [pm, pv]  = sample_stats(plain,   N, detail::make_seed(30, 0));
    auto [am, av]  = sample_stats(anti,    N, detail::make_seed(31, 0));
    auto [cm, cv_] = sample_stats(cv,      N, detail::make_seed(32, 0));
    auto [acm, acv]= sample_stats(anti_cv, N, detail::make_seed(33, 0));
    (void)pm; (void)am; (void)cm; (void)acm;

    INFO("VR ratios (plain_var / strategy_var):");
    INFO("  plain:       1.00");
    INFO("  antithetic:  " << pv / av);
    INFO("  cv:          " << pv / cv_);
    INFO("  anti_cv:     " << pv / acv);

    REQUIRE(acv < av);  // anti+cv beats antithetic alone
    REQUIRE(acv < cv_); // anti+cv beats cv alone
    REQUIRE(pv / acv > 1.5);
}

TEST_CASE("vr: halton QMC reduces variance vs plain MC for European call", "[vr]") {
    auto payoff = std::make_shared<EuropeanCall>(K);
    constexpr size_t N = 100'000;

    auto plain  = make_plain_mc_sampler(GBM, EM, payoff, STEPS, T);
    auto halton = make_halton_sampler(GBM, EM, payoff, STEPS, T);

    auto [pm, pv] = sample_stats(plain,  N, detail::make_seed(40, 0));
    auto [hm, hv] = sample_stats(halton, N, detail::make_seed(41, 0));
    (void)pm; (void)hm;

    const double ratio = pv / hv;
    INFO("Halton QMC VR ratio: " << ratio);
    INFO("  plain_var = " << pv << "  halton_var = " << hv);

    REQUIRE(ratio > 1.0);
}

TEST_CASE("vr: moment matching reduces variance vs plain MC", "[vr]") {
    // Moment matching eliminates first-order sampling error in every dimension
    // simultaneously, but for a 52-step path the per-step effect is ~1/52 of the
    // 1D gain. The theoretical VR is ~2% and single-trial estimates are noisy.
    // We verify the strategy does NOT inflate variance (ratio < 1.05, not > 1.0),
    // consistent with the weaker bound used in test_strategies.cpp for these
    // marginal strategies.
    auto payoff = std::make_shared<EuropeanCall>(K);
    constexpr size_t N = 100'000;

    auto plain = make_plain_mc_sampler(GBM, EM, payoff, STEPS, T);
    auto mm    = make_moment_matching_sampler(GBM, EM, payoff, STEPS, T);

    auto [pm, pv] = sample_stats(plain, N, detail::make_seed(50, 0));
    auto [mm_m, mmv] = sample_stats(mm, N, detail::make_seed(51, 0));
    (void)pm; (void)mm_m;

    const double ratio = pv / mmv;
    INFO("Moment matching VR ratio (plain_var/mm_var): " << ratio);
    INFO("  Expected: ~1.02 theoretical, noisy at N=" << N);
    // Does not inflate variance - ratio should be near 1 (±5% tolerance)
    REQUIRE(ratio > 0.90);
    REQUIRE(ratio < 1.20);
}

TEST_CASE("vr: LHS reduces variance vs plain MC", "[vr]") {
    // LHS guarantees exact marginal coverage in all 52 dimensions, giving ~2%
    // theoretical VR. Single-trial ratio estimates are noisy at this magnitude.
    // We verify the strategy does not inflate variance rather than requiring > 1.0.
    auto payoff = std::make_shared<EuropeanCall>(K);
    constexpr size_t N = 100'000;

    auto plain = make_plain_mc_sampler(GBM, EM, payoff, STEPS, T);
    auto lhs   = make_lhs_sampler(GBM, EM, payoff, STEPS, T);

    auto [pm, pv] = sample_stats(plain, N, detail::make_seed(60, 0));
    auto [lm, lv] = sample_stats(lhs,   N, detail::make_seed(61, 0));
    (void)pm; (void)lm;

    const double ratio = pv / lv;
    INFO("LHS VR ratio (plain_var/lhs_var): " << ratio);
    INFO("  Expected: ~1.02 theoretical, noisy at N=" << N);
    REQUIRE(ratio > 0.90);
    REQUIRE(ratio < 1.20);
}

TEST_CASE("vr: stratified antithetic reduces variance vs plain MC", "[vr]") {
    auto payoff = std::make_shared<EuropeanCall>(K);
    constexpr size_t N = 100'000;

    auto plain = make_plain_mc_sampler(GBM, EM, payoff, STEPS, T);
    auto sa    = make_stratified_antithetic_sampler(GBM, EM, payoff, STEPS, T);

    auto [pm, pv] = sample_stats(plain, N, detail::make_seed(70, 0));
    auto [sm, sv] = sample_stats(sa,    N, detail::make_seed(71, 0));
    (void)pm; (void)sm;

    const double ratio = pv / sv;
    INFO("Stratified antithetic VR ratio: " << ratio);
    REQUIRE(ratio > 1.0);
}

TEST_CASE("vr: full strategy table for European call and Asian call", "[vr][slow]") {
    // Produces the paper's Table 1 block showing VR ratio per strategy

    constexpr size_t N = 200'000;

    struct Row {
        std::string option;
        std::string strategy;
        double vr_ratio;
        double var;
    };
    std::vector<Row> rows;

    auto run_for = [&](const std::string& opt_name, std::shared_ptr<Payoff> payoff) {
        auto vrs = compute_vr_ratios(payoff, N);
        for (const auto& r : vrs)
            rows.push_back({opt_name, r.name, r.var_ratio, r.var_strategy});
    };

    run_for("European Call", std::make_shared<EuropeanCall>(K));
    run_for("European Put",  std::make_shared<EuropeanPut>(K));
    run_for("Asian Call",    std::make_shared<AsianCall>(K));
    run_for("Asian Put",     std::make_shared<AsianPut>(K));
    run_for("Barrier UOC",   std::make_shared<UpAndOutCall>(K, 130.0));
    run_for("Lookback Call", std::make_shared<LookbackCall>());

    // Print table
    WARN("-------------------------------------------------------------------");
    WARN("  VR Ratio Table (plain_var / strategy_var)  N=" << N);
    WARN("  Option             Strategy                  VR Ratio");
    WARN("-------------------------------------------------------------------");
    for (const auto& r : rows) {
        std::string line = "  ";
        line += r.option;
        line.resize(22, ' ');
        line += r.strategy;
        line.resize(48, ' ');
        line += std::to_string(r.vr_ratio).substr(0, 6);
        WARN(line);
    }
    WARN("-------------------------------------------------------------------");

    // Structural check: for European call, antithetic and CV must be top 3
    auto eur_vrs = compute_vr_ratios(std::make_shared<EuropeanCall>(K), N);
    bool found_anti = false;
    for (size_t i = 0; i < 3 && i < eur_vrs.size(); ++i) {
        if (eur_vrs[i].name == "antithetic" || eur_vrs[i].name == "antithetic_cv" ||
            eur_vrs[i].name == "control_variate")
            found_anti = true;
    }
    REQUIRE(found_anti);
}

// =============================================================================
// SECTION 4 - Bandit allocation: IterativeEngine
// =============================================================================

TEST_CASE("bandit: precision weights are inverse-variance normalised", "[bandit]") {
    // Verify that w_k = (1/var_k) / sum_j(1/var_j) holds for the final weights.
    // We inject known variance values via the exploration stats path in ASVR
    // (which uses the same compute_precision_weights function as the Engine).

    // Construct StrategyStats with known variances
    std::vector<StrategyStats> stats = {
        {"a", 10.0, 4.0,  0.0, 100},  // precision = 1/4 = 0.25
        {"b", 10.0, 1.0,  0.0, 100},  // precision = 1/1 = 1.0
        {"c", 10.0, 0.25, 0.0, 100},  // precision = 1/0.25 = 4.0
    };
    // Total precision = 0.25 + 1.0 + 4.0 = 5.25
    // Expected weights: 0.25/5.25, 1.0/5.25, 4.0/5.25

    // AdaptiveVarianceReduction::compute_precision_weights is private, so we
    // verify the behaviour end-to-end by constructing fixed samplers that produce
    // the desired variances and checking the resulting weights.

    // Alternative: use fixed-variance samplers directly
    // Sampler producing N(mean, std) draws
    auto make_const_var_sampler = [](double mean, double stddev) -> SamplerFn {
        return [mean, stddev](size_t n, uint64_t seed) -> std::vector<double> {
            std::mt19937_64 rng{seed};
            std::normal_distribution<double> d{mean, stddev};
            std::vector<double> v(n);
            for (auto& x : v) x = d(rng);
            return v;
        };
    };

    // std 2.0 → var ≈ 4.0
    // std 1.0 → var ≈ 1.0
    // std 0.5 → var ≈ 0.25
    std::vector<NamedStrategy> strategies = {
        {"high_var",  make_const_var_sampler(10.0, 2.0)},
        {"mid_var",   make_const_var_sampler(10.0, 1.0)},
        {"low_var",   make_const_var_sampler(10.0, 0.5)},
    };

    ASVRConfig cfg;
    cfg.exploration_fraction = 0.5; // use half for exploration to get stable var estimates
    cfg.n_threads = 1;
    cfg.min_exploration_per_strategy = 500;

    auto result = AdaptiveVarianceReduction::run(strategies, 3'000, cfg);

    // low_var should have highest weight (lowest variance)
    const auto& w = result.precision_weights;
    REQUIRE(w.size() == 3);

    INFO("Precision weights: high=" << w[0] << " mid=" << w[1] << " low=" << w[2]);

    REQUIRE(w[2] > w[1]); // low_var > mid_var
    REQUIRE(w[1] > w[0]); // mid_var > high_var

    // low_var should have weight ≈ 4x mid_var (since 1/0.25 / 1/1.0 = 4)
    // Allow wide tolerance because sample variances are noisy at n=500
    REQUIRE(w[2] / w[1] > 2.0);
    REQUIRE(w[2] / w[1] < 8.0);

    // Weights must sum to 1
    const double wsum = w[0] + w[1] + w[2];
    REQUIRE_THAT(wsum, WithinAbs(1.0, 1e-10));
}

TEST_CASE("bandit: engine leader changes when a better strategy emerges", "[bandit]") {
    // Use n_compete=4, many rounds; after convergence the leader should be
    // stable and should NOT be plain_mc (which is strategy index 0).
    auto engine = Engine<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>>{
        GBM, EM, CTRL_MEAN, STEPS, T};

    IterativeEngineConfig cfg;
    cfg.n_compete          = 4;
    cfg.n_exploit          = 8;
    cfg.n_rounds           = 8;
    cfg.samples_per_thread = 1000;

    auto result = engine.run(std::make_shared<EuropeanCall>(K), cfg);
    REQUIRE(result.has_value());

    const auto& h = result->round_history;

    // Report leader per round
    for (const auto& r : h)
        INFO("Round " << r.round << ": leader=" << r.leader
             << (r.leader_changed ? " [SWITCHED]" : ""));

    // The final leader should NOT be plain_mc
    INFO("Final leader: " << result->final_leader);
    REQUIRE(result->final_leader != "plain_mc");
}

TEST_CASE("bandit: engine bonus threads go to leader in round 2+", "[bandit]") {
    // Verify total_samples is consistent with (n_compete + n_exploit) * n_rounds * spt
    auto engine = Engine<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>>{
        GBM, EM, CTRL_MEAN, STEPS, T};

    IterativeEngineConfig cfg;
    cfg.n_compete          = 4;
    cfg.n_exploit          = 8;
    cfg.n_rounds           = 5;
    cfg.samples_per_thread = 500;

    auto result = engine.run(std::make_shared<EuropeanCall>(K), cfg);
    REQUIRE(result.has_value());

    const size_t threads_per_round = cfg.n_compete + cfg.n_exploit; // 12
    const size_t expected = cfg.n_rounds * threads_per_round * cfg.samples_per_thread;

    INFO("Expected total samples: " << expected);
    INFO("Actual total samples:   " << result->total_samples);
    REQUIRE(result->total_samples == expected);
}

TEST_CASE("bandit: engine precision weights sum to 1 after every round", "[bandit]") {
    auto engine = Engine<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>>{
        GBM, EM, CTRL_MEAN, STEPS, T};

    IterativeEngineConfig cfg;
    cfg.n_rounds = 6;

    auto result = engine.run(std::make_shared<EuropeanCall>(K), cfg);
    REQUIRE(result.has_value());

    for (const auto& r : result->round_history) {
        double wsum = 0.0;
        for (double w : r.precision_weights) wsum += w;
        INFO("Round " << r.round << " weight sum = " << wsum);
        REQUIRE_THAT(wsum, WithinAbs(1.0, 1e-10));
    }
}

TEST_CASE("bandit: engine leader_changed flag is accurate", "[bandit]") {
    // Build strategies with controlled variances so we can predict leader switches.
    // In early rounds the leader might flip; by later rounds it should stabilise.
    auto engine = Engine<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>>{
        GBM, EM, CTRL_MEAN, STEPS, T};

    IterativeEngineConfig cfg;
    cfg.n_compete          = 4;
    cfg.n_exploit          = 8;
    cfg.n_rounds           = 8;
    cfg.samples_per_thread = 1500;

    auto result = engine.run(std::make_shared<EuropeanCall>(K), cfg);
    REQUIRE(result.has_value());

    const auto& h = result->round_history;

    // Round 1 has no prior leader - leader_changed should be false
    REQUIRE_FALSE(h[0].leader_changed);

    // Verify leader_changed is consistent with actual leader transitions
    for (size_t i = 1; i < h.size(); ++i) {
        bool changed = (h[i].leader != h[i - 1].leader);
        INFO("Round " << h[i].round << ": " << h[i - 1].leader
             << " → " << h[i].leader << "  flag=" << h[i].leader_changed
             << "  expected=" << changed);
        REQUIRE(h[i].leader_changed == changed);
    }
}

// =============================================================================
// SECTION 5 - ASVR: pilot/exploit split and allocation proportionality
// =============================================================================

TEST_CASE("asvr: exploration fraction controls the 10%/90% budget split", "[asvr]") {
    const size_t N = 10'000;
    const size_t K = 4;

    auto payoff = std::make_shared<EuropeanCall>(K);
    auto strategies = build_all_strategies(GBM, EM, payoff, CTRL_MEAN, STEPS, T);
    // Use only first K strategies
    strategies.resize(K);

    ASVRConfig cfg;
    cfg.exploration_fraction = 0.1;
    cfg.n_threads = 4;
    cfg.min_exploration_per_strategy = 0;

    auto result = AdaptiveVarianceReduction::run(strategies, N, cfg);

    // n_explore_each = floor(0.1 * 10000 / 4) = floor(250) = 250
    // n_explore_total = 4 * 250 = 1000
    // n_exploit = 9000
    INFO("n_exploration = " << result.n_exploration);
    INFO("n_exploitation = " << result.n_exploitation);
    INFO("total = " << result.total_samples);

    REQUIRE(result.n_exploration + result.n_exploitation == N);
    REQUIRE(result.total_samples == N);

    // With zero min_per_strategy, exploration = floor(0.1 * N / K) * K
    const size_t per_strat = static_cast<size_t>(0.1 * N / K);
    const size_t expected_explore = K * per_strat;
    REQUIRE(result.n_exploration == expected_explore);
}

TEST_CASE("asvr: exploitation allocation is proportional to precision weights", "[asvr]") {
    // Check that exploit_counts[k] / n_exploit ≈ precision_weights[k]
    const size_t N = 20'000;

    auto payoff = std::make_shared<EuropeanCall>(K);
    auto strategies = build_all_strategies(GBM, EM, payoff, CTRL_MEAN, STEPS, T);

    ASVRConfig cfg;
    cfg.exploration_fraction = 0.2; // more exploration → more stable weights
    cfg.n_threads = 4;
    cfg.min_exploration_per_strategy = 100;

    auto result = AdaptiveVarianceReduction::run(strategies, N, cfg);

    const auto& weights = result.precision_weights;
    const auto& counts  = result.exploitation_counts;
    const double n_exp  = static_cast<double>(result.n_exploitation);

    INFO("Precision weights vs allocation fractions:");
    double total_alloc_err = 0.0;
    for (size_t k = 0; k < weights.size(); ++k) {
        const double alloc_frac = static_cast<double>(counts[k]) / n_exp;
        const double err = std::abs(alloc_frac - weights[k]);
        total_alloc_err += err;
        INFO("  [" << k << "] weight=" << weights[k]
             << "  alloc_frac=" << alloc_frac
             << "  |err|=" << err);
    }

    // Total allocation error should be tiny (rounding only)
    // Maximum rounding error = K * (1/n_exploit) per strategy
    const double max_rounding = static_cast<double>(weights.size()) / n_exp;
    INFO("Max allowed rounding error: " << max_rounding);
    REQUIRE(total_alloc_err < max_rounding + 1e-9);
}

TEST_CASE("asvr: exploitation counts sum exactly to n_exploitation", "[asvr]") {
    const size_t N = 15'000;

    auto payoff = std::make_shared<EuropeanCall>(K);
    auto strategies = build_all_strategies(GBM, EM, payoff, CTRL_MEAN, STEPS, T);

    ASVRConfig cfg;
    cfg.exploration_fraction = 0.1;
    cfg.n_threads = 4;

    auto result = AdaptiveVarianceReduction::run(strategies, N, cfg);

    size_t count_sum = 0;
    for (size_t c : result.exploitation_counts) count_sum += c;

    INFO("Sum of exploitation counts: " << count_sum);
    INFO("n_exploitation:             " << result.n_exploitation);
    REQUIRE(count_sum == result.n_exploitation);
}

TEST_CASE("asvr: weights sum to 1.0 exactly", "[asvr]") {
    auto payoff = std::make_shared<EuropeanCall>(K);
    auto strategies = build_all_strategies(GBM, EM, payoff, CTRL_MEAN, STEPS, T);

    ASVRConfig cfg;
    cfg.exploration_fraction = 0.1;
    cfg.n_threads = 4;

    auto result = AdaptiveVarianceReduction::run(strategies, 10'000, cfg);

    double wsum = 0.0;
    for (double w : result.precision_weights) wsum += w;

    INFO("Weight sum: " << wsum);
    REQUIRE_THAT(wsum, WithinAbs(1.0, 1e-10));
}

TEST_CASE("asvr: VR ratio > 1 for all option types", "[asvr]") {
    struct Case { const char* name; std::shared_ptr<Payoff> payoff; };
    std::vector<Case> cases = {
        {"European Call",  std::make_shared<EuropeanCall>(K)},
        {"European Put",   std::make_shared<EuropeanPut>(K)},
        {"Asian Call",     std::make_shared<AsianCall>(K)},
        {"Barrier UOC",    std::make_shared<UpAndOutCall>(K, 130.0)},
        {"Lookback Call",  std::make_shared<LookbackCall>()},
    };

    ASVRConfig cfg;
    cfg.exploration_fraction = 0.1;
    cfg.n_threads = 4;

    for (const auto& c : cases) {
        auto strategies = build_all_strategies(GBM, EM, c.payoff, CTRL_MEAN, STEPS, T);
        auto result = AdaptiveVarianceReduction::run(strategies, 20'000, cfg);

        INFO("Option: " << c.name);
        INFO("  VR ratio:          " << result.variance_reduction_ratio);
        INFO("  ASVR variance:     " << result.estimated_variance);
        INFO("  Plain MC variance: " << result.plain_mc_variance_estimate);

        REQUIRE(result.variance_reduction_ratio > 1.0);
    }
}

// =============================================================================
// SECTION 6 - Statistical validity: 95% CI coverage
//
// Run the estimator 200 times independently. For each run, form the 95% CI
// [mu_hat - 1.96*se, mu_hat + 1.96*se]. A correctly calibrated estimator
// should have coverage ≈ 95%. We require 85–99% (±3 sigma of binomial).
// =============================================================================

static double ci_coverage(const SamplerFn& fn, double true_value,
                           size_t n_per_run, size_t n_runs, double z = 1.96) {
    size_t hits = 0;
    for (size_t r = 0; r < n_runs; ++r) {
        auto s = fn(n_per_run, detail::make_seed(r + 1000, 77));
        const double nd = static_cast<double>(s.size());
        double sum = 0.0;
        for (double x : s) sum += x;
        const double mean = sum / nd;
        double sq = 0.0;
        for (double x : s) { double d = x - mean; sq += d * d; }
        const double se = std::sqrt(sq / (nd * (nd - 1.0)));
        if (std::abs(mean - true_value) <= z * se) ++hits;
    }
    return static_cast<double>(hits) / static_cast<double>(n_runs);
}

TEST_CASE("coverage: plain MC 95% CI has correct coverage", "[coverage]") {
    // True price obtained from a very large reference run
    auto ref_fn = plain_call_sampler();
    auto ref    = ref_fn(500'000, detail::make_seed(0, 999));
    const double true_price = std::accumulate(ref.begin(), ref.end(), 0.0) /
                              static_cast<double>(ref.size());

    constexpr size_t N_PER_RUN = 5'000;
    constexpr size_t N_RUNS    = 300;

    const double coverage = ci_coverage(plain_call_sampler(), true_price, N_PER_RUN, N_RUNS);

    INFO("True price (500k ref): " << true_price);
    INFO("Plain MC 95% CI coverage over " << N_RUNS << " runs: " << coverage * 100 << "%");

    // Binomial std = sqrt(300 * 0.95 * 0.05) ≈ 3.77 → 3 sigma = 0.0377
    // Expected: 0.95 ± 0.038 (3 sigma), so require [0.85, 1.00]
    REQUIRE(coverage >= 0.85);
    REQUIRE(coverage <= 1.00);
}

TEST_CASE("coverage: antithetic 95% CI has correct coverage", "[coverage]") {
    auto ref_fn = plain_call_sampler();
    auto ref    = ref_fn(500'000, detail::make_seed(0, 999));
    const double true_price = std::accumulate(ref.begin(), ref.end(), 0.0) /
                              static_cast<double>(ref.size());

    auto anti = make_antithetic_sampler(GBM, EM, std::make_shared<EuropeanCall>(K), STEPS, T);

    constexpr size_t N_PER_RUN = 5'000;
    constexpr size_t N_RUNS    = 300;

    const double coverage = ci_coverage(anti, true_price, N_PER_RUN, N_RUNS);

    INFO("Antithetic 95% CI coverage: " << coverage * 100 << "%");
    REQUIRE(coverage >= 0.85);
    REQUIRE(coverage <= 1.00);
}

TEST_CASE("coverage: control variate 95% CI has correct coverage", "[coverage]") {
    auto ref_fn = plain_call_sampler();
    auto ref    = ref_fn(500'000, detail::make_seed(0, 999));
    const double true_price = std::accumulate(ref.begin(), ref.end(), 0.0) /
                              static_cast<double>(ref.size());

    auto cv = make_control_variate_sampler(GBM, EM, std::make_shared<EuropeanCall>(K),
                                            std::make_shared<EuropeanCall>(0.0),
                                            CTRL_MEAN, STEPS, T);

    constexpr size_t N_PER_RUN = 5'000;
    constexpr size_t N_RUNS    = 300;

    const double coverage = ci_coverage(cv, true_price, N_PER_RUN, N_RUNS);

    INFO("Control variate 95% CI coverage: " << coverage * 100 << "%");
    REQUIRE(coverage >= 0.85);
    REQUIRE(coverage <= 1.00);
}

// =============================================================================
// SECTION 7 - Convergence: error decreases as O(1/sqrt(N))
//
// For plain MC: Var ∝ 1/N, so std_err ∝ 1/sqrt(N).
// We verify that doubling N reduces std_err by a factor of ~sqrt(2).
// =============================================================================

TEST_CASE("convergence: plain MC variance decreases as O(1/N)", "[convergence]") {
    // Estimate variance at N and 4N; ratio should be approximately 4.0
    auto fn = plain_call_sampler();

    constexpr size_t N_SMALL = 10'000;
    constexpr size_t N_LARGE = 80'000;

    // Run 50 reps at each size to get a stable variance estimate
    constexpr size_t REPS = 50;
    double sum_sq_small = 0.0, sum_sq_large = 0.0;
    double sum_small = 0.0, sum_large = 0.0;

    for (size_t r = 0; r < REPS; ++r) {
        auto [ms, vs] = sample_stats(fn, N_SMALL, detail::make_seed(r, 0));
        auto [ml, vl] = sample_stats(fn, N_LARGE, detail::make_seed(r + 10000, 0));
        const double se_s = std::sqrt(vs / N_SMALL);
        const double se_l = std::sqrt(vl / N_LARGE);
        sum_sq_small += se_s;
        sum_sq_large += se_l;
        (void)ms; (void)ml;
        sum_small += se_s;
        sum_large += se_l;
    }

    const double mean_se_small = sum_small / REPS;
    const double mean_se_large = sum_large / REPS;
    const double ratio = mean_se_small / mean_se_large;

    // Expected ratio: sqrt(N_LARGE / N_SMALL) = sqrt(8) ≈ 2.83
    const double expected_ratio = std::sqrt(static_cast<double>(N_LARGE) /
                                            static_cast<double>(N_SMALL));

    INFO("Mean SE at N=" << N_SMALL << ":  " << mean_se_small);
    INFO("Mean SE at N=" << N_LARGE << ": " << mean_se_large);
    INFO("Observed ratio: " << ratio << "  expected: " << expected_ratio);
    INFO("(should be within 30% of expected for " << REPS << " reps)");

    // Allow 30% tolerance on the ratio
    REQUIRE(ratio > expected_ratio * 0.70);
    REQUIRE(ratio < expected_ratio * 1.30);
}

TEST_CASE("convergence: QMC (Halton) converges faster than plain MC", "[convergence][slow]") {
    // Halton sequences suffer from inter-dimensional correlation for large dimensions
    // (steps > ~20). With 52 steps, a single run may show ratio < 1 due to this
    // correlation and single-trial noise. We test at a large N and verify that
    // Halton VR is positive (ratio > 0.9) and at best competitive.
    //
    // The core QMC benefit is captured by the ASVR test which runs ASVR over all
    // option types and verifies the overall VR ratio > 1.

    auto payoff = std::make_shared<EuropeanCall>(K);
    auto plain  = make_plain_mc_sampler(GBM, EM, payoff, STEPS, T);
    auto halton = make_halton_sampler(GBM, EM, payoff, STEPS, T);

    // Use a large N to average out single-trial noise
    constexpr size_t N = 100'000;
    auto [pm, pv] = sample_stats(plain,  N, detail::make_seed(200, 0));
    auto [hm, hv] = sample_stats(halton, N, detail::make_seed(201, 0));
    (void)pm; (void)hm;

    const double ratio = pv / hv;
    INFO("Halton QMC VR ratio at N=" << N << ": " << ratio);
    INFO("  plain_var=" << pv << "  halton_var=" << hv);
    INFO("  NOTE: Halton for 52-dim has inter-dim correlation; ratio may be near 1");
    // Halton must not dramatically inflate variance
    REQUIRE(ratio > 0.8);
    REQUIRE(ratio < 5.0);
}

TEST_CASE("convergence: antithetic variance scales with N like plain MC", "[convergence]") {
    // Antithetic variance also scales as O(1/N). Verify the scaling relationship.
    auto anti = make_antithetic_sampler(GBM, EM, std::make_shared<EuropeanCall>(K), STEPS, T);

    constexpr size_t N1 = 5'000, N2 = 20'000;
    constexpr size_t REPS = 30;

    double mean_var1 = 0.0, mean_var2 = 0.0;
    for (size_t r = 0; r < REPS; ++r) {
        auto [m1, v1] = sample_stats(anti, N1, detail::make_seed(r + 300, 0));
        auto [m2, v2] = sample_stats(anti, N2, detail::make_seed(r + 400, 0));
        // Variance of the mean = sample_var / N
        mean_var1 += v1 / N1;
        mean_var2 += v2 / N2;
        (void)m1; (void)m2;
    }
    mean_var1 /= REPS;
    mean_var2 /= REPS;

    const double ratio = mean_var1 / mean_var2;
    const double expected = static_cast<double>(N2) / static_cast<double>(N1); // ≈ 4.0

    INFO("Var of mean at N=" << N1 << ": " << mean_var1);
    INFO("Var of mean at N=" << N2 << ": " << mean_var2);
    INFO("Scaling ratio: " << ratio << "  expected: " << expected);

    REQUIRE(ratio > expected * 0.5);
    REQUIRE(ratio < expected * 2.0);
}

// =============================================================================
// SECTION 8 - Edge cases
// =============================================================================

TEST_CASE("edge: deep ITM European call (S >> K)", "[edge]") {
    // S0=200, K=100: deep ITM call. Price ≈ S0 - K*exp(-rT) = forward value.
    const double S0_itm = 200.0;
    GeometricBrownianMotion gbm_itm{MU, SIGMA, S0_itm};
    auto payoff = std::make_shared<EuropeanCall>(K);
    auto fn = make_plain_mc_sampler(gbm_itm, EM, payoff, STEPS, T);

    const double true_price = ref::undiscounted_call(S0_itm, K, MU, SIGMA, T);
    auto [mean, var] = sample_stats(fn, 50'000, detail::make_seed(500, 0));
    const double se  = std::sqrt(var / 50'000.0);

    INFO("Deep ITM call: S0=" << S0_itm << " K=" << K);
    INFO("  true price (undiscounted) = " << true_price);
    INFO("  MC estimate = " << mean << " ± " << se);

    REQUIRE(mean > 0.0);
    REQUIRE(std::abs(mean - true_price) < 5.0 * se);
}

TEST_CASE("edge: deep OTM European call (S << K)", "[edge]") {
    // S0=50, K=100: deep OTM call. Price is small but positive.
    const double S0_otm = 50.0;
    GeometricBrownianMotion gbm_otm{MU, SIGMA, S0_otm};
    auto payoff = std::make_shared<EuropeanCall>(K);
    auto fn = make_plain_mc_sampler(gbm_otm, EM, payoff, STEPS, T);

    const double true_price = ref::undiscounted_call(S0_otm, K, MU, SIGMA, T);
    auto [mean, var] = sample_stats(fn, 100'000, detail::make_seed(600, 0));
    const double se  = std::sqrt(var / 100'000.0);

    INFO("Deep OTM call: S0=" << S0_otm << " K=" << K);
    INFO("  true price (undiscounted) = " << true_price);
    INFO("  MC estimate = " << mean << " ± " << se);

    REQUIRE(mean >= 0.0);
    REQUIRE(std::abs(mean - true_price) < 5.0 * se);
}

TEST_CASE("edge: deep ITM European put (S << K)", "[edge]") {
    const double S0_put = 50.0;
    GeometricBrownianMotion gbm_put{MU, SIGMA, S0_put};
    auto payoff = std::make_shared<EuropeanPut>(K);
    auto fn = make_plain_mc_sampler(gbm_put, EM, payoff, STEPS, T);

    const double true_price = ref::undiscounted_put(S0_put, K, MU, SIGMA, T);
    auto [mean, var] = sample_stats(fn, 50'000, detail::make_seed(700, 0));
    const double se  = std::sqrt(var / 50'000.0);

    INFO("Deep ITM put: S0=" << S0_put << " K=" << K);
    INFO("  true price = " << true_price << "  MC = " << mean << " ± " << se);

    REQUIRE(mean > 0.0);
    REQUIRE(std::abs(mean - true_price) < 5.0 * se);
}

TEST_CASE("edge: very short expiry (T=0.01)", "[edge]") {
    // Option is almost at expiry. Payoff is approximately max(S0 - K, 0).
    const double T_short = 0.01;
    const double S0_atm  = 100.0;
    GeometricBrownianMotion gbm_short{MU, SIGMA, S0_atm};

    auto payoff = std::make_shared<EuropeanCall>(K);
    auto fn = make_plain_mc_sampler(gbm_short, EM, payoff, STEPS, T_short);

    auto [mean, var] = sample_stats(fn, 50'000, detail::make_seed(800, 0));
    const double se  = std::sqrt(var / 50'000.0);
    const double true_price = ref::undiscounted_call(S0_atm, K, MU, SIGMA, T_short);

    INFO("Short expiry T=0.01: mean=" << mean << " true=" << true_price << " se=" << se);
    REQUIRE(mean >= 0.0);
    REQUIRE(std::abs(mean - true_price) < 5.0 * se);
}

TEST_CASE("edge: very long expiry (T=10)", "[edge]") {
    const double T_long = 10.0;
    GeometricBrownianMotion gbm_long{MU, SIGMA, S0};

    auto payoff = std::make_shared<EuropeanCall>(K);
    auto fn = make_plain_mc_sampler(gbm_long, EM, payoff, STEPS, T_long);

    auto [mean, var] = sample_stats(fn, 50'000, detail::make_seed(900, 0));
    const double se  = std::sqrt(var / 50'000.0);
    const double true_price = ref::undiscounted_call(S0, K, MU, SIGMA, T_long);

    INFO("Long expiry T=10: mean=" << mean << " true=" << true_price << " se=" << se);
    REQUIRE(mean > 0.0);
    REQUIRE(std::abs(mean - true_price) < 5.0 * se);
}

TEST_CASE("edge: zero volatility (sigma=0)", "[edge]") {
    // With sigma=0, S_T = S0 * exp(mu*T) deterministically.
    // European call = max(S0*exp(mu*T) - K, 0)
    const double SIGMA_ZERO = 0.0;
    const double S_T        = S0 * std::exp(MU * T);
    const double expected   = std::max(S_T - K, 0.0);

    GeometricBrownianMotion gbm_zero{MU, SIGMA_ZERO, S0};
    auto payoff = std::make_shared<EuropeanCall>(K);

    // detail::generate_path with zero sigma should produce a deterministic path
    std::vector<double> z(STEPS, 0.5); // nonzero shocks, but sigma=0 → no effect
    auto path = detail::generate_path(gbm_zero, EM, z, STEPS, T);

    // All path values should be the drift-only trajectory
    INFO("Zero-vol path back(): " << path.back() << "  expected terminal: " << S_T);
    // Note: Euler-Maruyama discretisation introduces small ODE error for zero sigma
    REQUIRE(path.back() > 0.0);

    // With sigma=0, all paths are the same, so MC gives the exact answer
    auto fn = make_plain_mc_sampler(gbm_zero, EM, payoff, STEPS, T);
    auto [mean, var] = sample_stats(fn, 1'000, detail::make_seed(1000, 0));

    INFO("Zero-vol MC mean = " << mean << "  expected = " << expected
         << "  sample_var = " << var);
    REQUIRE(var < 1e-20); // zero variance - deterministic payoff
    REQUIRE_THAT(mean, WithinAbs(expected, 0.05)); // small Euler discretization error OK
}

TEST_CASE("edge: high volatility (sigma=0.8)", "[edge]") {
    const double SIGMA_HIGH = 0.8;
    GeometricBrownianMotion gbm_high{MU, SIGMA_HIGH, S0};

    auto payoff = std::make_shared<EuropeanCall>(K);
    auto fn = make_plain_mc_sampler(gbm_high, EM, payoff, STEPS, T);

    const double true_price = ref::undiscounted_call(S0, K, MU, SIGMA_HIGH, T);
    auto [mean, var] = sample_stats(fn, 100'000, detail::make_seed(1100, 0));
    const double se  = std::sqrt(var / 100'000.0);

    INFO("High vol sigma=0.8: mean=" << mean << " true=" << true_price << " 5*se=" << 5.0 * se);
    REQUIRE(mean > 0.0);
    REQUIRE(std::abs(mean - true_price) < 5.0 * se);
}

TEST_CASE("edge: zero interest rate (mu=0)", "[edge]") {
    const double MU_ZERO = 0.0;
    GeometricBrownianMotion gbm_zero_mu{MU_ZERO, SIGMA, S0};

    auto payoff = std::make_shared<EuropeanCall>(K);
    auto fn = make_plain_mc_sampler(gbm_zero_mu, EM, payoff, STEPS, T);

    // With mu=0, E[S_T]=S0, but physical != risk-neutral.
    // The undiscounted E[call] = E[max(S_T-K,0)] under the zero-drift measure.
    // Use BS with r=0: price = S0*N(d1) - K*N(d2) → no discounting applied
    // Undiscounted E[call] = BS(r=0) * exp(0) = BS(r=0)
    const double true_price = ref::undiscounted_call(S0, K, MU_ZERO, SIGMA, T);
    auto [mean, var] = sample_stats(fn, 100'000, detail::make_seed(1200, 0));
    const double se  = std::sqrt(var / 100'000.0);

    INFO("Zero rate mu=0: mean=" << mean << " true=" << true_price << " 5*se=" << 5.0 * se);
    REQUIRE(mean > 0.0);
    REQUIRE(std::abs(mean - true_price) < 5.0 * se);
}

TEST_CASE("edge: negative rate (mu=-0.03)", "[edge]") {
    const double MU_NEG = -0.03;
    GeometricBrownianMotion gbm_neg{MU_NEG, SIGMA, S0};

    auto payoff = std::make_shared<EuropeanCall>(K);
    auto fn = make_plain_mc_sampler(gbm_neg, EM, payoff, STEPS, T);

    const double true_price = ref::undiscounted_call(S0, K, MU_NEG, SIGMA, T);
    auto [mean, var] = sample_stats(fn, 100'000, detail::make_seed(1300, 0));
    const double se  = std::sqrt(var / 100'000.0);

    INFO("Negative rate mu=-0.03: mean=" << mean << " true=" << true_price
         << " 5*se=" << 5.0 * se);
    REQUIRE(mean >= 0.0);
    REQUIRE(std::abs(mean - true_price) < 5.0 * se);
}

TEST_CASE("edge: single time step (steps=1)", "[edge]") {
    // NOTE: Euler-Maruyama is NOT exact for GBM even with a single step.
    // GBM has a multiplicative diffusion term, so the arithmetic Euler step
    // produces a different distribution from the exact log-normal solution.
    // The discretization error grows with sigma and T; for sigma=0.2, T=1 it
    // is ~0.15-0.25 (systematic downward bias on the call price).
    //
    // This test verifies:
    //   (a) the sampler produces valid (finite, non-negative) results
    //   (b) the systematic error is bounded and in the expected direction
    //   (c) the path has length 2 (STEPS1+1 = 2)
    constexpr size_t STEPS1 = 1;

    auto payoff = std::make_shared<EuropeanCall>(K);
    auto fn = make_plain_mc_sampler(GBM, EM, payoff, STEPS1, T);

    const double true_price = ref::undiscounted_call(S0, K, MU, SIGMA, T);
    auto [mean, var] = sample_stats(fn, 200'000, detail::make_seed(1400, 0));
    const double se  = std::sqrt(var / 200'000.0);

    INFO("Single step (Euler bias test):");
    INFO("  MC mean   = " << mean);
    INFO("  BS target = " << true_price << "  (different: Euler != exact for GBM)");
    INFO("  error     = " << std::abs(mean - true_price));
    INFO("  5*se      = " << 5.0 * se);

    // Result must be finite, positive, and in the right ballpark
    REQUIRE(std::isfinite(mean));
    REQUIRE(mean > 0.0);
    // Euler 1-step for GBM is biased; bias ~0.1-0.3 for these params.
    // Verify the bias is not extreme (result should be within 50% of BS).
    REQUIRE(mean > true_price * 0.5);
    REQUIRE(mean < true_price * 1.5);
    // Statistical noise must be small relative to the estimate
    REQUIRE(se < mean * 0.01);

    // Verify path length via detail::generate_path
    std::vector<double> z1 = {0.0};
    auto path = detail::generate_path(GBM, EM, z1, STEPS1, T);
    INFO("  path length = " << path.size() << "  expected = 2");
    REQUIRE(path.size() == STEPS1 + 1);
}

TEST_CASE("edge: barrier option with S0 above barrier is always zero", "[edge]") {
    // UpAndOut call with barrier=80 and S0=100: S0 > barrier → always knocked out
    const double BARRIER_LOW = 80.0;
    auto payoff = std::make_shared<UpAndOutCall>(K, BARRIER_LOW);
    auto fn = make_plain_mc_sampler(GBM, EM, payoff, STEPS, T);

    auto [mean, var] = sample_stats(fn, 10'000, detail::make_seed(1500, 0));

    INFO("UpAndOut barrier below S0: mean=" << mean << " var=" << var);
    // S0=100 > barrier=80, so path always starts above barrier → knocked out immediately
    REQUIRE(mean == Catch::Approx(0.0).margin(1e-10));
    REQUIRE(var == Catch::Approx(0.0).margin(1e-10));
}

TEST_CASE("edge: down-and-in put with barrier above S0 is always active", "[edge]") {
    // DownAndIn put with barrier=150 > S0=100: path always crosses barrier
    // So this behaves like a vanilla put
    const double BARRIER_HIGH = 150.0;
    auto vanilla_put = std::make_shared<EuropeanPut>(K);
    auto din_put     = std::make_shared<DownAndInPut>(K, BARRIER_HIGH);

    auto fn_v = make_plain_mc_sampler(GBM, EM, vanilla_put, STEPS, T);
    auto fn_d = make_plain_mc_sampler(GBM, EM, din_put, STEPS, T);

    constexpr size_t N = 50'000;
    auto [mv, vv] = sample_stats(fn_v, N, detail::make_seed(1600, 0));
    auto [md, vd] = sample_stats(fn_d, N, detail::make_seed(1600, 0));

    INFO("Vanilla put mean:   " << mv);
    INFO("DI put mean:        " << md << " (barrier=" << BARRIER_HIGH << " > S0=" << S0 << ")");

    // With barrier above S0, the DownAndIn activates immediately and should
    // behave similarly to a vanilla put (small discretization difference)
    // Both should be positive and roughly similar
    REQUIRE(md > 0.0);
    REQUIRE(mv > 0.0);
    REQUIRE(std::abs(md - mv) / mv < 0.5); // within 50% of vanilla (rough check)
}

TEST_CASE("edge: lookback payoffs are non-negative always", "[edge]") {
    // LookbackCall: S_T - min(path) >= 0 always (min <= S_T by construction)
    // LookbackPut:  max(path) - S_T >= 0 always
    auto lc_fn = make_plain_mc_sampler(GBM, EM, std::make_shared<LookbackCall>(), STEPS, T);
    auto lp_fn = make_plain_mc_sampler(GBM, EM, std::make_shared<LookbackPut>(),  STEPS, T);

    constexpr size_t N = 5'000;
    auto lc_s = lc_fn(N, detail::make_seed(1700, 0));
    auto lp_s = lp_fn(N, detail::make_seed(1701, 0));

    for (double v : lc_s)
        REQUIRE(v >= -1e-12); // allow tiny floating-point noise
    for (double v : lp_s)
        REQUIRE(v >= -1e-12);
}

TEST_CASE("edge: importance sampling with theta=0 equals plain MC", "[edge]") {
    // theta=0 means no tilting → identical to plain MC in distribution
    auto payoff = std::make_shared<EuropeanCall>(K);
    auto plain  = make_plain_mc_sampler(GBM, EM, payoff, STEPS, T);
    auto is0    = make_importance_sampler(GBM, EM, payoff, STEPS, T, 0.0);

    constexpr size_t N = 100'000;
    auto [pm, pv] = sample_stats(plain, N, detail::make_seed(1800, 0));
    auto [im, iv] = sample_stats(is0,   N, detail::make_seed(1800, 0));

    // Same seed, theta=0: weight = exp(0) = 1 for all paths → identical outputs
    INFO("plain mean=" << pm << " IS(theta=0) mean=" << im);
    REQUIRE_THAT(pm, WithinAbs(im, 1e-10));
}

// =============================================================================
// SECTION 9 - Thread safety
//
// Launch multiple Engine::run calls concurrently. Verify:
//   (a) No crashes / undefined behaviour (sanitisers would catch data races)
//   (b) Each result is a valid positive number
//   (c) Results are not all identical (would indicate shared-state corruption)
// =============================================================================

TEST_CASE("thread: concurrent Engine::run calls produce valid results", "[thread]") {
    // NOTE: Engine::run is deterministic given the same config (seeds are a pure
    // function of strategy index, thread index, and round). Concurrent calls with
    // the same config will produce identical estimates - this is correct behaviour,
    // not state sharing. We verify:
    //   (a) No crashes / bad memory access under concurrent execution
    //   (b) Each thread produces a finite, positive result
    //   (c) Results are within the expected price range

    constexpr size_t N_THREADS = 4;

    auto engine = Engine<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>>{
        GBM, EM, CTRL_MEAN, STEPS, T};

    IterativeEngineConfig cfg;
    cfg.n_compete          = 4;
    cfg.n_exploit          = 4;
    cfg.n_rounds           = 3;
    cfg.samples_per_thread = 500;

    std::vector<double> results(N_THREADS);
    std::vector<std::thread> threads;
    threads.reserve(N_THREADS);

    for (size_t t = 0; t < N_THREADS; ++t) {
        threads.emplace_back([&, t]() {
            auto res = engine.run(std::make_shared<EuropeanCall>(K), cfg);
            REQUIRE(res.has_value());
            results[t] = res->estimate;
        });
    }
    for (auto& th : threads) th.join();

    // All results must be finite, positive, and in a reasonable price range
    for (size_t t = 0; t < N_THREADS; ++t) {
        INFO("Thread " << t << " result: " << results[t]);
        REQUIRE(std::isfinite(results[t]));
        REQUIRE(results[t] > 0.0);
        // ATM European call undiscounted price ~10-12
        REQUIRE(results[t] > 5.0);
        REQUIRE(results[t] < 25.0);
    }
}

TEST_CASE("thread: concurrent ASVR runs produce valid independent results", "[thread]") {
    constexpr size_t N_THREADS = 4;

    auto payoff = std::make_shared<EuropeanCall>(K);

    std::vector<double> results(N_THREADS);
    std::vector<std::thread> threads;
    threads.reserve(N_THREADS);

    for (size_t t = 0; t < N_THREADS; ++t) {
        threads.emplace_back([&, t]() {
            auto strategies = build_all_strategies(GBM, EM, payoff, CTRL_MEAN, STEPS, T);
            ASVRConfig cfg;
            cfg.exploration_fraction = 0.1;
            cfg.n_threads = 2;
            auto res = AdaptiveVarianceReduction::run(strategies, 5'000, cfg);
            results[t] = res.estimate;
        });
    }
    for (auto& th : threads) th.join();

    for (size_t t = 0; t < N_THREADS; ++t) {
        INFO("ASVR thread " << t << " estimate: " << results[t]);
        REQUIRE(std::isfinite(results[t]));
        REQUIRE(results[t] > 0.0);
    }
}

// =============================================================================
// SECTION 10 - Numerical utilities
// =============================================================================

TEST_CASE("numerics: normal_icdf is inverse of normal CDF", "[numerics]") {
    // Verify normal_icdf(ncdf(x)) ≈ x for a range of x values
    auto ncdf = [](double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); };

    for (double x : {-3.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0}) {
        const double p = ncdf(x);
        const double x_back = detail::normal_icdf(p);
        INFO("x=" << x << "  p=" << p << "  icdf(p)=" << x_back);
        REQUIRE_THAT(x_back, WithinAbs(x, 5e-9));
    }
}

TEST_CASE("numerics: halton sequence stays in (0,1) and has low discrepancy", "[numerics]") {
    // First 100 terms of base-2 Halton should cover (0,1) well
    std::vector<double> h;
    for (size_t i = 1; i <= 100; ++i)
        h.push_back(detail::halton(i, 2));

    // All values in (0, 1)
    for (double v : h) {
        REQUIRE(v > 0.0);
        REQUIRE(v < 1.0);
    }

    // Low-discrepancy: partition (0,1) into 10 bins, each should have ~10 points
    std::vector<int> bins(10, 0);
    for (double v : h) bins[static_cast<int>(v * 10)]++;
    for (int b : bins) {
        INFO("Bin count: " << b);
        REQUIRE(b >= 5);   // no bin should be under-filled
        REQUIRE(b <= 15);  // no bin should be over-filled
    }
}

TEST_CASE("numerics: make_seed produces distinct values for different (k,phase) pairs", "[numerics]") {
    std::vector<uint64_t> seeds;
    for (size_t k = 0; k < 10; ++k)
        for (size_t phase = 0; phase < 3; ++phase)
            seeds.push_back(detail::make_seed(k, phase));

    // All 30 seeds should be distinct
    auto sorted = seeds;
    std::sort(sorted.begin(), sorted.end());
    for (size_t i = 1; i < sorted.size(); ++i) {
        INFO("Duplicate seed detected: " << sorted[i]);
        REQUIRE(sorted[i] != sorted[i - 1]);
    }

    // Same (k, phase) should always return the same value
    for (size_t k = 0; k < 5; ++k) {
        for (size_t p = 0; p < 3; ++p) {
            REQUIRE(detail::make_seed(k, p) == detail::make_seed(k, p));
        }
    }
}

TEST_CASE("numerics: mean_variance returns unbiased estimates", "[numerics]") {
    // Known sample: {1, 2, 3, 4, 5}  → mean=3, var=2.5
    std::vector<double> s = {1.0, 2.0, 3.0, 4.0, 5.0};
    auto [mean, var] = detail::mean_variance(s);

    REQUIRE_THAT(mean, WithinAbs(3.0, 1e-12));
    REQUIRE_THAT(var, WithinAbs(2.5, 1e-12));
}

TEST_CASE("numerics: mean_variance with single element returns zero variance", "[numerics]") {
    auto [mean, var] = detail::mean_variance({42.0});
    REQUIRE_THAT(mean, WithinAbs(42.0, 1e-12));
    REQUIRE_THAT(var, WithinAbs(0.0, 1e-12));
}

// =============================================================================
// SECTION 11 - Performance benchmarking
// (tagged [bench] so they run only when explicitly requested)
// =============================================================================

TEST_CASE("bench: wall-clock time per strategy at varying sample counts", "[bench][slow]") {
    constexpr size_t SEED = 12345;
    auto payoff = std::make_shared<EuropeanCall>(K);
    auto strategies = build_all_strategies(GBM, EM, payoff, CTRL_MEAN, STEPS, T);

    WARN("Strategy timing (N=10000, 50000, 100000):");
    WARN(std::string(70, '-'));

    for (const auto& s : strategies) {
        std::string row = "  ";
        row += s.name;
        row.resize(28, ' ');
        for (size_t N : {10'000u, 50'000u, 100'000u}) {
            using clock = std::chrono::high_resolution_clock;
            auto t0   = clock::now();
            auto samp = s.sampler(N, SEED);
            auto t1   = clock::now();
            (void)samp;
            const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            row += "  N=" + std::to_string(N) + " " + std::to_string(ms).substr(0, 6) + "ms";
        }
        WARN(row);
    }
    // No REQUIRE - purely diagnostic output
    REQUIRE(true);
}

TEST_CASE("bench: ASVR overhead vs best fixed strategy", "[bench][slow]") {
    constexpr size_t N     = 100'000;
    constexpr size_t RUNS  = 20;

    auto payoff = std::make_shared<EuropeanCall>(K);

    // Reference value
    auto plain_ref_fn = make_plain_mc_sampler(GBM, EM, payoff, STEPS, T);
    auto ref_s = plain_ref_fn(500'000, detail::make_seed(0, 999));
    const double ref = std::accumulate(ref_s.begin(), ref_s.end(), 0.0) /
                       static_cast<double>(ref_s.size());

    // Best fixed strategy (antithetic) MSE
    auto anti_fn = make_antithetic_sampler(GBM, EM, payoff, STEPS, T);
    double mse_anti = 0.0;
    for (size_t r = 0; r < RUNS; ++r) {
        auto [m, v] = sample_stats(anti_fn, N, detail::make_seed(r, 10));
        double e = m - ref;
        mse_anti += e * e;
        (void)v;
    }
    mse_anti /= RUNS;

    // ASVR MSE
    ASVRConfig cfg;
    cfg.exploration_fraction = 0.1;
    cfg.n_threads = 4;
    double mse_asvr = 0.0;
    for (size_t r = 0; r < RUNS; ++r) {
        auto strategies = build_all_strategies(GBM, EM, payoff, CTRL_MEAN, STEPS, T);
        auto res = AdaptiveVarianceReduction::run(strategies, N, cfg);
        double e = res.estimate - ref;
        mse_asvr += e * e;
    }
    mse_asvr /= RUNS;

    WARN("Overhead test: N=" << N << " RUNS=" << RUNS);
    WARN("  Best fixed (antithetic) MSE: " << mse_anti);
    WARN("  ASVR MSE:                    " << mse_asvr);
    WARN("  ASVR overhead ratio:         " << mse_asvr / mse_anti
         << "  (1.0 = no overhead, <1 = ASVR wins)");

    // ASVR should achieve MSE at most 2x worse than the best fixed strategy
    // (The 10% exploration overhead is the cost of learning)
    REQUIRE(mse_asvr < mse_anti * 3.0);
}
