#include "lfmc/adaptive/adaptive_estimator.hpp"
#include "lfmc/engine/engine.hpp"
#include "lfmc/numerical_scheme/euler_maruyama.hpp"
#include "lfmc/payoff/asian_payoffs.hpp"
#include "lfmc/payoff/barrier_payoffs.hpp"
#include "lfmc/payoff/european_payoffs.hpp"
#include "lfmc/payoff/lookback_payoffs.hpp"
#include "lfmc/stochastic_process/geometric_brownian_motion.hpp"
#include "lfmc/strategy/strategies.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numeric>

using namespace lfmc;

static constexpr double S0 = 100.0;
static constexpr double MU = 0.05;
static constexpr double SIGMA = 0.2;
static constexpr double K = 100.0;
static constexpr double T = 1.0;
static constexpr size_t STEPS = 52;
static const double CONTROL_MEAN = S0 * std::exp(MU * T);

static Engine<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>> make_engine() {
    return {GeometricBrownianMotion{MU, SIGMA, S0}, EulerMaruyama<GeometricBrownianMotion>{},
            CONTROL_MEAN, STEPS, T};
}

TEST_CASE("Engine: returns a valid positive estimate for European call") {
    auto engine = make_engine();
    auto result = engine.run(std::make_shared<EuropeanCall>(K));

    REQUIRE(result.has_value());
    REQUIRE(result->estimate > 0.0);
    REQUIRE(result->estimated_stderr > 0.0);
}

TEST_CASE("Engine: estimate is within reasonable range of Black-Scholes") {
    auto engine = make_engine();
    auto result = engine.run(std::make_shared<EuropeanCall>(K));

    REQUIRE(result.has_value());
    WARN("Estimate: " << result->estimate << " ± " << result->estimated_stderr);
    WARN("Leader:   " << result->final_leader);
    REQUIRE(std::abs(result->estimate - 10.99) < 2.0);
}

TEST_CASE("Engine: round history has n_rounds entries") {
    auto engine = make_engine();

    IterativeEngineConfig cfg;
    cfg.n_rounds = 3;

    auto result = engine.run(std::make_shared<EuropeanCall>(K), cfg);

    REQUIRE(result.has_value());
    REQUIRE(result->round_history.size() == 3);
    for (size_t i = 0; i < 3; ++i)
        REQUIRE(result->round_history[i].round == i + 1);
}

TEST_CASE("Engine: precision weights sum to 1 after each round") {
    auto engine = make_engine();

    IterativeEngineConfig cfg;
    cfg.n_rounds = 4;

    auto result = engine.run(std::make_shared<EuropeanCall>(K), cfg);
    REQUIRE(result.has_value());

    for (const auto& r : result->round_history) {
        double sum = 0.0;
        for (double w : r.precision_weights)
            sum += w;
        REQUIRE(sum == Catch::Approx(1.0).epsilon(1e-10));
    }
}

TEST_CASE("Engine: total samples consistent with config") {
    auto engine = make_engine();

    IterativeEngineConfig cfg;
    cfg.n_compete = 4;
    cfg.n_exploit = 8;
    cfg.n_rounds = 3;
    cfg.samples_per_thread = 500;

    auto result = engine.run(std::make_shared<EuropeanCall>(K), cfg);
    REQUIRE(result.has_value());

    // Round 1: (4 compete + 8 exploit spread equally) = 4*(1+2)=12 threads * 500
    //          + leftover 8%4=0, so 3 per strategy = 12 threads * 500 = 6000
    // Rounds 2-3: 4 compete (1 each) + 8 exploit on leader = 12 threads * 500 = 6000 each
    // Total = 3 * 12 * 500 = 18000
    const size_t expected = cfg.n_rounds * (cfg.n_compete + cfg.n_exploit) * cfg.samples_per_thread;
    REQUIRE(result->total_samples == expected);
}

TEST_CASE("Engine: final stats cover all competing strategies") {
    auto engine = make_engine();
    auto result = engine.run(std::make_shared<EuropeanCall>(K));

    REQUIRE(result.has_value());
    REQUIRE(result->final_stats.size() == 10); // default n_compete = 10 (all strategies)
    for (const auto& s : result->final_stats) {
        REQUIRE(s.n_samples > 0);
        REQUIRE(s.sample_variance >= 0.0);
    }
}

TEST_CASE("Engine: works with Asian call payoff") {
    auto engine = make_engine();
    auto result = engine.run(std::make_shared<AsianCall>(K));

    REQUIRE(result.has_value());
    REQUIRE(result->estimate > 0.0);
}

TEST_CASE("Engine: works with barrier option payoff") {
    auto engine = make_engine();
    auto result = engine.run(std::make_shared<UpAndOutCall>(K, 130.0));

    REQUIRE(result.has_value());
    REQUIRE(result->estimate >= 0.0);
}

TEST_CASE("Engine: leader stabilises over rounds") {
    auto engine = make_engine();

    IterativeEngineConfig cfg;
    cfg.n_compete = 4;
    cfg.n_exploit = 8;
    cfg.n_rounds = 6;
    cfg.samples_per_thread = 1000;

    auto result = engine.run(std::make_shared<EuropeanCall>(K), cfg);
    REQUIRE(result.has_value());

    // Count leader changes in the last half of rounds
    size_t late_changes = 0;
    const auto& h = result->round_history;
    for (size_t i = h.size() / 2; i < h.size(); ++i)
        if (h[i].leader_changed)
            ++late_changes;

    WARN("Leader per round:");
    for (const auto& r : h)
        WARN("  Round " << r.round << ": " << r.leader << (r.leader_changed ? " (SWITCHED)" : ""));

    // Allow at most 1 late switch - the leader might flip once as estimates refine
    REQUIRE(late_changes <= 1);
}

TEST_CASE("Engine: empirical MSE vs fixed strategies", "[bench]") {
    static constexpr size_t TOTAL = 48000;
    static constexpr size_t RUNS = 10;

    GeometricBrownianMotion gbm{MU, SIGMA, S0};
    EulerMaruyama<GeometricBrownianMotion> euler;

    struct Row {
        std::string name;
        double mse;
    };

    auto run_comparison = [&](const char* opt_name, std::shared_ptr<Payoff> payoff) {
        auto plain_ref = make_plain_mc_sampler(gbm, euler, payoff, STEPS, T);
        auto ref_s = plain_ref(300'000, detail::make_seed(0, 77));
        double ref =
            std::accumulate(ref_s.begin(), ref_s.end(), 0.0) / static_cast<double>(ref_s.size());

        auto strategies = build_all_strategies(gbm, euler, payoff, CONTROL_MEAN, STEPS, T);

        std::vector<Row> rows;

        for (const auto& s : strategies) {
            double mse = 0.0;
            for (size_t r = 0; r < RUNS; ++r) {
                auto samps = s.sampler(TOTAL, detail::make_seed(r + 1000, 5));
                double est = std::accumulate(samps.begin(), samps.end(), 0.0) /
                             static_cast<double>(samps.size());
                double e = est - ref;
                mse += e * e;
            }
            rows.push_back({s.name, mse / RUNS});
        }

        // Engine: 10 compete + 20 exploit, 2 rounds, 800 samples/thread
        // = 2 * (10 + 20) * 800 = 48000 total
        Engine<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>> engine{
            gbm, euler, CONTROL_MEAN, STEPS, T};
        IterativeEngineConfig cfg;
        cfg.n_compete = 10;
        cfg.n_exploit = 20;
        cfg.n_rounds = 2;
        cfg.samples_per_thread = 800;
        {
            double mse = 0.0;
            for (size_t r = 0; r < RUNS; ++r) {
                cfg.run_index = r;
                auto res = engine.run(payoff, cfg);
                REQUIRE(res.has_value());
                double e = res->estimate - ref;
                mse += e * e;
            }
            rows.push_back({"[Engine]", mse / RUNS});
        }

        std::sort(rows.begin(), rows.end(), [](auto& a, auto& b) { return a.mse < b.mse; });
        double plain_mse = 0.0;
        for (auto& r : rows)
            if (r.name == "plain_mc")
                plain_mse = r.mse;

        WARN(opt_name << "  (ref=" << ref << ", budget=" << TOTAL << ", runs=" << RUNS << ")");
        for (auto& r : rows)
            WARN("  " << r.name << ": MSE=" << r.mse << "  (" << (plain_mse / r.mse)
                      << "x vs plain_mc)");

        double engine_mse = 0.0;
        for (auto& r : rows)
            if (r.name == "[Engine]")
                engine_mse = r.mse;
        REQUIRE(engine_mse < plain_mse);
    };

    run_comparison("European Call", std::make_shared<EuropeanCall>(K));
    run_comparison("Asian Call", std::make_shared<AsianCall>(K));
    run_comparison("Barrier (Up-Out B=130)", std::make_shared<UpAndOutCall>(K, 130.0));
    run_comparison("Lookback Call", std::make_shared<LookbackCall>());
}
