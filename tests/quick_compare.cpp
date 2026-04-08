// Quick comparison: Engine (winner-takes-all) vs single fixed strategies
// Uses same total sample budget across all approaches.

#include "lfmc/adaptive_estimator.hpp"
#include "lfmc/engine.hpp"
#include "lfmc/numerical_scheme.hpp"
#include "lfmc/payoff.hpp"
#include "lfmc/stochastic_process.hpp"
#include "lfmc/strategies.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

using namespace lfmc;

static constexpr double S0 = 100.0;
static constexpr double MU = 0.05;
static constexpr double SIGMA = 0.2;
static constexpr double K = 100.0;
static constexpr double T = 1.0;
static constexpr size_t STEPS = 52;
static const double CONTROL_MEAN = S0 * std::exp(MU * T);

static GeometricBrownianMotion GBM{MU, SIGMA, S0};
static EulerMaruyama<GeometricBrownianMotion> EULER;

// Total budget: Engine uses n_rounds*(n_compete+n_exploit)*spt = 2*30*800 = 48000 samples.
// Fixed strategies get the same 48000 samples for a fair comparison.
static constexpr size_t TOTAL = 48000;
static constexpr size_t RUNS = 50;

static double ref_mean(std::shared_ptr<Payoff> payoff) {
    auto fn = make_plain_mc_sampler(GBM, EULER, payoff, STEPS, T);
    auto s = fn(500'000, detail::make_seed(0, 77));
    double sum = 0.0;
    for (double x : s)
        sum += x;
    return sum / static_cast<double>(s.size());
}

static double bench_fixed(const std::string& name, SamplerFn fn, double ref, size_t runs) {
    double mse = 0.0;
    for (size_t r = 0; r < runs; ++r) {
        auto s = fn(TOTAL, detail::make_seed(r + 1000, 5));
        double est = std::accumulate(s.begin(), s.end(), 0.0) / static_cast<double>(s.size());
        double e = est - ref;
        mse += e * e;
    }
    return mse / static_cast<double>(runs);
}

static double bench_engine(std::shared_ptr<Payoff> payoff, double ref, size_t runs) {
    Engine<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>> engine{
        GBM, EULER, CONTROL_MEAN, STEPS, T};

    // IterativeEngineConfig: n_rounds*(n_compete+n_exploit)*samples_per_thread = 2*30*800 = 48000
    IterativeEngineConfig cfg;
    cfg.n_compete = 10;
    cfg.n_exploit = 20;
    cfg.n_rounds = 2;
    cfg.samples_per_thread = 800;

    double mse = 0.0;
    for (size_t r = 0; r < runs; ++r) {
        cfg.run_index = r;
        auto res = engine.run(payoff, cfg);
        if (!res)
            continue;
        double e = res->estimate - ref;
        mse += e * e;
    }
    return mse / static_cast<double>(runs);
}

void compare(const char* opt_name, std::shared_ptr<Payoff> payoff) {
    printf("\n=== %s ===\n", opt_name);
    const double ref = ref_mean(payoff);
    printf("Reference: %.4f  (budget per run: %zu samples)\n\n", ref, TOTAL);

    auto strategies = build_all_strategies(GBM, EULER, payoff, CONTROL_MEAN, STEPS, T);

    printf("  %-28s  %10s  %8s\n", "Method", "MSE", "vs PlainMC");
    printf("  %-28s  %10s  %8s\n", "----------------------------", "----------", "--------");

    double plain_mse = 0.0;
    std::vector<std::pair<std::string, double>> rows;

    for (const auto& s : strategies) {
        double mse = bench_fixed(s.name, s.sampler, ref, RUNS);
        rows.push_back({s.name, mse});
        if (s.name == "plain_mc")
            plain_mse = mse;
    }

    double engine_mse = bench_engine(payoff, ref, RUNS);
    rows.push_back({"[Engine]", engine_mse});

    std::sort(rows.begin(), rows.end(), [](auto& a, auto& b) { return a.second < b.second; });

    for (auto& [name, mse] : rows) {
        double ratio = plain_mse / mse;
        printf("  %-28s  %10.6f  %8.2fx\n", name.c_str(), mse, ratio);
    }
}

int main() {
    printf("Engine vs Fixed Strategy Comparison\n");
    printf("Same total sample budget (%zu), %zu independent runs\n\n", TOTAL, RUNS);

    compare("European Call", std::make_shared<EuropeanCall>(K));
    compare("Asian Call", std::make_shared<AsianCall>(K));
    compare("Barrier (Up-Out Call B=130)", std::make_shared<UpAndOutCall>(K, 130.0));
    compare("Lookback Call", std::make_shared<LookbackCall>());

    printf("\nDone.\n");
}
