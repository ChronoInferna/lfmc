// Stability sweep: Engine MSE vs fixed strategies across multiple run counts.
// For each run count in {15, 30, 50, 100} and each option type, reports:
//   - Engine ratio vs plain_mc
//   - Best fixed strategy name and ratio
//   - Gap between engine and best fixed
// Flags any option type whose engine ratio varies >20% between Runs=15 and Runs=100.

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
#include <string>
#include <vector>

using namespace lfmc;

static constexpr double S0    = 100.0;
static constexpr double MU    = 0.05;
static constexpr double SIGMA = 0.2;
static constexpr double K     = 100.0;
static constexpr double T     = 1.0;
static constexpr size_t STEPS = 52;
static constexpr size_t TOTAL = 48000; // budget per run
static const double CONTROL_MEAN = S0 * std::exp(MU * T);

static GeometricBrownianMotion GBM{MU, SIGMA, S0};
static EulerMaruyama<GeometricBrownianMotion> EULER;

static double ref_mean(std::shared_ptr<Payoff> payoff) {
    auto fn = make_plain_mc_sampler(GBM, EULER, payoff, STEPS, T);
    auto s = fn(1'000'000, detail::make_seed(0, 77));
    double sum = 0.0;
    for (double x : s) sum += x;
    return sum / static_cast<double>(s.size());
}

struct RunResult {
    double engine_mse;
    double plain_mse;
    double best_fixed_mse;
    std::string best_fixed_name;
};

static RunResult run_sweep(std::shared_ptr<Payoff> payoff, double ref, size_t n_runs) {
    // Fixed strategies
    auto strategies = build_all_strategies(GBM, EULER, payoff, CONTROL_MEAN, STEPS, T);

    double plain_mse = 0.0;
    double best_mse  = 1e300;
    std::string best_name;

    for (const auto& s : strategies) {
        double mse = 0.0;
        for (size_t r = 0; r < n_runs; ++r) {
            auto samps = s.sampler(TOTAL, detail::make_seed(r + 1000, 5));
            double est = std::accumulate(samps.begin(), samps.end(), 0.0) /
                         static_cast<double>(samps.size());
            double e = est - ref;
            mse += e * e;
        }
        mse /= static_cast<double>(n_runs);

        if (s.name == "plain_mc") plain_mse = mse;
        if (mse < best_mse) {
            best_mse  = mse;
            best_name = s.name;
        }
    }

    // Engine: 2 rounds * (10 compete + 20 exploit) * 800 spt = 48000
    Engine<GeometricBrownianMotion, EulerMaruyama<GeometricBrownianMotion>> engine{
        GBM, EULER, CONTROL_MEAN, STEPS, T};
    IterativeEngineConfig cfg;
    cfg.n_compete          = 10;
    cfg.n_exploit          = 20;
    cfg.n_rounds           = 2;
    cfg.samples_per_thread = 800;

    double engine_mse = 0.0;
    for (size_t r = 0; r < n_runs; ++r) {
        cfg.run_index = r;
        auto res = engine.run(payoff, cfg);
        if (!res) continue;
        double e = res->estimate - ref;
        engine_mse += e * e;
    }
    engine_mse /= static_cast<double>(n_runs);

    return {engine_mse, plain_mse, best_mse, best_name};
}

// Abbreviate long strategy names for the table
static std::string abbrev(const std::string& name) {
    if (name == "plain_mc")              return "pmc";
    if (name == "antithetic")            return "ant";
    if (name == "control_variate")       return "cv";
    if (name == "antithetic_cv")         return "acv";
    if (name == "stratified")            return "str";
    if (name == "halton_qmc")            return "hal";
    if (name == "importance_sampling")   return "imp";
    if (name == "moment_matching")       return "mm";
    if (name == "lhs")                   return "lhs";
    if (name == "stratified_antithetic") return "sat";
    return name;
}

struct OptionResult {
    const char*              name;
    std::vector<RunResult>   results; // one per run_count
};

int main() {
    static constexpr size_t RUN_COUNTS[] = {15, 30, 50, 100};
    static constexpr size_t N_COUNTS = sizeof(RUN_COUNTS) / sizeof(RUN_COUNTS[0]);

    struct OptSpec { const char* name; std::shared_ptr<Payoff> payoff; };
    std::vector<OptSpec> opts = {
        {"European Call", std::make_shared<EuropeanCall>(K)},
        {"Asian Call",    std::make_shared<AsianCall>(K)},
        {"Barrier",       std::make_shared<UpAndOutCall>(K, 130.0)},
        {"Lookback",      std::make_shared<LookbackCall>()},
    };

    printf("Computing reference means (1M samples each)...\n");
    std::vector<double> refs;
    for (auto& o : opts) {
        refs.push_back(ref_mean(o.payoff));
        printf("  %-16s ref = %.6f\n", o.name, refs.back());
    }
    printf("\n");

    // Gather all results
    std::vector<OptionResult> all_results;
    for (size_t oi = 0; oi < opts.size(); ++oi) {
        printf("Running %s...\n", opts[oi].name);
        OptionResult or_;
        or_.name = opts[oi].name;
        for (size_t ci = 0; ci < N_COUNTS; ++ci) {
            printf("  runs=%zu ...\n", RUN_COUNTS[ci]);
            fflush(stdout);
            or_.results.push_back(run_sweep(opts[oi].payoff, refs[oi], RUN_COUNTS[ci]));
        }
        all_results.push_back(std::move(or_));
    }

    // --- Print stability table ---------------------------------------------------
    printf("\n");
    printf("Stability Analysis (budget=%zu per run)\n", TOTAL);
    printf("%-16s", "Option");
    for (size_t ci = 0; ci < N_COUNTS; ++ci)
        printf("  Runs=%-3zu              ", RUN_COUNTS[ci]);
    printf("\n");

    printf("%-16s", "");
    for (size_t ci = 0; ci < N_COUNTS; ++ci)
        printf("  %-9s %-12s", "Eng", "Best");
    printf("\n");

    printf("%-16s", "");
    for (size_t ci = 0; ci < N_COUNTS; ++ci)
        printf("  %-9s %-12s", "---------", "------------");
    printf("\n");

    for (auto& or_ : all_results) {
        printf("%-16s", or_.name);
        for (size_t ci = 0; ci < N_COUNTS; ++ci) {
            const auto& r = or_.results[ci];
            double eng_ratio  = r.plain_mse / r.engine_mse;
            double best_ratio = r.plain_mse / r.best_fixed_mse;
            char eng_buf[16], best_buf[24];
            snprintf(eng_buf,  sizeof(eng_buf),  "%.1fx", eng_ratio);
            snprintf(best_buf, sizeof(best_buf), "%s:%.1fx",
                     abbrev(r.best_fixed_name).c_str(), best_ratio);
            printf("  %-9s %-12s", eng_buf, best_buf);
        }
        printf("\n");
    }

    // --- Stability flags -------------------------------------------------------
    printf("\n--- Stability flags (engine ratio, >20%% variance between Runs=15 and Runs=100) ---\n");
    for (auto& or_ : all_results) {
        const auto& r15  = or_.results[0];
        const auto& r100 = or_.results[N_COUNTS - 1];
        double ratio15  = r15.plain_mse  / r15.engine_mse;
        double ratio100 = r100.plain_mse / r100.engine_mse;
        double pct_diff = std::abs(ratio100 - ratio15) / ((ratio15 + ratio100) / 2.0) * 100.0;
        const char* flag = (pct_diff > 20.0) ? "UNSTABLE" : "stable";
        printf("  %-16s  %.1fx → %.1fx  (%+.1f%%)  %s\n",
               or_.name, ratio15, ratio100, ratio100 - ratio15, flag);
    }

    // --- Per-run-count detail (engine MSE and best fixed MSE) -----------------
    printf("\n--- Detailed MSE values ---\n");
    for (auto& or_ : all_results) {
        printf("\n%s\n", or_.name);
        printf("  %-8s  %12s  %12s  %12s  %8s  %-16s  %8s\n",
               "Runs", "Engine MSE", "PlainMC MSE", "BestFix MSE", "Eng/plain", "BestFix", "Gap");
        for (size_t ci = 0; ci < N_COUNTS; ++ci) {
            const auto& r = or_.results[ci];
            double eng_ratio  = r.plain_mse / r.engine_mse;
            double best_ratio = r.plain_mse / r.best_fixed_mse;
            // gap: how much better best-fixed is vs engine (positive = engine loses)
            double gap = best_ratio - eng_ratio;
            printf("  %-8zu  %12.6f  %12.6f  %12.6f  %8.2fx  %-16s  %+.2fx\n",
                   RUN_COUNTS[ci],
                   r.engine_mse, r.plain_mse, r.best_fixed_mse,
                   eng_ratio,
                   (abbrev(r.best_fixed_name) + ":" +
                    std::to_string(best_ratio).substr(0, 4) + "x").c_str(),
                   gap);
        }
    }

    // --- Analysis questions ---------------------------------------------------
    printf("\n=== Analysis ===\n");

    // Q1: Are Runs=15 numbers representative?
    printf("\nQ1: Are Runs=15 numbers representative of true engine performance?\n");
    for (auto& or_ : all_results) {
        double r15  = or_.results[0].plain_mse / or_.results[0].engine_mse;
        double r100 = or_.results[N_COUNTS-1].plain_mse / or_.results[N_COUNTS-1].engine_mse;
        double pct = (r15 - r100) / r100 * 100.0;
        printf("  %-16s  15-run=%.2fx  100-run=%.2fx  bias=%+.1f%%\n",
               or_.name, r15, r100, pct);
    }

    // Q2: Variance in engine MSE across run counts
    printf("\nQ2: Engine ratio variance across run counts (min/max/range):\n");
    for (auto& or_ : all_results) {
        double mn = 1e300, mx = -1e300;
        for (auto& r : or_.results) {
            double ratio = r.plain_mse / r.engine_mse;
            mn = std::min(mn, ratio);
            mx = std::max(mx, ratio);
        }
        printf("  %-16s  min=%.2fx  max=%.2fx  range=%.2fx\n",
               or_.name, mn, mx, mx - mn);
    }

    // Q3: Asian Call - is 1.1x consistent or does it sometimes look better?
    printf("\nQ3: Asian Call engine ratio across run counts:\n");
    for (auto& or_ : all_results) {
        if (std::string(or_.name).find("Asian") == std::string::npos) continue;
        for (size_t ci = 0; ci < N_COUNTS; ++ci) {
            double ratio = or_.results[ci].plain_mse / or_.results[ci].engine_mse;
            printf("  Runs=%-3zu  %.2fx\n", RUN_COUNTS[ci], ratio);
        }
    }

    // Q4: Minimum run count for stable estimates
    printf("\nQ4: Minimum run count for stable estimates (engine ratio within 10%% of Runs=100):\n");
    for (auto& or_ : all_results) {
        double r100 = or_.results[N_COUNTS-1].plain_mse / or_.results[N_COUNTS-1].engine_mse;
        for (size_t ci = 0; ci < N_COUNTS; ++ci) {
            double r = or_.results[ci].plain_mse / or_.results[ci].engine_mse;
            double pct = std::abs(r - r100) / r100 * 100.0;
            if (pct <= 10.0) {
                printf("  %-16s  Runs=%-3zu  (%.1f%% from Runs=100 estimate)\n",
                       or_.name, RUN_COUNTS[ci], pct);
                break;
            }
        }
    }

    printf("\nDone.\n");
    return 0;
}
