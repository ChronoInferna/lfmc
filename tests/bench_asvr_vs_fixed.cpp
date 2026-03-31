// bench_asvr_vs_fixed.cpp
//
// Critical publishability benchmark:
// Compare ASVR against EVERY fixed strategy (not just plain MC)
// using the same total sample budget.
//
// If ASVR beats the *best fixed strategy known in advance*, that is a strong result.
// If ASVR only beats plain MC but not antithetic+CV, that is a weak result.
//
// Run standalone: compile and link with lfmc, then ./bench

#include "lfmc/adaptive_estimator.hpp"
#include "lfmc/numerical_scheme.hpp"
#include "lfmc/payoff.hpp"
#include "lfmc/stochastic_process.hpp"
#include "lfmc/strategies.hpp"

#include <cmath>
#include <cstdio>
#include <memory>
#include <numeric>
#include <string>
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

// Run RUNS trials of each fixed strategy and ASVR; report average squared error
struct BenchResult {
    std::string name;
    double mean_sq_error;
    double mean_estimate;
    double std_error_of_mean;
};

BenchResult bench_fixed(const std::string& name, SamplerFn sampler, double ref_mean, size_t N,
                        size_t RUNS) {
    std::vector<double> errors;
    errors.reserve(RUNS);
    double sum_est = 0.0;
    for (size_t r = 0; r < RUNS; ++r) {
        auto samples = sampler(N, detail::make_seed(r + 9999, 44));
        double est = std::accumulate(samples.begin(), samples.end(), 0.0) /
                     static_cast<double>(samples.size());
        double e = est - ref_mean;
        errors.push_back(e * e);
        sum_est += est;
    }
    double mse = std::accumulate(errors.begin(), errors.end(), 0.0) / static_cast<double>(RUNS);
    // std error of the MSE estimate (rough)
    double mean_err2 = mse;
    double sq_sum = 0.0;
    for (double e2 : errors) {
        double d = e2 - mean_err2;
        sq_sum += d * d;
    }
    double sem = std::sqrt(sq_sum / static_cast<double>(RUNS * (RUNS - 1)));
    return {name, mse, sum_est / RUNS, sem};
}

BenchResult bench_asvr(std::shared_ptr<Payoff> payoff, double ref_mean, size_t N, size_t RUNS) {
    std::vector<double> errors;
    errors.reserve(RUNS);
    double sum_est = 0.0;
    ASVRConfig cfg;
    cfg.exploration_fraction = 0.1;
    cfg.n_threads = 4;
    for (size_t r = 0; r < RUNS; ++r) {
        auto strategies = build_all_strategies(GBM, EULER, payoff, CONTROL_MEAN, STEPS, T);
        auto res = AdaptiveVarianceReduction::run(std::move(strategies), N, cfg);
        double e = res.estimate - ref_mean;
        errors.push_back(e * e);
        sum_est += res.estimate;
    }
    double mse = std::accumulate(errors.begin(), errors.end(), 0.0) / static_cast<double>(RUNS);
    double mean_err2 = mse;
    double sq_sum = 0.0;
    for (double e2 : errors) {
        double d = e2 - mean_err2;
        sq_sum += d * d;
    }
    double sem = std::sqrt(sq_sum / static_cast<double>(RUNS * (RUNS - 1)));
    return {"ASVR", mse, sum_est / RUNS, sem};
}

void run_option_bench(const std::string& opt_name, std::shared_ptr<Payoff> payoff, size_t N,
                      size_t RUNS) {
    printf("\n=== %-35s (N=%zu, runs=%zu) ===\n", opt_name.c_str(), N, RUNS);

    // Reference from 500k plain MC samples
    auto plain_ref_fn = make_plain_mc_sampler(GBM, EULER, payoff, STEPS, T);
    auto ref_samples = plain_ref_fn(500'000, detail::make_seed(0, 99));
    double ref_mean = std::accumulate(ref_samples.begin(), ref_samples.end(), 0.0) /
                      static_cast<double>(ref_samples.size());
    printf("  Reference mean: %.6f\n", ref_mean);

    auto strategies = build_all_strategies(GBM, EULER, payoff, CONTROL_MEAN, STEPS, T);
    std::vector<BenchResult> results;

    // Benchmark each fixed strategy
    for (auto& [name, sampler] : strategies) {
        results.push_back(bench_fixed(name, sampler, ref_mean, N, RUNS));
    }

    // Benchmark ASVR
    results.push_back(bench_asvr(payoff, ref_mean, N, RUNS));

    // Sort by MSE ascending
    std::sort(results.begin(), results.end(), [](const BenchResult& a, const BenchResult& b) {
        return a.mean_sq_error < b.mean_sq_error;
    });

    double plain_mse = 0.0;
    for (auto& r : results)
        if (r.name == "plain_mc")
            plain_mse = r.mean_sq_error;

    printf("  %-28s  %12s  %8s  %8s\n", "Strategy", "MSE", "vs PlainMC", "vs ASVR");
    printf("  %-28s  %12s  %8s  %8s\n", std::string(28, '-').c_str(), std::string(12, '-').c_str(),
           std::string(8, '-').c_str(), std::string(8, '-').c_str());

    double asvr_mse = 0.0;
    for (auto& r : results)
        if (r.name == "ASVR")
            asvr_mse = r.mean_sq_error;

    for (auto& r : results) {
        double vs_plain = plain_mse / r.mean_sq_error;
        double vs_asvr = asvr_mse / r.mean_sq_error;
        printf("  %-28s  %12.6f  %8.3fx  %8.3fx\n", r.name.c_str(), r.mean_sq_error, vs_plain,
               vs_asvr);
    }
}

int main() {
    const size_t N = 10'000; // sample budget per trial
    const size_t RUNS = 200; // trials for MSE estimate

    printf("ASVR vs Fixed Strategy Benchmark\n");
    printf("=================================\n");
    printf("Budget N=%zu per trial, %zu independent runs\n", N, RUNS);
    printf("Strategies compared: all 10 fixed + ASVR (adaptive)\n");

    run_option_bench("European Call", std::make_shared<EuropeanCall>(K), N, RUNS);
    run_option_bench("European Put", std::make_shared<EuropeanPut>(K), N, RUNS);
    run_option_bench("Asian Call", std::make_shared<AsianCall>(K), N, RUNS);
    run_option_bench("Asian Put", std::make_shared<AsianPut>(K), N, RUNS);
    run_option_bench("Up-And-Out Call B=130", std::make_shared<UpAndOutCall>(K, 130.0), N, RUNS);
    run_option_bench("Down-And-In Put B=70", std::make_shared<DownAndInPut>(K, 70.0), N, RUNS);
    run_option_bench("Lookback Call", std::make_shared<LookbackCall>(), N, RUNS);
    run_option_bench("Lookback Put", std::make_shared<LookbackPut>(), N, RUNS);

    printf("\nDone.\n");
    return 0;
}
