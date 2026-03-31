#pragma once

// All 6 new variance reduction strategy factories.
// Each factory returns a SamplerFn: (n_samples, seed) -> vector<double>.
// The seed guarantees RNG independence between strategies and ASVR phases.

#include "lfmc/adaptive_estimator.hpp"
#include "lfmc/numerical_scheme.hpp"
#include "lfmc/payoff.hpp"
#include "lfmc/stochastic_process.hpp"
#include "lfmc/types.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

namespace lfmc::detail {

// ─── Acklam's inverse normal CDF ─────────────────────────────────────────────
// Rational approximation accurate to ~5e-9 over (0, 1).
// Reference: Peter J. Acklam, "An algorithm for computing the inverse normal
// cumulative distribution function", 2010.
inline double normal_icdf(double p) noexcept {
    static constexpr double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
                                   -2.759285104469687e+02, 1.383577518672690e+02,
                                   -3.066479806614716e+01, 2.506628277459239e+00};
    static constexpr double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
                                   -1.556989798598866e+02, 6.680131188771972e+01,
                                   -1.328068155288572e+01};
    static constexpr double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                                   -2.400758277161838e+00, -2.549732539343734e+00,
                                   4.374664141464968e+00,  2.938163982698783e+00};
    static constexpr double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
                                   2.445134137142996e+00, 3.754408661907416e+00};

    if (p < 0.02425) {
        const double q = std::sqrt(-2.0 * std::log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    if (p <= 0.97575) {
        const double q = p - 0.5;
        const double r = q * q;
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
               (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    }
    const double q = std::sqrt(-2.0 * std::log(1.0 - p));
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
}

// ─── Halton radical inverse ───────────────────────────────────────────────────
// Returns the i-th term (1-indexed) of the Van der Corput sequence in `base`.
// Output in (0, 1). Caller must pass i >= 1.
inline double halton(size_t i, int base) noexcept {
    double result = 0.0;
    double f = 1.0 / static_cast<double>(base);
    size_t idx = i;
    while (idx > 0) {
        result += f * static_cast<double>(idx % static_cast<size_t>(base));
        idx /= static_cast<size_t>(base);
        f /= static_cast<double>(base);
    }
    return result;
}

// ─── First n primes (sieve of Eratosthenes) ───────────────────────────────────
// Used by Halton QMC: dimension d gets the d-th prime as its base.
inline std::vector<int> first_n_primes(size_t n) {
    std::vector<int> primes;
    primes.reserve(n);
    for (int candidate = 2; primes.size() < n; ++candidate) {
        bool is_prime = true;
        for (int p : primes) {
            if (p * p > candidate)
                break;
            if (candidate % p == 0) {
                is_prime = false;
                break;
            }
        }
        if (is_prime)
            primes.push_back(candidate);
    }
    return primes;
}

} // namespace lfmc::detail

namespace lfmc {

// ─── Strategy 5: Stratified sampling ─────────────────────────────────────────
// Stratifies only the first Brownian increment dimension (the one with highest
// correlation to the terminal price). Remaining dimensions are pseudo-random.
// The n strata are [(i/n, (i+1)/n)] for i=0..n-1; one uniform is drawn from
// each stratum and converted to N(0,1) via normal_icdf.
template <StochasticProcess SP, NumericalScheme<SP> NS>
SamplerFn make_stratified_sampler(SP process, NS scheme, std::shared_ptr<Payoff> payoff,
                                   size_t steps, double T) {
    return [process, scheme, payoff, steps, T](size_t n, uint64_t seed) -> std::vector<double> {
        std::mt19937_64 rng{seed};
        std::uniform_real_distribution<double> udist{0.0, 1.0};
        std::normal_distribution<double> ndist{0.0, 1.0};

        std::vector<double> samples;
        samples.reserve(n);

        for (size_t i = 0; i < n; ++i) {
            // First dimension: one draw from stratum i
            const double u = (static_cast<double>(i) + udist(rng)) / static_cast<double>(n);
            const double z0 = detail::normal_icdf(u);

            // Remaining dimensions: pseudo-random
            Normals normals(steps);
            normals[0] = z0;
            for (size_t j = 1; j < steps; ++j)
                normals[j] = ndist(rng);

            auto path = detail::generate_path(process, scheme, normals, steps, T);
            auto result = payoff->generate_payoffs({path});
            if (result && !result->empty() && !(*result)[0].empty())
                samples.push_back((*result)[0][0]);
        }
        return samples;
    };
}

// ─── Strategy 6: Halton QMC with random shift ─────────────────────────────────
// Generates a `steps`-dimensional Halton sequence. Each dimension is shifted
// by a random uniform (derived from the seed) before converting to N(0,1).
// The random shift randomises the sequence per ASVR phase while preserving
// the low-discrepancy property within each call.
//
// Known limitation: Halton sequences have inter-dimensional correlation for
// large dimension counts (steps > ~20). ASVR will automatically assign low
// weight to this strategy for path-dependent options where all steps matter.
template <StochasticProcess SP, NumericalScheme<SP> NS>
SamplerFn make_halton_sampler(SP process, NS scheme, std::shared_ptr<Payoff> payoff,
                               size_t steps, double T) {
    return [process, scheme, payoff, steps, T](size_t n, uint64_t seed) -> std::vector<double> {
        std::mt19937_64 rng{seed};
        std::uniform_real_distribution<double> udist{0.0, 1.0};

        // Random shift per dimension (Owen-style randomisation simplified to a shift)
        const auto primes = detail::first_n_primes(steps);
        std::vector<double> shifts(steps);
        for (double& s : shifts)
            s = udist(rng);

        std::vector<double> samples;
        samples.reserve(n);

        for (size_t i = 0; i < n; ++i) {
            Normals normals(steps);
            for (size_t d = 0; d < steps; ++d) {
                // Clamp away from the boundary to avoid normal_icdf singularity
                double u = std::fmod(detail::halton(i + 1, primes[d]) + shifts[d], 1.0);
                u = std::clamp(u, 1e-10, 1.0 - 1e-10);
                normals[d] = detail::normal_icdf(u);
            }
            auto path = detail::generate_path(process, scheme, normals, steps, T);
            auto result = payoff->generate_payoffs({path});
            if (result && !result->empty() && !(*result)[0].empty())
                samples.push_back((*result)[0][0]);
        }
        return samples;
    };
}

// ─── Strategy 7: Importance sampling (exponential tilting) ───────────────────
// Tilts all normal draws using a TOTAL drift parameter `theta`, distributed
// evenly as theta/sqrt(steps) per time step. This parameterisation makes theta
// dimensionally consistent regardless of the number of steps:
//
//   Z_t ~ N(theta/sqrt(steps), 1)  under importance measure Q
//   per-step shift: delta = theta / sqrt(steps)
//
// The Radon-Nikodym correction weight is:
//   w = exp(-delta * sum_t(Z_t) + 0.5 * delta^2 * steps)
//     = exp(-theta/sqrt(steps) * sum_t(Z_t) + 0.5 * theta^2)
//
// This is equivalent to shifting the terminal log-price by sigma*theta*sqrt(dt)*steps
// = sigma*theta*sqrt(T). It matches the single-normal exact simulation where shifting
// Z by theta yields the same log-price change.
//
// Choosing theta:
//   theta = 0        → plain MC
//   theta > 0        → tilts S_T upward; beneficial for OTM calls
//   theta < 0        → tilts S_T downward; beneficial for OTM puts / barrier knock-ins
//   practical range  → |theta| in [0, 2]; beyond 2 the Radon-Nikodym weight
//                      exp(-0.5*theta^2) shrinks too fast and variance increases
template <StochasticProcess SP, NumericalScheme<SP> NS>
SamplerFn make_importance_sampler(SP process, NS scheme, std::shared_ptr<Payoff> payoff,
                                   size_t steps, double T, double theta) {
    return [process, scheme, payoff, steps, T, theta](size_t n,
                                                       uint64_t seed) -> std::vector<double> {
        std::mt19937_64 rng{seed};
        // Per-step shift: distribute the total theta equally across sqrt(steps) normal draws
        const double delta = theta / std::sqrt(static_cast<double>(steps));
        std::normal_distribution<double> ndist{delta, 1.0}; // tilted draw

        // Log-weight constant: 0.5 * delta^2 * steps = 0.5 * theta^2
        const double lw_const = 0.5 * theta * theta;

        std::vector<double> samples;
        samples.reserve(n);

        for (size_t i = 0; i < n; ++i) {
            Normals normals(steps);
            double z_sum = 0.0;
            for (size_t j = 0; j < steps; ++j) {
                normals[j] = ndist(rng);
                z_sum += normals[j];
            }
            // Likelihood ratio dP/dQ = exp(-delta * sum(Z') + 0.5 * delta^2 * steps)
            //                        = exp(-theta/sqrt(steps) * sum(Z') + 0.5 * theta^2)
            const double weight = std::exp(-delta * z_sum + lw_const);

            auto path = detail::generate_path(process, scheme, normals, steps, T);
            auto result = payoff->generate_payoffs({path});
            if (result && !result->empty() && !(*result)[0].empty())
                samples.push_back((*result)[0][0] * weight);
        }
        return samples;
    };
}

// ─── Strategy 8: Moment matching ─────────────────────────────────────────────
// Generates all n * steps normals up front (step-major layout for cache
// efficiency in the rescaling pass), then rescales each time step's n draws
// to have exact sample mean = 0 and sample std = 1. This eliminates first-
// order sampling error in every dimension simultaneously.
//
// Memory: O(n * steps) doubles. For n=10,000 and steps=52 this is ~4 MB.
template <StochasticProcess SP, NumericalScheme<SP> NS>
SamplerFn make_moment_matching_sampler(SP process, NS scheme, std::shared_ptr<Payoff> payoff,
                                        size_t steps, double T) {
    return [process, scheme, payoff, steps, T](size_t n, uint64_t seed) -> std::vector<double> {
        std::mt19937_64 rng{seed};
        std::normal_distribution<double> ndist{0.0, 1.0};

        // Step-major storage: z[j][i] = step j, path i
        std::vector<std::vector<double>> z(steps, std::vector<double>(n));
        for (size_t j = 0; j < steps; ++j)
            for (size_t i = 0; i < n; ++i)
                z[j][i] = ndist(rng);

        // Rescale each step independently to mean=0, std=1
        for (size_t j = 0; j < steps; ++j) {
            const double nd = static_cast<double>(n);
            double mean = 0.0;
            for (double v : z[j])
                mean += v;
            mean /= nd;

            double var = 0.0;
            for (double v : z[j]) {
                double d = v - mean;
                var += d * d;
            }
            var /= (nd - 1.0);
            const double sd = std::sqrt(var);

            if (sd > 1e-12) {
                for (double& v : z[j])
                    v = (v - mean) / sd;
            }
        }

        // Simulate paths from the moment-matched normals
        std::vector<double> samples;
        samples.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            Normals normals(steps);
            for (size_t j = 0; j < steps; ++j)
                normals[j] = z[j][i];
            auto path = detail::generate_path(process, scheme, normals, steps, T);
            auto result = payoff->generate_payoffs({path});
            if (result && !result->empty() && !(*result)[0].empty())
                samples.push_back((*result)[0][0]);
        }
        return samples;
    };
}

// ─── Strategy 9: Latin Hypercube Sampling ────────────────────────────────────
// In each of the `steps` dimensions independently: generates n stratified
// uniform samples (one per stratum), randomly permutes them (so different
// dimensions are uncorrelated), then converts to N(0,1) via normal_icdf.
// Guarantees exact coverage of every marginal distribution.
template <StochasticProcess SP, NumericalScheme<SP> NS>
SamplerFn make_lhs_sampler(SP process, NS scheme, std::shared_ptr<Payoff> payoff,
                             size_t steps, double T) {
    return [process, scheme, payoff, steps, T](size_t n, uint64_t seed) -> std::vector<double> {
        std::mt19937_64 rng{seed};
        std::uniform_real_distribution<double> udist{0.0, 1.0};

        const double nd = static_cast<double>(n);

        // Build the LHS matrix: z[j][i] = normal for path i, step j
        // Column (step) j is a shuffled stratified sample
        std::vector<std::vector<double>> z(steps, std::vector<double>(n));
        for (size_t j = 0; j < steps; ++j) {
            // Stratified uniforms
            for (size_t i = 0; i < n; ++i) {
                double u = (static_cast<double>(i) + udist(rng)) / nd;
                u = std::clamp(u, 1e-10, 1.0 - 1e-10);
                z[j][i] = detail::normal_icdf(u);
            }
            // Fisher-Yates shuffle to decorrelate dimensions
            for (size_t i = n - 1; i > 0; --i) {
                std::uniform_int_distribution<size_t> idist{0, i};
                std::swap(z[j][i], z[j][idist(rng)]);
            }
        }

        std::vector<double> samples;
        samples.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            Normals normals(steps);
            for (size_t j = 0; j < steps; ++j)
                normals[j] = z[j][i];
            auto path = detail::generate_path(process, scheme, normals, steps, T);
            auto result = payoff->generate_payoffs({path});
            if (result && !result->empty() && !(*result)[0].empty())
                samples.push_back((*result)[0][0]);
        }
        return samples;
    };
}

// ─── Strategy 10: Stratified + Antithetic ────────────────────────────────────
// Combines stratified sampling on the first dimension with antithetic variates.
// For n requested samples, generates n pairs:
//   - Stratum i: u_i = (i + U) / n → z_1 = normal_icdf(u_i)
//   - Positive path: (z_1, z_2, ..., z_steps) with z_{2..steps} ~ N(0,1)
//   - Negative path: (-z_1, z_2, ..., z_steps) with the SAME remaining draws
//   - Returned sample i = 0.5 * (payoff(pos_path) + payoff(neg_path))
//
// Sharing z_{2..steps} between both paths of a pair gives variance reduction
// not only from the first-dimension antithetic pairing but also from the
// correlated higher dimensions within the pair.
template <StochasticProcess SP, NumericalScheme<SP> NS>
SamplerFn make_stratified_antithetic_sampler(SP process, NS scheme, std::shared_ptr<Payoff> payoff,
                                              size_t steps, double T) {
    return [process, scheme, payoff, steps, T](size_t n, uint64_t seed) -> std::vector<double> {
        std::mt19937_64 rng{seed};
        std::uniform_real_distribution<double> udist{0.0, 1.0};
        std::normal_distribution<double> ndist{0.0, 1.0};

        std::vector<double> samples;
        samples.reserve(n);

        for (size_t i = 0; i < n; ++i) {
            // Stratified first dimension
            const double u = (static_cast<double>(i) + udist(rng)) / static_cast<double>(n);
            const double z0 = detail::normal_icdf(u);

            // Shared higher dimensions
            Normals normals_pos(steps), normals_neg(steps);
            normals_pos[0] = z0;
            normals_neg[0] = -z0;
            for (size_t j = 1; j < steps; ++j) {
                const double zj = ndist(rng);
                normals_pos[j] = zj;
                normals_neg[j] = zj; // shared (not antithetic in remaining dims)
            }

            auto path_pos = detail::generate_path(process, scheme, normals_pos, steps, T);
            auto path_neg = detail::generate_path(process, scheme, normals_neg, steps, T);

            auto r_pos = payoff->generate_payoffs({path_pos});
            auto r_neg = payoff->generate_payoffs({path_neg});

            if (r_pos && r_neg && !r_pos->empty() && !r_neg->empty() &&
                !(*r_pos)[0].empty() && !(*r_neg)[0].empty())
                samples.push_back(0.5 * ((*r_pos)[0][0] + (*r_neg)[0][0]));
        }
        return samples;
    };
}

// ─── Convenience: build all 10 standard strategies ───────────────────────────
// Assembles the canonical 10-strategy vector for a given payoff on GBM +
// Euler-Maruyama. The control variate strategies use EuropeanCall(0.0) as the
// control (= S_T) with control_mean = S0 * exp(mu * T).
// For path-dependent options where S_T is a poor control, the CV strategies
// will receive low precision weights from ASVR automatically.
template <StochasticProcess SP, NumericalScheme<SP> NS>
std::vector<std::pair<std::string, SamplerFn>>
build_all_strategies(SP process, NS scheme, std::shared_ptr<Payoff> payoff, double control_mean,
                     size_t steps, double T, double is_theta = 0.5) {
    // The control payoff is S_T: EuropeanCall(0) gives max(S_T - 0, 0) = S_T for GBM.
    auto control_payoff = []() { return std::make_shared<EuropeanCall>(0.0); };

    return {
        {"plain_mc", make_plain_mc_sampler(process, scheme, payoff, steps, T)},
        {"antithetic", make_antithetic_sampler(process, scheme, payoff, steps, T)},
        {"control_variate",
         make_control_variate_sampler(process, scheme, payoff, control_payoff(), control_mean,
                                       steps, T)},
        {"antithetic_cv",
         make_antithetic_cv_sampler(process, scheme, payoff, control_payoff(), control_mean, steps,
                                     T)},
        {"stratified", make_stratified_sampler(process, scheme, payoff, steps, T)},
        {"halton_qmc", make_halton_sampler(process, scheme, payoff, steps, T)},
        {"importance_sampling",
         make_importance_sampler(process, scheme, payoff, steps, T, is_theta)},
        {"moment_matching", make_moment_matching_sampler(process, scheme, payoff, steps, T)},
        {"lhs", make_lhs_sampler(process, scheme, payoff, steps, T)},
        {"stratified_antithetic",
         make_stratified_antithetic_sampler(process, scheme, payoff, steps, T)},
    };
}

} // namespace lfmc
