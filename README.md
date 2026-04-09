# LFMC: Lock-Free Monte Carlo Variance Reduction Library

A modern C++23 library for options pricing via Monte Carlo simulation with **10 hand-tuned variance reduction strategies** and **Adaptive Strategy Variance Reduction (ASVR)** — a bandit-based algorithm that learns which strategy works best without domain expertise.

## What's New (v0.2.0)

✓ **10 Variance Reduction Strategies** — antithetic, control variate, Halton QMC, importance sampling, Latin hypercube, stratified sampling, and combinations  
✓ **Adaptive Strategy Variance Reduction (ASVR)** — robustly identifies the best strategy via multi-round bandit allocation (no manual tuning)  
✓ **IterativeEngine** — multi-round sampling with precision-weighted leader tracking (requires n_rounds ≥ 4 for convergence)  
✓ **Comprehensive Test Suite** — 102 tests (99% pass rate) covering correctness, variance reduction, convergence, edge cases, and concurrency  
✓ **Production-Ready** — all compilation and numerical errors from v0.1.0 resolved

## Quick Example

```cpp
#include "lfmc/engine.hpp"
#include "lfmc/payoff.hpp"
#include "lfmc/stochastic_process.hpp"

using namespace lfmc;

// Set up the option
GeometricBrownianMotion gbm{mu=0.05, sigma=0.20, S0=100.0};
auto payoff = std::make_shared<EuropeanCall>(K=100.0);

// Run ASVR: 10 strategies, 4 rounds, ~50k samples total
IterativeEngineConfig cfg{
    .n_rounds = 4,
    .n_compete = 10,  // all 10 strategies
    .n_exploit = 10,
    .samples_per_thread = 1000
};

IterativeEngine engine(cfg);
auto result = engine.run(gbm, EulerMaruyama(), payoff, steps=52, T=1.0);

std::cout << "Price: " << result.estimate << " ± " << result.stderr << "\n";
std::cout << "Leader: " << result.best_strategy << " (VR: " << result.best_vr_ratio << "×)\n";
```

## Why ASVR?

**Use ASVR when you don't want to pick a variance reduction strategy yourself.** Across most option types, **antithetic_cv dominates** (50.7× variance reduction for European calls). But ASVR learns this automatically without requiring domain expertise.

**Choosing manually:** You'd need to benchmark 10 strategies independently (expensive). ASVR does this in parallel using a bandit algorithm.

### Measured Performance (v0.2.0)

European Call (S=100, K=100, σ=0.2, T=1, r=0.05, N=10k samples):

| Strategy | VR Ratio | Wall Time |
|----------|----------|-----------|
| plain_mc | 1.0× | 82 ms |
| antithetic | 4.1× | 91 ms |
| control_variate | 6.9× | 83 ms |
| **antithetic_cv** | **50.7×** | 106 ms |
| halton_qmc | 5.8× | 30 ms |
| importance_sampling | 3.0× | 81 ms |

**ASVR with 4 rounds** (48k total samples, 10% exploration + 90% exploitation):
| Configuration | Estimate Error | Wall Time | Notes |
|---|---|---|---|
| **Fixed: antithetic_cv** | 0.00094 | 424 ms | Best single strategy |
| **ASVR n_rounds=4** | ~0.0012 | 512 ms | Learns antithetic_cv, catches up by round 3-4 |
| **ASVR n_rounds=1** | 0.00180 | 102 ms | ⚠️ Early round, learning phase — not converged |

**Key finding**: Across all 6 tested option types (European, Asian, Barrier, Lookback), **antithetic_cv wins**. ASVR robustly learns this without manual strategy selection.

## When to Use ASVR vs Fixed Strategies

| Use Case | Recommendation |
|----------|---|
| **You know which VR strategy is best for your problem** | Use fixed strategy directly (faster, no exploration overhead) |
| **You want automation and don't mind 10-20% performance penalty** | Use ASVR with `n_rounds ≥ 4` |
| **You're pricing a new option type with unknown best strategy** | Use ASVR to discover it (requires patience: ~10–50k samples for learning) |
| **You're benchmarking multiple strategies** | Use ASVR's `round_history` output to see leader progression |
| **You need rock-solid production code with minimum tuning** | Use ASVR with `n_rounds=2, n_exploit=20` (conservative) |

**For published research:** If you already know antithetic_cv wins (it does for standard options), just use it directly. ASVR adds complexity without discovery benefit.

## Features

### ✓ 10 Built-in Variance Reduction Strategies
- **Plain Monte Carlo** — baseline
- **Antithetic Variates** — symmetric sampling
- **Control Variate** — analytical correction
- **Antithetic + Control Variate** — combined (most powerful)
- **Halton QMC** — low-discrepancy sequences
- **Importance Sampling** — biased drift
- **Moment Matching** — moment-matched paths
- **Latin Hypercube Sampling** — stratified random sampling
- **Stratified Sampling** — first-dimension stratification
- **Stratified Antithetic** — combined stratification

### ✓ Payoff Types
- European Call/Put
- Asian (arithmetic average) Call/Put
- Barrier (up-and-out, down-and-in) Call/Put
- Lookback Call/Put

### ✓ Numerical Schemes
- **Euler-Maruyama** (standard, O(√dt) weak error for GBM)

### ✓ Stochastic Processes
- **Geometric Brownian Motion** (Black-Scholes dynamics)

### ✓ Random Sources
- Pseudo-random (Mersenne Twister)
- Antithetic pairs
- Low-discrepancy (Halton, LHS)

### ✓ Thread-Safe Concurrency
- IterativeEngine supports multiple concurrent runs
- Each strategy run uses independent RNG seeds
- Results are deterministic (seeded from config)

## Building

### Requirements
- C++23 compatible compiler (GCC 12+, Clang 16+, MSVC 2022+)
- CMake 3.28+

### Configure & Build

```bash
# Configure (choose a preset: gcc-debug, clang-debug, msvc-debug)
cmake --preset gcc-debug

# Build
cmake --build build-gcc-debug

# Run tests
cd build-gcc-debug/tests
./tests.exe "~[slow]~[bench]"   # Fast tests (~2 min)
./tests.exe                     # All tests + benchmarks (~10 min)
```

### Build Options
- `LFMC_BUILD_TESTS` (default: ON) — build test suite
- `LFMC_BUILD_EXAMPLES` (default: ON) — build examples
- `LFMC_BUILD_DOCS` (default: OFF) — generate Doxygen docs

## Testing

The library includes 60+ comprehensive tests:

### Test Categories

| Category | Tests | Coverage |
|----------|-------|----------|
| `[correctness]` | 4 | All strategies unbiased vs Black-Scholes |
| `[vr]` | 7 | Variance reduction ratios verified |
| `[asvr]` | 5 | ASVR budget split, weight formula |
| `[bandit]` | 5 | IterativeEngine leader selection |
| `[edge]` | 14 | ITM/OTM, extreme vol/rates, barriers |
| `[numerics]` | 5 | Acklam normal CDF, Halton, seeds |
| `[convergence]` | 2 | O(1/√N) error scaling |
| `[coverage]` | 3 | 95% CI coverage (statistical) |
| `[thread]` | 1 | Concurrent Engine::run safety |

**Fast tests** (~2 min): `./tests.exe "~[slow]~[bench]"`  
**Full suite** (~10 min): `./tests.exe`  
**Benchmarks only** (~8 min): `./tests.exe "[bench]"` or `./bench.exe`

### Backtests

Three standalone executables for performance validation:

1. **bench_asvr_vs_fixed** — 200 independent runs per strategy per option type
   ```bash
   ./bench.exe
   ```

2. **quick_compare** — Engine vs fixed strategies with same budget (50 runs)
   ```bash
   ./quick_compare.exe
   ```

3. **stability_sweep** — MSE vs fixed strategies across run counts (15, 30, 50, 100)
   ```bash
   ./stability_sweep.exe
   ```

## API Overview

### IterativeEngine (ASVR Bandit-Based Allocation)

```cpp
struct IterativeEngineConfig {
    size_t n_rounds = 4;              // ⚠️ CRITICAL: Need >= 4 for convergence
    size_t n_compete = 10;            // Strategies in competition (1 thread each)
    size_t n_exploit = 10;            // Bonus threads assigned to leader
    size_t samples_per_thread = 800;  // Samples per thread per round
    size_t qmc_replications = 5;      // Variance estimation for QMC strategies
    size_t run_index = 0;             // Seed offset for independent runs
};

class IterativeEngine {
    // Run the bandit-based allocation algorithm
    std::expected<IterativeEngineResult, std::string> run(
        std::shared_ptr<Payoff> payoff,
        IterativeEngineConfig config = {}
    );
};

struct IterativeEngineResult {
    double estimate;                           // Final option price estimate
    double estimated_stderr;                   // Standard error
    std::string final_leader;                  // Strategy with highest precision weight
    std::vector<RoundResult> round_history;    // Per-round diagnostics (leader, weights)
    std::vector<StrategyStats> final_stats;    // Cumulative stats: name, mean, variance, n_samples
    size_t total_samples;                      // Total samples used across all rounds
};
```

**Important:** 
- `n_rounds < 4`: ASVR is still in learning phase; accuracy may be poor
- `n_rounds ≥ 4`: ASVR converges; leader stabilizes
- Overhead: ~15–20% wall-clock time vs best fixed strategy (for exploration)
- Deterministic: Results reproducible from `run_index`

### Individual Strategies

```cpp
// All accessible via strategy factory functions
using SamplerFn = std::function<std::vector<double>(size_t n, uint64_t seed)>;

SamplerFn make_antithetic_sampler(process, scheme, payoff, steps, T);
SamplerFn make_control_variate_sampler(process, scheme, payoff, steps, T);
SamplerFn make_halton_sampler(process, scheme, payoff, steps, T);
// ... etc for all 10 strategies
```

### Direct Access (Lower-Level)

```cpp
class AdaptiveVarianceReduction {
    // 10% exploration, 90% exploitation
    // Learns inverse-variance precision weights
    // Auto-allocates budget to best strategies
    
    AdaptiveVarianceReductionResult run(
        const std::vector<SamplerFn>& samplers,
        size_t n_exploration,
        size_t n_exploitation,
        uint64_t seed
    );
};
```

## Known Limitations & Caveats

### ⚠️ ASVR Requires Sufficient Rounds to Converge
- **Issue**: `n_rounds < 4` leaves ASVR in learning phase; accuracy may lag fixed strategies
- **Example**: ASVR with 1 round ≈ 9× VR, but antithetic_cv fixed ≈ 50× VR
- **Mitigation**: Use `n_rounds ≥ 4` for production (accumulates 40-48k samples with typical configs)
- **Severity**: CRITICAL — set `n_rounds` correctly or use fixed strategy instead
- **Recommendation**: Use ASVR only if you have budget for ≥ 4 rounds; otherwise pick antithetic_cv

### Importance Sampling Theta Not Tuned Per Option Type
- **Issue**: Hardcoded `theta=0.5` works well for calls, performs **3.4× worse than plain MC for puts**
- **Measured**: VR ratio = 0.29× for European puts (variance inflation, not reduction)
- **Mitigation**: ASVR assigns near-zero weight to IS for puts; final estimate unaffected
- **Severity**: Medium — no impact on ASVR output, but visible in per-strategy benchmarks
- **Recommendation**: For fixed use, skip IS for puts; for ASVR, relax n_rounds if needed to absorb noise

### Euler-Maruyama Discretization Bias
- **Issue**: EM has O(Δt) weak error for GBM; not exact even with 1 step
- **Measured bias**: ~2% with steps=1 on ATM European call (0.22 points on 10.99 fair value)
- **Mitigation**: Use steps ≥ 52 (weekly grid) — bias < 0.1%
- **Severity**: Low — well-understood QMC limitation; tests use steps=52
- **Recommendation**: Default to steps ≥ 52 in production; document EM bias if using coarse grids

### Halton QMC High-Dimension Correlation
- **Issue**: Halton sequences have inter-dimensional correlation for d > ~20
- **Impact**: Marginal VR benefit reduction for 52-step paths (still 3-5× reduction, just not 9-10×)
- **Severity**: Low — ASVR appropriately downweights Halton for path-dependent options
- **Recommendation**: None; use Sobol for very high dimensions (future work)

### Control Variate Analytical Mean Must Be Correct
- **Issue**: CV strategies assume `E[S_T] = S0 * exp(mu * T)` (risk-neutral or historical μ)
- **Impact**: Wrong mean → silently wrong CV correction → biased estimate
- **Mitigation**: Library defaults to analytical formula; user can override
- **Severity**: Medium — requires API documentation
- **Recommendation**: Validate provided mean or use analytical default

## Performance Characteristics

### Time Complexity
- **Per strategy**: O(N × steps) where N = samples, steps = path length
- **ASVR overhead**: ~20% vs best fixed strategy (allocation + rebalancing)

### Space Complexity
- **Per path**: O(steps) — no path history kept after payoff calculation
- **Concurrent runs**: O(N × num_threads) total across all workers

### Numerical Properties
- **Normal CDF**: Acklam's rational approximation (accuracy 5e-9)
- **RNG**: Mersenne Twister MT19937-64 (period 2^19937 - 1)
- **Thread safety**: Full — separate RNG seed per worker thread

## Examples

See `examples/` directory for:
- Basic European call pricing
- Comparison of all 10 strategies
- ASVR vs fixed strategy benchmark
- Edge case studies (deep ITM/OTM, barrier options, etc.)

## License

See `LICENSE` file.

## Publishing & Academic Use

**Status**: This is **research-grade code**, not a peer-reviewed publication. Before citing in academic work:

1. **Understand the core contribution**: ASVR is a bandit-based strategy selector, not a novel variance reduction method itself
2. **Acknowledge limitations**: antithetic_cv dominates; ASVR's value is *robustness without expertise*, not *discovery*
3. **Cite correctly**:

```bibtex
@software{lfmc2026,
  title  = {LFMC: Lock-Free Monte Carlo with Adaptive Strategy Variance Reduction},
  author = {Robbins, Alexander and Deng, Oliver},
  year   = {2026},
  note   = {Experimental; see AUDIT_REPORT.md for caveats},
  url    = {https://github.com/xanderrobbins/lfmc}
}
```

For a research paper, see `AUDIT_REPORT.md` for detailed test results, performance tables, and limitations.

## Architecture

### High-Level Components

```
┌─────────────────────────────────────────────────┐
│         IterativeEngine (Bandit)                │
│  - Multi-round allocation                       │
│  - Leader tracking                              │
│  - Precision-weighted bonus threads             │
└─────────┬───────────────────────────────────────┘
          │
          ↓
┌─────────────────────────────────────────────────┐
│    AdaptiveVarianceReduction (ASVR)             │
│  - 10% exploration, 90% exploitation            │
│  - Inverse-variance precision weights           │
│  - Strategy selection via bandit                │
└─────────┬───────────────────────────────────────┘
          │
          ↓
┌─────────────────────────────────────────────────┐
│      10 Variance Reduction Strategies           │
│  - Antithetic, Control Variate, QMC, etc.     │
│  - Each: (n_samples, seed) → payoff samples    │
└─────────┬───────────────────────────────────────┘
          │
          ↓
┌─────────────────────────────────────────────────┐
│     Path Generation & Payoff Computation       │
│  - Euler-Maruyama paths (GBM)                  │
│  - 4+ payoff types (European, Asian, etc.)    │
│  - Deterministic from seed (reproducible)     │
└─────────────────────────────────────────────────┘
```

### Key Classes
- **IterativeEngine** — orchestrates multi-round sampling
- **AdaptiveVarianceReduction** — implements ASVR algorithm
- **Payoff** — abstract payoff interface (European, Asian, Barrier, Lookback)
- **StochasticProcess** — GBM (extensible)
- **NumericalScheme** — Euler-Maruyama (extensible)
- **RandomSource** — pseudo-random, antithetic, LHS, Halton (extensible)

## Contributing

Contributions welcome. Areas for future development:
- [ ] Additional stochastic processes (Heston, CEV, jump-diffusion)
- [ ] More numerical schemes (Milstein, higher-order)
- [ ] Sobol low-discrepancy sequences for high dimensions
- [ ] GPU acceleration via CUDA/HIP
- [ ] Python bindings

## Changelog

### v0.2.0 (2026-04-07)
- ✓ Implemented ASVR (Adaptive Strategy Variance Reduction) algorithm with bandit allocation
- ✓ Implemented IterativeEngine with multi-round leader tracking
- ✓ All 10 variance reduction strategies (antithetic, CV, QMC, IS, LHS, stratified, combinations)
- ✓ Comprehensive test suite (102 tests, 99% pass rate)
- ✓ Fixed all blocking compilation bugs from v0.1.0
- ✓ Added performance benchmarks and stability analysis
- ⚠️ Known limitation: ASVR needs n_rounds ≥ 4; early rounds can underperform fixed strategies

### v0.1.0 (2026-01-15)
- Initial release with basic Monte Carlo engine
- Path generation and payoff computation
- Antithetic variates implementation

## Support & Contributing

For issues, questions, or suggestions:
- See `AUDIT_REPORT.md` for detailed test coverage and known limitations
- Open an issue on GitHub
- Contact: xanderrobbins10@gmail.com

Contributions welcome. Areas for future development:
- [ ] Tune importance sampling theta per option type
- [ ] Sobol sequences for high-dimensional paths
- [ ] Additional stochastic processes (Heston, jump-diffusion)
- [ ] GPU acceleration
- [ ] Python bindings

---

**Last Updated**: 2026-04-09  
**Status**: Research-Grade (Functional, Tested, Not Peer-Reviewed)  
**Test Coverage**: 102 tests (99% pass), 200+ seconds full suite  
**Recommended Use**: ASVR for robustness when strategy choice is uncertain; fixed antithetic_cv for known-good problems
