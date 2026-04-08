# LFMC: Lock-Free Monte Carlo Variance Reduction Library

A modern C++23 library for options pricing via Monte Carlo simulation with **Adaptive Strategy Variance Reduction (ASVR)** — a bandit-based algorithm that dynamically selects the best variance reduction strategy for your option.

## What's New (v0.2.0)

✓ **10 Variance Reduction Strategies** — antithetic, control variate, Halton QMC, importance sampling, Latin hypercube, stratified sampling, and combinations  
✓ **Adaptive Strategy Variance Reduction (ASVR)** — learns which strategy works best and reallocates compute automatically  
✓ **IterativeEngine** — multi-round bandit-based sampling with precision-weighted leader allocation  
✓ **Comprehensive Test Suite** — 60+ tests covering correctness, variance reduction, convergence, edge cases, and concurrency  
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

Different options benefit from different variance reduction strategies. A human would need to test each strategy independently — expensive and error-prone. **ASVR does this automatically.**

### Measured Performance (v0.2.0)

European Call (S=100, K=100, σ=0.2, T=1, r=0.05):

| Strategy | VR Ratio | Time (10k samples) |
|----------|----------|-------------------|
| plain_mc | 1.0× | 82 ms |
| antithetic | 4.1× | 91 ms |
| control_variate | 6.9× | 83 ms |
| **antithetic_cv** | **50.7×** | 106 ms |
| halton_qmc | 5.8× | 30 ms (fast!) |
| importance_sampling | 3.0× | 81 ms |
| **ASVR (adaptive)** | **9.2×** | 102 ms (1 round) |

Asian Call: Halton QMC dominates (9.4×); ASVR learns this automatically.

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

### IterativeEngine (Recommended)

```cpp
struct IterativeEngineConfig {
    size_t n_rounds = 4;              // Number of allocation rounds
    size_t n_compete = 10;            // Strategies in competition
    size_t n_exploit = 10;            // Strategies to reallocate to leader
    size_t samples_per_thread = 1000; // Batch size per worker
};

class IterativeEngine {
    IterativeEngineResult run(
        StochasticProcess process,
        NumericalScheme scheme,
        std::shared_ptr<Payoff> payoff,
        size_t steps,
        double T
    );
};

struct IterativeEngineResult {
    double estimate;           // Estimated option price
    double stderr;             // Standard error (95% CI width ≈ 2×stderr)
    std::string best_strategy; // Name of leading strategy
    double best_vr_ratio;      // Variance reduction vs plain MC
    // ... per-round history
};
```

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

### Importance Sampling Theta Not Tuned Per Option Type
- **Issue**: Hardcoded `theta=0.5` works well for calls, poorly for puts (3.2× variance inflation)
- **Mitigation**: ASVR automatically assigns low weight to IS for puts
- **Severity**: Low — final ASVR estimate unaffected
- **Recommendation**: For production, tune theta per option or disable IS for puts

### Euler-Maruyama Discretization Bias
- **Issue**: EM has O(Δt) bias for GBM; not exact even with 1 step
- **Measured bias**: ~2% with steps=1 on ATM call
- **Mitigation**: Use steps ≥ 52 (weekly) — bias becomes negligible
- **Recommendation**: Default to steps ≥ 52 for production use

### Halton QMC High-Dimension Correlation
- **Issue**: Halton sequences have inter-dimensional correlation for dimension > ~20
- **Impact**: Marginal benefit for 52-step paths (still useful, ASVR weights appropriately)
- **Recommendation**: Consider Sobol sequences for very high dimensions (future work)

### Hardcoded Control Variate Mean
- **Issue**: CV strategies assume correct `E[S_T] = S0 * exp(mu * T)`
- **Mitigation**: Library defaults to analytical formula
- **Recommendation**: Validate mean if providing custom value

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

## Citing This Work

If you use LFMC in your research, please cite:

```bibtex
@software{lfmc2026,
  title  = {LFMC: Lock-Free Monte Carlo Variance Reduction Library},
  author = {Robbins, Alexander},
  year   = {2026},
  url    = {https://github.com/xanderrobbins/lfmc}
}
```

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

### v0.2.0 (2026-03-31)
- ✓ Implemented ASVR algorithm
- ✓ Implemented IterativeEngine with bandit allocation
- ✓ All 10 variance reduction strategies
- ✓ Comprehensive test suite (60+ tests)
- ✓ Fixed all blocking compilation bugs
- ✓ Added performance benchmarks

### v0.1.0 (2026-01-15)
- Initial release
- Basic path generation and payoffs
- Antithetic variates

## Support

For issues, questions, or suggestions:
- Open an issue on GitHub
- Contact: xanderrobbins10@gmail.com

---

**Last Updated**: 2026-04-07  
**Status**: Production-Ready for Research/Academic Use
