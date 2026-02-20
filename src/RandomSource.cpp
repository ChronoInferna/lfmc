#include "lfmc/RandomSource.hpp"

#include "lfmc/types.hpp"

namespace lfmc {

PseudoRandomSource::PseudoRandomSource(unsigned seed) : rng_(seed), dist_(0.0, 1.0) {}

Normals PseudoRandomSource::generate(size_t n) {
    Normals normals(n);
    for (size_t i = 0; i < n; ++i) {
        normals[i] = dist_(rng_);
    }
    return normals;
}

void PseudoRandomSource::seed(unsigned seed) {
    rng_.seed(seed);
}

} // namespace lfmc
