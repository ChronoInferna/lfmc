#include "lfmc/random_source.hpp"
#include "lfmc/types.hpp"

#include <expected>

namespace lfmc {

PseudoRandomSource::PseudoRandomSource(unsigned seed) : rng_(seed), dist_(0.0, 1.0) {}

std::expected<std::vector<Normals>, std::string>
PseudoRandomSource::generate_normals(size_t steps, size_t samples) {
    std::vector<Normals> result(samples, Normals(steps));
    for (size_t i = 0; i < samples; ++i) {
        Normals normals(steps);
        for (size_t j = 0; j < steps; ++j) {
            normals[j] = dist_(rng_);
        }
        result[i] = std::move(normals);
    }
    return result;
}

void PseudoRandomSource::seed(unsigned seed) {
    rng_.seed(seed);
}

} // namespace lfmc
