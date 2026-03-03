#include "lfmc/random_source.hpp"

namespace lfmc {
AntitheticRandomSource::AntitheticRandomSource(unsigned seed) : rng_(seed), dist_(0.0, 1.0) {}

std::expected<std::vector<Normals>, std::string>
AntitheticRandomSource::generate_normals(size_t steps, size_t samples) {
    std::vector<Normals> result(samples, Normals(steps));
    for (size_t i = 0; i < samples; ++i) {
        Normals normals(steps);
        for (size_t j = 0; j < steps; ++j) {
            double z = dist_(rng_);
            normals[j] = toggle_ ? -z : z;
        }
        result[i] = std::move(normals);
        toggle_ = !toggle_; // Alternate between normal and antithetic
    }
    return result;
}

void AntitheticRandomSource::seed(unsigned seed) {
    rng_.seed(seed);
    toggle_ = false; // Reset toggle when reseeding
}

} // namespace lfmc
