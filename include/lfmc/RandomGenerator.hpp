#pragma once

#include "types.hpp"

#include <random>

namespace lfmc {

template <class RNG>
concept RandomGenerator = requires(RNG generator, size_t n) {
    { generator.generate(n) } -> std::same_as<Normals>;
};

struct PseudoRandom {
    std::mt19937 rng_;
    std::normal_distribution<double> dist_;

    PseudoRandom(unsigned seed = std::random_device{}()) : rng_(seed), dist_(0.0, 1.0) {}

    Normals generate(size_t n) {
        Normals normals(n);
        for (size_t i = 0; i < n; ++i) {
            normals[i] = dist_(rng_);
        }
        return normals;
    }

    void seed(unsigned s) {
        rng_.seed(s);
    }
};

} // namespace lfmc
