#pragma once

#include <random>
#include <vector>

namespace lfmc {

class RandomGenerator {
  public:
    explicit RandomGenerator(unsigned seed = std::random_device{}())
        : rng_(seed), normalDist_(0.0, 1.0) {}

    std::vector<double> generateNormals(size_t n) {
        std::vector<double> randoms(n);
        for (auto& r : randoms)
            r = normalDist_(rng_);
        return randoms;
    }

    double generateSingle() {
        return normalDist_(rng_);
    }

    void seed(unsigned s) {
        rng_.seed(s);
    }

  private:
    std::mt19937 rng_;
    std::normal_distribution<double> normalDist_;
};

} // namespace lfmc
