#pragma once
#include <random>
#include <vector>

namespace lfmc {

class RandomGenerator {
public:
    explicit RandomGenerator(unsigned seed = std::random_device{}())
        : rng_(seed), normal_dist_(0.0, 1.0) {}
    

    std::vector<double> generate_normals(size_t n) {
        std::vector<double> randoms(n);
        for (auto& r : randoms) {
            r = normal_dist_(rng_);
        }
        return randoms;
    }

    double generate_single() {
        return normal_dist_(rng_);
    }
    
    void seed(unsigned s) {
        rng_.seed(s);
    }

private:
    std::mt19937 rng_;
    std::normal_distribution<double> normal_dist_;
};

} // namespace lfmc