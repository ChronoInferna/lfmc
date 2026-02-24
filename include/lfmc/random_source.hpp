#pragma once

#include "lfmc/types.hpp"

#include <random>
#include <vector>

namespace lfmc {

class RandomSource {
  public:
    virtual ~RandomSource() = default;
    virtual std::vector<Normals> generate_normals(size_t steps, size_t n = 1) = 0;
};

class PseudoRandomSource : public RandomSource {
  private:
    std::mt19937 rng_;
    std::normal_distribution<double> dist_;

  public:
    PseudoRandomSource(unsigned seed = std::random_device{}());

    std::vector<Normals> generate_normals(size_t steps, size_t) override;

    void seed(unsigned seed);
};

} // namespace lfmc
