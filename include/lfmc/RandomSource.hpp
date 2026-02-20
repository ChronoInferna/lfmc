#pragma once

#include "lfmc/types.hpp"

#include <random>

namespace lfmc {

class RandomSource {
  public:
    virtual ~RandomSource() = default;
    virtual Normals generate(size_t n) = 0;
};

class PseudoRandomSource : public RandomSource {
  private:
    std::mt19937 rng_;
    std::normal_distribution<double> dist_;

  public:
    PseudoRandomSource(unsigned seed = std::random_device{}());

    Normals generate(size_t n) override;

    void seed(unsigned seed);
};

} // namespace lfmc
