#pragma once

#include "lfmc/random_source/random_source.hpp"

#include <random>

namespace lfmc {

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
