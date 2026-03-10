#pragma once

#include "lfmc/core/types.hpp"
#include "lfmc/random_source/random_source.hpp"

#include <expected>
#include <random>
#include <vector>

namespace lfmc {

class PseudoRandomSource : public RandomSource {
  private:
    std::mt19937 rng_;
    std::normal_distribution<double> dist_;

  public:
    PseudoRandomSource(unsigned seed = std::random_device{}());

    std::expected<std::vector<Normals>, std::string> generate_normals(size_t steps,
                                                                      size_t) override;

    void seed(unsigned seed);
};

} // namespace lfmc
