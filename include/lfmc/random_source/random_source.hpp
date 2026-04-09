#pragma once

#include "lfmc/core/types.hpp"

#include <vector>

namespace lfmc {

class RandomSource {
  public:
    virtual ~RandomSource() = default;
    virtual std::vector<Normals> generate_normals(size_t steps, size_t n = 1) = 0;
};

} // namespace lfmc
