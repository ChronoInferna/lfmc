#pragma once

#include "lfmc/core/types.hpp"

#include <expected>
#include <string>
#include <vector>

namespace lfmc {

class RandomSource {
  public:
    virtual ~RandomSource() = default;
    virtual std::expected<std::vector<Normals>, std::string> generate_normals(size_t steps,
                                                                              size_t n = 1) = 0;
};

} // namespace lfmc
