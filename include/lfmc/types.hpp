#pragma once

#include <vector>

namespace lfmc {

struct State {
    double initialValue;
    double timeToMaturity;
    size_t steps;
};

using Path = std::vector<double>;
using Normals = std::vector<double>;

} // namespace lfmc
