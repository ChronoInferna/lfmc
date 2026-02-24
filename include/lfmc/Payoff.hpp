#pragma once

#include <concepts>
#include <vector>

namespace lfmc {

template <typename P>
concept TerminalPayoff = requires(P const& p, double terminal) {
    { p(terminal) } -> std::convertible_to<double>;
};

template <typename P>
concept PathPayoff = requires(P const& p, std::vector<double> const& path) {
    { p(path) } -> std::convertible_to<double>;
};

} // namespace lfmc
