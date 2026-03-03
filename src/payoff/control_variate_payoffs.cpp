#include "lfmc/payoff.hpp"

#include <expected>
#include <vector>

namespace lfmc {

ControlVariatePayoff::ControlVariatePayoff(std::unique_ptr<Payoff> target_payoff,
                                           std::unique_ptr<Payoff> control_payoff)
    : target_payoff_(std::move(target_payoff)), control_payoff_(std::move(control_payoff)) {}

std::expected<std::vector<Payoffs>, std::string>
ControlVariatePayoff::generate_payoffs(const std::vector<Path>& paths) const {
    auto target_payoffs = target_payoff_->generate_payoffs(paths);
    auto control_payoffs = control_payoff_->generate_payoffs(paths);

    if (!target_payoffs) {
        return std::unexpected("Failed to generate target payoffs: " + target_payoffs.error());
    } else if (!control_payoffs) {
        return std::unexpected("Failed to generate control payoffs: " + control_payoffs.error());
    }

    std::vector<Payoffs> result;
    for (size_t i = 0; i < target_payoffs.value().size(); ++i) {
        if (target_payoffs.value()[i].size() != 1 || control_payoffs.value()[i].size() != 1) {
            return std::unexpected("Each row of target payoffs must have exactly one element");
        }

        Payoffs combined_row;
        combined_row.push_back(target_payoffs.value()[i][0]);  // Original payoff
        combined_row.push_back(control_payoffs.value()[i][0]); // Control variate payoff
        result.push_back(std::move(combined_row));
    }

    return result;
}

} // namespace lfmc
