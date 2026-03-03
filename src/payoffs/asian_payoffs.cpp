#pragma once

#include "lfmc/payoff.hpp"
#include "lfmc/types.hpp"

#include <expected>
#include <numeric>
#include <algorithm>
#include <vector>
#include <string>

namespace lfmc {




class AsianCall : public Payoff {
private:
    double strike_;

public:
    explicit AsianCall(double strike) : strike_(strike) {}

    std::expected<std::vector<Payoffs>, std::string>
    generate_payoffs(const std::vector<Path>& paths) const override {
        Payoffs payoffs;
        
        payoffs.reserve(paths.size());

        for (const auto& path : paths) {
            
            if (path.empty())
                return std::unexpected("Empty path encountered in AsianCall");

            double mean = std::reduce(path.begin(), path.end(), 0.0) 
                          / static_cast<double>(path.size());
            payoffs.push_back(std::max(mean - strike_, 0.0));
        }

        return std::vector<Payoffs>{payoffs};
    }
};


class AsianPut : public Payoff {
private:
    double strike_;

public:
    explicit AsianPut(double strike) : strike_(strike) {}

    std::expected<std::vector<Payoffs>, std::string>
    generate_payoffs(const std::vector<Path>& paths) const override {
        
        
        Payoffs payoffs;
        
        
        
        payoffs.reserve(paths.size());

        for (const auto& path : paths) {
            if (path.empty())
                return std::unexpected("Empty path encountered in AsianPut");

            
            double mean = std::reduce(path.begin(), path.end(), 0.0)
                          / static_cast<double>(path.size());
            payoffs.push_back(std::max(strike_ - mean, 0.0));
        }

        return std::vector<Payoffs>{payoffs};
    }
};

}