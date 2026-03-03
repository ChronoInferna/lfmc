#pragma once

#include "lfmc/payoff.hpp"
#include "lfmc/types.hpp"

#include <expected>
#include <numeric>
#include <algorithm>
#include <vector>
#include <string>

namespace lfmc {


    class LookbackCall : public Payoff {
public:
    std::expected<std::vector<Payoffs>, std::string>



    generate_payoffs(const std::vector<Path>& paths) const override {
        Payoffs payoffs;
        payoffs.reserve(paths.size());

        for (const auto& path : paths) {
            if (path.empty())
                return std::unexpected("Empty path encountered in LookbackCall");



            double min_price = *std::min_element(path.begin(), path.end());
            payoffs.push_back(path.back() - min_price);
        }

        return std::vector<Payoffs>{payoffs};
    }
};



class LookbackPut : public Payoff {
public:

    std::expected<std::vector<Payoffs>, std::string>
    generate_payoffs(const std::vector<Path>& paths) const override {
        Payoffs payoffs;
        payoffs.reserve(paths.size());

        for (const auto& path : paths) {
            if (path.empty())
                return std::unexpected("Empty path encountered in LookbackPut");

            double max_price = *std::max_element(path.begin(), path.end());
    
    
    
            payoffs.push_back(max_price - path.back());
        }

        return std::vector<Payoffs>{payoffs};
    }
};

}