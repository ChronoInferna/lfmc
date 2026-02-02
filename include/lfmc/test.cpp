#define _USE_MATH_DEFINES  
#include <cmath>
#include "StochasticProcess.hpp"
#include "NumericalScheme.hpp"
#include "Payoff.hpp"
#include "Simulator.hpp"
#include <iostream>
#include <iomanip>

//for comp 
double black_scholes_call(double S, double K, double r, double sigma, double T) {
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    
    // NormCDF approx
    auto norm_cdf = [](double x) {
        return 0.5 * std::erfc(-x * M_SQRT1_2);
    };
    
    return S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
}

int main() {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Monte Carlo Option Pricing\n\n";
    
    // Parameters
    double S0 = 100.0;       // Initial stock price
    double K = 100.0;        // Strike price
    double r = 0.05;         // Risk-free rate (5%)
    double sigma = 0.20;     // Volatility (20%)
    double T = 1.0;          // Time to maturity (1 year)
    size_t n_steps = 252;    // Daily steps
    size_t n_sims = 1000;   // Number of simulations
    
    lfmc::GeometricBrownianMotion gbm(r, sigma);
    lfmc::EulerMaruyama<lfmc::GeometricBrownianMotion> euler(gbm);
    
    // Create payoff (European Call)
    lfmc::EuropeanCall call(K);
    
    // Create simulator
    lfmc::Simulator sim(gbm, euler, call, S0, T, n_steps, n_sims);
    
    // Run sims
    std::cout << "Running " << n_sims << " Monte Carlo simulations\n";
    auto [mc_price, std_error] = sim.run_simulations_with_error();
    
    // Discount to present value
    double discounted_price = mc_price * std::exp(-r * T);
    double discounted_error = std_error * std::exp(-r * T);
    
    // Calculate Black-Scholes for comparison
    double bs_price = black_scholes_call(S0, K, r, sigma, T);
    
    // Results
    std::cout << "\n--- Results ---\n";
    std::cout << "Monte Carlo Price:    " << discounted_price << " + or - " << discounted_error << "\n";
    std::cout << "Black-Scholes Price:  " << bs_price << "\n";
    std::cout << "Absolute Error:       " << std::abs(discounted_price - bs_price) << "\n";
    std::cout << "Relative Error:       " << std::abs(discounted_price - bs_price) / bs_price * 100 << "%\n";
    
    // 95% confidence interval
    double ci_lower = discounted_price - 1.96 * discounted_error;
    double ci_upper = discounted_price + 1.96 * discounted_error;
    std::cout << "95% CI:               [" << ci_lower << ", " << ci_upper << "]\n";
    
    if (bs_price >= ci_lower && bs_price <= ci_upper) {
        std::cout << "Black-Scholes price is within confidence interval\n";
    } else {
        std::cout << "Black-Scholes price is outside confidence interval\n";
    }
    
    return 0;
}