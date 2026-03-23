#pragma once

#include "lfmc/estimator/estimator.hpp"
#include "lfmc/estimator/monte_carlo_estimator.hpp"
#include "lfmc/numerical_scheme/numerical_scheme.hpp"
#include "lfmc/path_generator/path_generator.hpp"
#include "lfmc/payoff/payoff.hpp"
#include "lfmc/pipeline/pipeline.hpp"
#include "lfmc/random_source/pseudo_random_source.hpp"
#include "lfmc/random_source/random_source.hpp"
#include "lfmc/stochastic_process/stochastic_process.hpp"

namespace lfmc {

template <StochasticProcess SP, NumericalScheme<SP> NS> class PipelineBuilder {
  private:
    std::unique_ptr<RandomSource> random_source_;
    std::unique_ptr<PathGenerator<SP, NS>> path_generator_;
    std::unique_ptr<Payoff> payoff_;
    std::unique_ptr<Estimator> estimator_;

  public:
    PipelineBuilder&
    random_source(std::unique_ptr<RandomSource> rs = std::make_unique<PseudoRandomSource>()) {
        random_source_ = std::move(rs);
        return *this;
    }

    PipelineBuilder& path_generator(std::unique_ptr<PathGenerator<SP, NS>> pg) {
        path_generator_ = std::move(pg);
        return *this;
    }

    PipelineBuilder& payoff(std::unique_ptr<Payoff> p) {
        payoff_ = std::move(p);
        return *this;
    }

    PipelineBuilder&
    estimator(std::unique_ptr<Estimator> e = std::make_unique<MonteCarloEstimator>()) {
        estimator_ = std::move(e);
        return *this;
    }

    std::expected<Pipeline<SP, NS>, std::string> build() {
        if (!random_source_)
            return std::unexpected("RandomSource missing");

        if (!path_generator_)
            return std::unexpected("PathGenerator missing");

        if (!payoff_)
            return std::unexpected("Payoff missing");

        if (!estimator_)
            return std::unexpected("Estimator missing");

        return Pipeline<SP, NS>(std::move(random_source_), std::move(path_generator_),
                                std::move(payoff_), std::move(estimator_));
    }
};

} // namespace lfmc
