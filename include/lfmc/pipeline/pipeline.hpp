#pragma once

#include "lfmc/estimator/estimator.hpp"
#include "lfmc/path_generator/path_generator.hpp"
#include "lfmc/payoff/payoff.hpp"
#include "lfmc/random_source/random_source.hpp"

#include <expected>
#include <memory>

namespace lfmc {

template <StochasticProcess SP, NumericalScheme<SP> NS> class Pipeline {
  private:
    std::unique_ptr<RandomSource> random_source_;
    std::unique_ptr<PathGenerator<SP, NS>> path_generator_;
    std::unique_ptr<Payoff> payoff_;
    std::unique_ptr<Estimator> estimator_;

  public:
    Pipeline(std::unique_ptr<RandomSource> random_source,
             std::unique_ptr<PathGenerator<SP, NS>> path_generator, std::unique_ptr<Payoff> payoff,
             std::unique_ptr<Estimator> estimator)
        : random_source_(std::move(random_source)), path_generator_(std::move(path_generator)),
          payoff_(std::move(payoff)), estimator_(std::move(estimator)) {}

    std::expected<double, std::string> run(size_t steps, double T) {
        while (!estimator_->converged()) {
            // TODO make all the return types std::expected? IDK if it's a good idea for hot loops
            // like this but it would make error handling easier and more consistent. For now just
            // return
            auto normals = random_source_->generate_normals(steps);
            if (normals.empty()) {
                return std::unexpected("Failed to generate random normals");
            }

            auto paths = path_generator_->generate_paths(normals, steps, T);
            if (paths.empty()) {
                return std::unexpected("Failed to generate paths");
            }

            auto payoffs = payoff_->generate_payoffs(paths);
            if (!payoffs) {
                return std::unexpected("Failed to generate payoffs");
            }

            auto result = estimator_->add_payoffs(payoffs.value());
            if (!result) {
                return std::unexpected(result.error());
            }
        }

        return estimator_->result();
    }
};

} // namespace lfmc
