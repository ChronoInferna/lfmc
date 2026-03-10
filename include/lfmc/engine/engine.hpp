#pragma once

#include "lfmc/pipeline/pipeline.hpp"
#include "lfmc/pipeline/pipeline_builder.hpp"

namespace lfmc {

// class Engine {
//   private:
//     Pipeline<StochasticProcess auto, NumericalScheme<StochasticProcess auto>> pipeline_;
//
//   public:
//     Engine(Pipeline<StochasticProcess auto, NumericalScheme<StochasticProcess auto>> pipeline)
//         : pipeline_(std::move(pipeline)) {}
//
//     std::expected<double, std::string> run() {
//         return pipeline_.run();
//     }
// };

} // namespace lfmc
