#pragma once

namespace lfmc {

// TODO maybe not the best way to do this, but the easiest and cleanest for now
enum class Strategy { PseudoRandom, Antithetic, ControlVariate };

// Helper function to convert strategy enum to string for logging and metrics
inline const char* to_string(Strategy strategy) {
    switch (strategy) {
    case Strategy::PseudoRandom:
        return "PseudoRandom";
    case Strategy::Antithetic:
        return "Antithetic";
    case Strategy::ControlVariate:
        return "ControlVariate";
    default:
        return "Unknown";
    }
}

} // namespace lfmc
