#include "core/covalent_radii.h"

#include <array>

namespace sbox::chem {
namespace {

constexpr double kAngstromToBohr = 1.8897259886;
constexpr double kFallbackRadiusBohr = 1.5;

constexpr std::array<double, 119> kCovalentRadiiAngstrom = {{
    0.0,
    0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,
    1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
    2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22,
    1.22, 1.20, 1.19, 1.20, 1.20, 1.16,
    2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44,
    1.42, 1.39, 1.39, 1.38, 1.39, 1.40,
    2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87,
    1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50,
    2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69, 1.68, 1.68, 1.65, 1.67, 1.73, 1.76, 1.61,
    1.57, 1.49, 1.43, 1.41, 1.34, 1.29, 1.28, 1.21, 1.22, 1.36, 1.43, 1.62, 1.75, 1.65, 1.57,
}};

}  // namespace

double covalent_radius(int Z) {
    if (Z < 1 || Z >= static_cast<int>(kCovalentRadiiAngstrom.size())) {
        return kFallbackRadiusBohr;
    }

    return kCovalentRadiiAngstrom[static_cast<std::size_t>(Z)] * kAngstromToBohr;
}

}  // namespace sbox::chem
