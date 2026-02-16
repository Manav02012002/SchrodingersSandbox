#include "core/hydrogen.h"

#include "core/constants.h"
#include "core/special_functions.h"

#include <algorithm>
#include <cmath>

namespace sbox::hydrogen {

namespace {

constexpr double kMinR = 1e-15;
constexpr double kRhoCutoff = 100.0;

}  // namespace

double radial_wavefunction(int n, int l, double zeff, double r) {
    if (n <= 0 || l < 0 || l >= n || zeff <= 0.0 || r < 0.0) {
        return 0.0;
    }

    const double radial_scale = 2.0 * zeff / (static_cast<double>(n) * a0);
    const double log_n = 0.5 *
                         (3.0 * std::log(radial_scale) +
                          std::lgamma(static_cast<double>(n - l)) -
                          std::log(2.0 * static_cast<double>(n)) -
                          std::lgamma(static_cast<double>(n + l + 1)));
    const double normalization = std::exp(log_n);

    if (r < kMinR) {
        return l == 0 ? normalization : 0.0;
    }

    const double rho = radial_scale * r;
    if (rho > kRhoCutoff) {
        return 0.0;
    }

    const int laguerre_n = n - l - 1;
    const double laguerre_alpha = static_cast<double>(2 * l + 1);
    const double laguerre = sbox::math::associated_laguerre(laguerre_n, laguerre_alpha, rho);

    return normalization * std::exp(-0.5 * rho) * std::pow(rho, static_cast<double>(l)) * laguerre;
}

double orbital_value(int n, int l, int m, double zeff, double x, double y, double z) {
    const double r = std::sqrt(x * x + y * y + z * z);
    if (r <= 0.0) {
        return 0.0;
    }

    const double cos_theta = std::clamp(z / r, -1.0, 1.0);
    const double theta = std::acos(cos_theta);
    const double phi = std::atan2(y, x);

    const double radial = radial_wavefunction(n, l, zeff, r);
    const double angular = sbox::math::real_spherical_harmonic(l, m, theta, phi);
    return radial * angular;
}

double probability_density(int n, int l, int m, double zeff, double x, double y, double z) {
    const double psi = orbital_value(n, l, m, zeff, x, y, z);
    return psi * psi;
}

}  // namespace sbox::hydrogen
