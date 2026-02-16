#include "core/hydrogen.h"

#include <cmath>
#include <functional>

#include <gtest/gtest.h>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

double simpson_integral(const std::function<double(double)>& f, double a, double b, int steps) {
    if (steps % 2 != 0) {
        ++steps;
    }

    const double h = (b - a) / static_cast<double>(steps);
    double sum = f(a) + f(b);

    for (int i = 1; i < steps; ++i) {
        const double x = a + static_cast<double>(i) * h;
        sum += (i % 2 == 0 ? 2.0 : 4.0) * f(x);
    }

    return sum * h / 3.0;
}

}  // namespace

TEST(HydrogenTest, RadialWavefunctionValues) {
    EXPECT_NEAR(sbox::hydrogen::radial_wavefunction(1, 0, 1.0, 0.0), 2.0, 1e-6);
    EXPECT_NEAR(sbox::hydrogen::radial_wavefunction(1, 0, 1.0, 1.0), 2.0 * std::exp(-1.0), 1e-6);
    EXPECT_NEAR(sbox::hydrogen::radial_wavefunction(2, 0, 1.0, 2.0), 0.0, 1e-4);
    EXPECT_NEAR(sbox::hydrogen::radial_wavefunction(2, 1, 1.0, 0.0), 0.0, 1e-6);
}

TEST(HydrogenTest, RadialWavefunctionNormalization) {
    const auto r10 = [](double r) {
        const double v = sbox::hydrogen::radial_wavefunction(1, 0, 1.0, r);
        return r * r * v * v;
    };
    const auto r21 = [](double r) {
        const double v = sbox::hydrogen::radial_wavefunction(2, 1, 1.0, r);
        return r * r * v * v;
    };
    const auto r32 = [](double r) {
        const double v = sbox::hydrogen::radial_wavefunction(3, 2, 1.0, r);
        return r * r * v * v;
    };

    EXPECT_NEAR(simpson_integral(r10, 0.0, 50.0, 10000), 1.0, 1e-3);
    EXPECT_NEAR(simpson_integral(r21, 0.0, 80.0, 10000), 1.0, 1e-3);
    EXPECT_NEAR(simpson_integral(r32, 0.0, 150.0, 10000), 1.0, 1e-3);
}

TEST(HydrogenTest, OrbitalValues) {
    const double psi100 = sbox::hydrogen::orbital_value(1, 0, 0, 1.0, 0.0, 0.0, 1.0);
    const double expected_psi100 = (2.0 * std::exp(-1.0)) * (1.0 / (2.0 * std::sqrt(kPi)));
    EXPECT_NEAR(psi100, expected_psi100, 1e-6);

    const double psi210_z = sbox::hydrogen::orbital_value(2, 1, 0, 1.0, 0.0, 0.0, 1.0);
    EXPECT_NE(psi210_z, 0.0);

    const double psi210_x = sbox::hydrogen::orbital_value(2, 1, 0, 1.0, 1.0, 0.0, 0.0);
    EXPECT_NEAR(psi210_x, 0.0, 1e-6);
}
