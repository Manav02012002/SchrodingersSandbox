#include "core/special_functions.h"

#include <cmath>

#include <gtest/gtest.h>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

}  // namespace

TEST(SpecialFunctionsTest, AssociatedLaguerreValues) {
    EXPECT_NEAR(sbox::math::associated_laguerre(0, 0.0, 5.0), 1.0, 1e-10);
    EXPECT_NEAR(sbox::math::associated_laguerre(1, 1.0, 3.0), -1.0, 1e-10);
    EXPECT_NEAR(sbox::math::associated_laguerre(2, 0.0, 1.0), -0.5, 1e-10);
    EXPECT_NEAR(sbox::math::associated_laguerre(3, 1.0, 2.0), -4.0 / 3.0, 1e-10);
}

TEST(SpecialFunctionsTest, AssociatedLegendreValues) {
    const double x = 0.5;
    EXPECT_NEAR(sbox::math::associated_legendre(0, 0, x), 1.0, 1e-10);
    EXPECT_NEAR(sbox::math::associated_legendre(1, 0, x), 0.5, 1e-10);
    EXPECT_NEAR(sbox::math::associated_legendre(1, 1, x), -std::sqrt(0.75), 1e-10);
    EXPECT_NEAR(sbox::math::associated_legendre(2, 0, x), -0.125, 1e-10);
    EXPECT_NEAR(sbox::math::associated_legendre(2, 1, x), -3.0 * x * std::sqrt(0.75), 1e-10);
    EXPECT_NEAR(sbox::math::associated_legendre(2, 2, x), 2.25, 1e-10);
}

TEST(SpecialFunctionsTest, RealSphericalHarmonicValues) {
    EXPECT_NEAR(sbox::math::real_spherical_harmonic(0, 0, 0.3, 1.2), 1.0 / (2.0 * std::sqrt(kPi)), 1e-8);
    EXPECT_NEAR(sbox::math::real_spherical_harmonic(1, 0, 0.0, 0.0), std::sqrt(3.0 / (4.0 * kPi)), 1e-8);

    const double y11_equator = sbox::math::real_spherical_harmonic(1, 1, kPi / 2.0, 0.0);
    EXPECT_NE(y11_equator, 0.0);
    EXPECT_NEAR(sbox::math::real_spherical_harmonic(1, 1, 0.0, 0.0), 0.0, 1e-8);

    const double expected_y20 = -0.5 * std::sqrt(5.0 / (4.0 * kPi));
    EXPECT_NEAR(sbox::math::real_spherical_harmonic(2, 0, kPi / 2.0, 0.0), expected_y20, 1e-8);
}
