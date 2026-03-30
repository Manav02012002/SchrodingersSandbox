#include "core/spline.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

using sbox::math::CubicSpline;
constexpr double kPi = 3.14159265358979323846;

TEST(Spline, FitsSinAccurately) {
    CubicSpline spline;
    std::vector<double> x(10);
    std::vector<double> y(10);
    const double two_pi = 2.0 * kPi;
    for (int i = 0; i < 10; ++i) {
        x[static_cast<std::size_t>(i)] = two_pi * static_cast<double>(i) / 9.0;
        y[static_cast<std::size_t>(i)] = std::sin(x[static_cast<std::size_t>(i)]);
    }
    spline.fit(x, y);

    double max_error = 0.0;
    for (int i = 0; i < 100; ++i) {
        const double xi = two_pi * static_cast<double>(i) / 99.0;
        max_error = std::max(max_error, std::abs(spline.evaluate(xi) - std::sin(xi)));
    }
    EXPECT_LT(max_error, 0.01);
}

TEST(Spline, ReproducesLinearFunction) {
    CubicSpline spline;
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> y;
    for (double xi : x) {
        y.push_back(2.0 * xi + 1.0);
    }
    spline.fit(x, y);

    for (int i = 0; i <= 30; ++i) {
        const double xi = 3.0 * static_cast<double>(i) / 30.0;
        EXPECT_NEAR(spline.evaluate(xi), 2.0 * xi + 1.0, 1.0e-12);
    }
}

TEST(Spline, ParabolaApproximation) {
    CubicSpline spline;
    spline.fit({0.0, 1.0, 2.0}, {0.0, 1.0, 4.0});
    EXPECT_NEAR(spline.evaluate(0.5), 0.25, 0.1);
}

TEST(Spline, DerivativeMatchesCosine) {
    CubicSpline spline;
    std::vector<double> x(10);
    std::vector<double> y(10);
    const double two_pi = 2.0 * kPi;
    for (int i = 0; i < 10; ++i) {
        x[static_cast<std::size_t>(i)] = two_pi * static_cast<double>(i) / 9.0;
        y[static_cast<std::size_t>(i)] = std::sin(x[static_cast<std::size_t>(i)]);
    }
    spline.fit(x, y);
    EXPECT_NEAR(spline.derivative(kPi / 4.0), std::cos(kPi / 4.0), 0.05);
}

TEST(Spline, EvaluateGridProducesSmoothOutput) {
    CubicSpline spline;
    spline.fit({0.0, 1.0, 2.0, 3.0}, {1.0, 0.0, 1.0, 0.0});
    const auto [xs, ys] = spline.evaluate_grid(0.0, 3.0, 25);
    EXPECT_EQ(xs.size(), 25);
    EXPECT_EQ(ys.size(), 25);
    for (std::size_t i = 1; i < xs.size(); ++i) {
        EXPECT_GT(xs[i], xs[i - 1]);
        EXPECT_LT(std::abs(ys[i] - ys[i - 1]), 0.5);
    }
}

TEST(Spline, TwoPointCaseIsLinear) {
    CubicSpline spline;
    spline.fit({1.0, 3.0}, {2.0, 6.0});
    EXPECT_NEAR(spline.evaluate(2.0), 4.0, 1.0e-12);
    EXPECT_NEAR(spline.derivative(2.0), 2.0, 1.0e-12);
    EXPECT_TRUE(spline.is_fitted());
    EXPECT_EQ(spline.size(), 2);
}

}  // namespace
