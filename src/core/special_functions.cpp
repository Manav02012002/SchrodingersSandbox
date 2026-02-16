#include "core/special_functions.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace sbox::math {

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

const std::array<double, 171>& factorial_table() {
    static const std::array<double, 171> table = [] {
        std::array<double, 171> values{};
        values[0] = 1.0;
        for (int i = 1; i <= 170; ++i) {
            values[static_cast<std::size_t>(i)] =
                values[static_cast<std::size_t>(i - 1)] * static_cast<double>(i);
        }
        return values;
    }();
    return table;
}

}  // namespace

double factorial(int n) {
    if (n < 0) {
        return 1.0;
    }
    if (n <= 170) {
        return factorial_table()[static_cast<std::size_t>(n)];
    }
    return std::numeric_limits<double>::infinity();
}

double double_factorial(int n) {
    if (n <= 0) {
        return 1.0;
    }

    double result = 1.0;
    for (int k = n; k > 1; k -= 2) {
        result *= static_cast<double>(k);
    }
    return result;
}

double associated_laguerre(int n, double alpha, double x) {
    if (n < 0) {
        return 0.0;
    }
    if (n == 0) {
        return 1.0;
    }

    double l_nm2 = 1.0;
    double l_nm1 = 1.0 + alpha - x;
    if (n == 1) {
        return l_nm1;
    }

    for (int k = 2; k <= n; ++k) {
        const double kd = static_cast<double>(k);
        const double term1 = (2.0 * kd - 1.0 + alpha - x) * l_nm1;
        const double term2 = (kd - 1.0 + alpha) * l_nm2;
        const double l_n = (term1 - term2) / kd;
        l_nm2 = l_nm1;
        l_nm1 = l_n;
    }
    return l_nm1;
}

double associated_legendre(int l, int m, double x) {
    if (l < 0) {
        return 0.0;
    }

    const int m_abs = std::abs(m);
    if (m_abs > l) {
        return 0.0;
    }

    x = std::clamp(x, -1.0, 1.0);

    double p_mm = 1.0;
    if (m_abs > 0) {
        const double one_minus_x2 = std::max(0.0, 1.0 - x * x);
        p_mm = ((m_abs % 2) == 0 ? 1.0 : -1.0) *
               double_factorial(2 * m_abs - 1) *
               std::pow(one_minus_x2, 0.5 * static_cast<double>(m_abs));
    }

    if (l == m_abs) {
        return p_mm;
    }

    const double p_m1m = x * (2.0 * static_cast<double>(m_abs) + 1.0) * p_mm;
    if (l == m_abs + 1) {
        return p_m1m;
    }

    double p_lm2 = p_mm;
    double p_lm1 = p_m1m;
    for (int ell = m_abs + 2; ell <= l; ++ell) {
        const double ell_d = static_cast<double>(ell);
        const double numerator =
            (2.0 * ell_d - 1.0) * x * p_lm1 -
            static_cast<double>(ell + m_abs - 1) * p_lm2;
        const double p_lm = numerator / static_cast<double>(ell - m_abs);
        p_lm2 = p_lm1;
        p_lm1 = p_lm;
    }

    return p_lm1;
}

double real_spherical_harmonic(int l, int m, double theta, double phi) {
    if (l < 0 || std::abs(m) > l) {
        return 0.0;
    }

    const int m_abs = std::abs(m);
    const double x = std::cos(theta);
    const double p_lm = associated_legendre(l, m_abs, x);
    const double normalization = std::sqrt(
        ((2.0 * static_cast<double>(l) + 1.0) / (4.0 * kPi)) *
        (factorial(l - m_abs) / factorial(l + m_abs)));

    if (m > 0) {
        return std::sqrt(2.0) * normalization * p_lm * std::cos(static_cast<double>(m) * phi);
    }
    if (m < 0) {
        return std::sqrt(2.0) * normalization * p_lm * std::sin(static_cast<double>(m_abs) * phi);
    }
    return normalization * p_lm;
}

}  // namespace sbox::math
