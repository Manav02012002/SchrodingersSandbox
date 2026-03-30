#include "core/spline.h"

#include <algorithm>
#include <stdexcept>

namespace sbox::math {

void CubicSpline::fit(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("CubicSpline::fit requires x and y to have the same size");
    }
    if (x.size() < 2) {
        throw std::invalid_argument("CubicSpline::fit requires at least two data points");
    }

    for (std::size_t i = 1; i < x.size(); ++i) {
        if (x[i] <= x[i - 1]) {
            throw std::invalid_argument("CubicSpline::fit requires x to be strictly increasing");
        }
    }

    x_ = x;
    y_ = y;

    const std::size_t n = x_.size() - 1;
    a_.assign(n, 0.0);
    b_.assign(n, 0.0);
    c_.assign(n, 0.0);
    d_.assign(n, 0.0);

    if (n == 1) {
        a_[0] = y_[0];
        b_[0] = (y_[1] - y_[0]) / (x_[1] - x_[0]);
        c_[0] = 0.0;
        d_[0] = 0.0;
        fitted_ = true;
        return;
    }

    std::vector<double> h(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        h[i] = x_[i + 1] - x_[i];
    }

    std::vector<double> lower(n - 1, 0.0);
    std::vector<double> diag(n - 1, 0.0);
    std::vector<double> upper(n - 1, 0.0);
    std::vector<double> rhs(n - 1, 0.0);

    for (std::size_t i = 1; i < n; ++i) {
        lower[i - 1] = h[i - 1];
        diag[i - 1] = 2.0 * (h[i - 1] + h[i]);
        upper[i - 1] = h[i];
        rhs[i - 1] =
            6.0 * (((y_[i + 1] - y_[i]) / h[i]) - ((y_[i] - y_[i - 1]) / h[i - 1]));
    }

    for (std::size_t i = 1; i < diag.size(); ++i) {
        const double factor = lower[i] / diag[i - 1];
        diag[i] -= factor * upper[i - 1];
        rhs[i] -= factor * rhs[i - 1];
    }

    std::vector<double> m(n + 1, 0.0);
    m[n - 1] = rhs.back() / diag.back();
    for (std::size_t i = diag.size() - 1; i-- > 0;) {
        m[i + 1] = (rhs[i] - upper[i] * m[i + 2]) / diag[i];
    }

    for (std::size_t i = 0; i < n; ++i) {
        a_[i] = y_[i];
        b_[i] = (y_[i + 1] - y_[i]) / h[i] - h[i] * (2.0 * m[i] + m[i + 1]) / 6.0;
        c_[i] = m[i] / 2.0;
        d_[i] = (m[i + 1] - m[i]) / (6.0 * h[i]);
    }

    fitted_ = true;
}

int CubicSpline::find_segment(double x) const {
    if (x_.size() < 2) {
        throw std::runtime_error("CubicSpline::find_segment called before fit");
    }
    if (x <= x_.front()) {
        return 0;
    }
    if (x >= x_.back()) {
        return static_cast<int>(x_.size()) - 2;
    }
    const auto it = std::upper_bound(x_.begin(), x_.end(), x);
    return static_cast<int>(std::distance(x_.begin(), it) - 1);
}

double CubicSpline::evaluate(double x) const {
    if (!fitted_) {
        throw std::runtime_error("CubicSpline::evaluate called before fit");
    }

    const int segment = find_segment(x);
    if (x <= x_.front()) {
        return y_.front() + b_.front() * (x - x_.front());
    }
    if (x >= x_.back()) {
        const int last = static_cast<int>(x_.size()) - 2;
        const double dx_end = x_.back() - x_[last];
        const double y_end = a_[last] + b_[last] * dx_end + c_[last] * dx_end * dx_end
                             + d_[last] * dx_end * dx_end * dx_end;
        const double slope_end = b_[last] + 2.0 * c_[last] * dx_end + 3.0 * d_[last] * dx_end * dx_end;
        return y_end + slope_end * (x - x_.back());
    }

    const double dx = x - x_[static_cast<std::size_t>(segment)];
    return a_[static_cast<std::size_t>(segment)] + b_[static_cast<std::size_t>(segment)] * dx
           + c_[static_cast<std::size_t>(segment)] * dx * dx
           + d_[static_cast<std::size_t>(segment)] * dx * dx * dx;
}

double CubicSpline::derivative(double x) const {
    if (!fitted_) {
        throw std::runtime_error("CubicSpline::derivative called before fit");
    }

    const int segment = find_segment(x);
    if (x <= x_.front()) {
        return b_.front();
    }
    if (x >= x_.back()) {
        const int last = static_cast<int>(x_.size()) - 2;
        const double dx_end = x_.back() - x_[last];
        return b_[last] + 2.0 * c_[last] * dx_end + 3.0 * d_[last] * dx_end * dx_end;
    }

    const double dx = x - x_[static_cast<std::size_t>(segment)];
    return b_[static_cast<std::size_t>(segment)] + 2.0 * c_[static_cast<std::size_t>(segment)] * dx
           + 3.0 * d_[static_cast<std::size_t>(segment)] * dx * dx;
}

std::pair<std::vector<double>, std::vector<double>> CubicSpline::evaluate_grid(
    double x_min, double x_max, int n_points) const {
    if (!fitted_) {
        throw std::runtime_error("CubicSpline::evaluate_grid called before fit");
    }
    if (n_points <= 0) {
        throw std::invalid_argument("CubicSpline::evaluate_grid requires n_points > 0");
    }

    std::vector<double> xs(static_cast<std::size_t>(n_points), 0.0);
    std::vector<double> ys(static_cast<std::size_t>(n_points), 0.0);
    if (n_points == 1) {
        xs[0] = x_min;
        ys[0] = evaluate(x_min);
        return {xs, ys};
    }

    const double step = (x_max - x_min) / static_cast<double>(n_points - 1);
    for (int i = 0; i < n_points; ++i) {
        xs[static_cast<std::size_t>(i)] = x_min + step * static_cast<double>(i);
        ys[static_cast<std::size_t>(i)] = evaluate(xs[static_cast<std::size_t>(i)]);
    }
    return {xs, ys};
}

bool CubicSpline::is_fitted() const {
    return fitted_;
}

int CubicSpline::size() const {
    return static_cast<int>(x_.size());
}

}  // namespace sbox::math
