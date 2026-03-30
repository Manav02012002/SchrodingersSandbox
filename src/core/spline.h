#pragma once

#include <utility>
#include <vector>

namespace sbox::math {

class CubicSpline {
public:
    void fit(const std::vector<double>& x, const std::vector<double>& y);

    double evaluate(double x) const;
    double derivative(double x) const;

    std::pair<std::vector<double>, std::vector<double>> evaluate_grid(
        double x_min, double x_max, int n_points) const;

    bool is_fitted() const;
    int size() const;

private:
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> a_;
    std::vector<double> b_;
    std::vector<double> c_;
    std::vector<double> d_;
    bool fitted_ = false;

    int find_segment(double x) const;
};

}  // namespace sbox::math
