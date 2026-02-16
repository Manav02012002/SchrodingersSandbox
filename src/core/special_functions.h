#pragma once

namespace sbox::math {

double factorial(int n);
double double_factorial(int n);
double associated_laguerre(int n, double alpha, double x);
double associated_legendre(int l, int m, double x);
double real_spherical_harmonic(int l, int m, double theta, double phi);

}  // namespace sbox::math
