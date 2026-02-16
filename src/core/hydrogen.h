#pragma once

namespace sbox::hydrogen {

inline constexpr double a0 = 1.0;
inline constexpr double a0_angstrom = 0.529177;

double radial_wavefunction(int n, int l, double zeff, double r);
double orbital_value(int n, int l, int m, double zeff, double x, double y, double z);
double probability_density(int n, int l, int m, double zeff, double x, double y, double z);

}  // namespace sbox::hydrogen
