#pragma once

#include <Eigen/Core>

#include <vector>

namespace sbox::basis {

struct GaussianPrimitive {
    double exponent;
    double coefficient;  // contraction coefficient (already includes normalisation)
};

struct BasisShell {
    int atom_index;  // which atom this shell sits on
    int angular_momentum;  // 0=s, 1=p, 2=d, 3=f
    std::vector<GaussianPrimitive> primitives;
};

struct BasisSet {
    std::vector<BasisShell> shells;
    int num_basis_functions() const;  // total count: 1 per s, 3 per p, 6 per d (cartesian) or 5 (spherical)
    bool spherical = true;  // spherical (5d,7f) vs cartesian (6d,10f)
};

struct MOData {
    BasisSet basis;
    std::vector<int> atomic_numbers;
    std::vector<Eigen::Vector3d> atom_positions;
    Eigen::MatrixXd coefficients;  // num_basis × num_mo
    Eigen::VectorXd energies;  // num_mo
    Eigen::VectorXd occupations;  // num_mo
    double total_energy;
};

}  // namespace sbox::basis
