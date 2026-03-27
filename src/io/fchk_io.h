#pragma once

#include <Eigen/Core>

#include <string>
#include <vector>

namespace sbox::io {

struct FchkData {
    std::string title;
    std::string method;
    std::string basis_name;

    int num_atoms = 0;
    int num_basis = 0;
    int num_mo = 0;
    int charge = 0;
    int multiplicity = 1;
    double total_energy = 0.0;

    std::vector<int> atomic_numbers;
    std::vector<double> coordinates;
    std::vector<int> shell_types;
    std::vector<int> shell_to_atom_map;
    std::vector<int> primitives_per_shell;
    std::vector<double> primitive_exponents;
    std::vector<double> contraction_coefficients;
    std::vector<double> sp_contraction_coefficients;
    Eigen::VectorXd mo_energies;
    Eigen::MatrixXd mo_coefficients;
    Eigen::VectorXd occupations;
    Eigen::Vector3d dipole_moment = Eigen::Vector3d::Zero();

    std::vector<double> mulliken_charges;
};

FchkData read_fchk(const std::string& filepath);

}  // namespace sbox::io
