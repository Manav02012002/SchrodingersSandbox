#include "core/zmatrix.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace sbox::chem {
namespace {

void validate_ref(int ref, int max_index, const char* name) {
    if (ref < 0 || ref >= max_index) {
        throw std::invalid_argument(std::string("Invalid Z-matrix reference for ") + name);
    }
}

Eigen::Vector3d choose_perpendicular(const Eigen::Vector3d& axis) {
    const Eigen::Vector3d helper =
        std::abs(axis.dot(Eigen::Vector3d::UnitX())) < 0.9 ? Eigen::Vector3d::UnitX() : Eigen::Vector3d::UnitY();
    Eigen::Vector3d perp = helper - helper.dot(axis) * axis;
    const double norm = perp.norm();
    if (norm <= 1e-12) {
        throw std::invalid_argument("Could not construct perpendicular vector for Z-matrix");
    }
    return perp / norm;
}

}  // namespace

MolecularSystem zmatrix_to_cartesian(const std::vector<ZMatrixEntry>& zmat) {
    MolecularSystem mol;
    if (zmat.empty()) {
        return mol;
    }

    std::vector<Eigen::Vector3d> positions;
    positions.reserve(zmat.size());

    for (std::size_t i = 0; i < zmat.size(); ++i) {
        const ZMatrixEntry& entry = zmat[i];
        Eigen::Vector3d pos = Eigen::Vector3d::Zero();

        if (i == 0) {
            pos = Eigen::Vector3d::Zero();
        } else if (i == 1) {
            validate_ref(entry.bond_ref, static_cast<int>(i), "bond_ref");
            pos = positions[static_cast<std::size_t>(entry.bond_ref)] +
                  Eigen::Vector3d(0.0, 0.0, entry.bond_length);
        } else if (i == 2) {
            validate_ref(entry.bond_ref, static_cast<int>(i), "bond_ref");
            validate_ref(entry.angle_ref, static_cast<int>(i), "angle_ref");

            const Eigen::Vector3d& a = positions[static_cast<std::size_t>(entry.bond_ref)];
            const Eigen::Vector3d& b = positions[static_cast<std::size_t>(entry.angle_ref)];
            const Eigen::Vector3d axis = (b - a).normalized();
            const Eigen::Vector3d perp = choose_perpendicular(axis);
            pos = a + entry.bond_length *
                          (axis * std::cos(entry.bond_angle) + perp * std::sin(entry.bond_angle));
        } else {
            validate_ref(entry.bond_ref, static_cast<int>(i), "bond_ref");
            validate_ref(entry.angle_ref, static_cast<int>(i), "angle_ref");
            validate_ref(entry.dihedral_ref, static_cast<int>(i), "dihedral_ref");

            const Eigen::Vector3d& a = positions[static_cast<std::size_t>(entry.bond_ref)];
            const Eigen::Vector3d& b = positions[static_cast<std::size_t>(entry.angle_ref)];
            const Eigen::Vector3d& c = positions[static_cast<std::size_t>(entry.dihedral_ref)];

            const Eigen::Vector3d v1 = (b - a).normalized();
            const Eigen::Vector3d v2 = a - b;
            const Eigen::Vector3d v3 = b - c;
            const Eigen::Vector3d n = v2.cross(v3).normalized();
            const Eigen::Vector3d nn = n.cross(v1).normalized();

            pos = a + entry.bond_length *
                          (v1 * std::cos(entry.bond_angle) +
                           nn * std::sin(entry.bond_angle) * std::cos(entry.dihedral_angle) +
                           n * std::sin(entry.bond_angle) * std::sin(entry.dihedral_angle));
        }

        positions.push_back(pos);
        mol.add_atom({entry.Z, pos, "", 0});
    }

    mol.perceive_bonds();
    return mol;
}

}  // namespace sbox::chem
