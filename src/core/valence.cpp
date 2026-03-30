#include "core/valence.h"

#include "core/covalent_radii.h"

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace sbox::chem {

namespace {

double xh_bond_length(int Z) {
    switch (Z) {
    case 1: return 1.4;
    case 6: return 2.06;
    case 7: return 1.91;
    case 8: return 1.81;
    case 16: return 2.53;
    default: return 1.9;
    }
}

Eigen::Vector3d normalized_or(const Eigen::Vector3d& v, const Eigen::Vector3d& fallback) {
    if (v.norm() <= 1.0e-10) {
        return fallback.normalized();
    }
    return v.normalized();
}

Eigen::Vector3d any_perpendicular(const Eigen::Vector3d& axis) {
    Eigen::Vector3d ref = std::abs(axis.z()) < 0.9 ? Eigen::Vector3d::UnitZ() : Eigen::Vector3d::UnitX();
    Eigen::Vector3d perp = axis.cross(ref);
    if (perp.norm() <= 1.0e-10) {
        ref = Eigen::Vector3d::UnitY();
        perp = axis.cross(ref);
    }
    return perp.normalized();
}

}  // namespace

int default_valence(int Z) {
    switch (Z) {
    case 1: return 1;
    case 2: return 0;
    case 3: return 1;
    case 4: return 2;
    case 5: return 3;
    case 6: return 4;
    case 7: return 3;
    case 8: return 2;
    case 9: return 1;
    case 10: return 0;
    case 11: return 1;
    case 12: return 2;
    case 13: return 3;
    case 14: return 4;
    case 15: return 3;
    case 16: return 2;
    case 17: return 1;
    case 19: return 1;
    case 20: return 2;
    case 29: return 1;
    case 30: return 2;
    case 35: return 1;
    case 53: return 1;
    default: return -1;
    }
}

int current_valence(const MolecularSystem& mol, int atom_index) {
    int valence = std::abs(mol.atom(atom_index).formal_charge);
    for (const Bond& bond : mol.bonds()) {
        if (bond.atom_i != atom_index && bond.atom_j != atom_index) {
            continue;
        }
        switch (bond.order) {
        case BondOrder::Single: valence += 1; break;
        case BondOrder::Double: valence += 2; break;
        case BondOrder::Triple: valence += 3; break;
        case BondOrder::Aromatic: valence += 2; break;
        case BondOrder::Unknown: valence += 1; break;
        }
    }
    return valence;
}

int missing_hydrogens(const MolecularSystem& mol, int atom_index) {
    const int target = default_valence(mol.atom(atom_index).Z);
    if (target < 0) {
        return 0;
    }
    return std::max(0, target - current_valence(mol, atom_index));
}

std::vector<Eigen::Vector3d> distribute_around_axis(
    const Eigen::Vector3d& axis,
    double angle_from_axis,
    int count,
    const Eigen::Vector3d& reference_dir) {
    std::vector<Eigen::Vector3d> result;
    if (count <= 0) {
        return result;
    }

    const Eigen::Vector3d a = normalized_or(axis, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d ref = reference_dir;
    if (ref.norm() <= 1.0e-10 || std::abs(ref.normalized().dot(a)) > 0.999) {
        ref = any_perpendicular(a);
    } else {
        ref = (ref - ref.dot(a) * a).normalized();
        if (ref.norm() <= 1.0e-10) {
            ref = any_perpendicular(a);
        }
    }

    const Eigen::Vector3d base = Eigen::AngleAxisd(angle_from_axis, ref.cross(a).normalized()) * a;
    for (int i = 0; i < count; ++i) {
        const double phi = 2.0 * std::acos(-1.0) * static_cast<double>(i) / static_cast<double>(count);
        result.push_back((Eigen::AngleAxisd(phi, a) * base).normalized());
    }
    return result;
}

std::vector<Eigen::Vector3d> compute_hydrogen_positions(
    const MolecularSystem& mol,
    int atom_index,
    int num_hydrogens) {
    std::vector<Eigen::Vector3d> positions;
    if (num_hydrogens <= 0) {
        return positions;
    }

    const Atom& center_atom = mol.atom(atom_index);
    const Eigen::Vector3d center = center_atom.position;
    const double bond_length = xh_bond_length(center_atom.Z);
    const std::vector<int> neighbors = mol.neighbors(atom_index);

    std::vector<Eigen::Vector3d> neighbor_dirs;
    for (int neighbor : neighbors) {
        neighbor_dirs.push_back((mol.atom(neighbor).position - center).normalized());
    }

    const int total_connections = static_cast<int>(neighbors.size()) + num_hydrogens;
    std::vector<Eigen::Vector3d> dirs;

    if (neighbors.empty()) {
        if (num_hydrogens == 1) {
            dirs.push_back(Eigen::Vector3d::UnitZ());
        } else if (num_hydrogens == 2) {
            dirs = distribute_around_axis(Eigen::Vector3d::UnitZ(), std::acos(1.0 / std::sqrt(3.0)), 2, Eigen::Vector3d::UnitX());
        } else {
            const std::vector<Eigen::Vector3d> tetra = {
                Eigen::Vector3d(1.0, 1.0, 1.0).normalized(),
                Eigen::Vector3d(-1.0, -1.0, 1.0).normalized(),
                Eigen::Vector3d(-1.0, 1.0, -1.0).normalized(),
                Eigen::Vector3d(1.0, -1.0, -1.0).normalized(),
            };
            dirs.assign(tetra.begin(), tetra.begin() + std::min<int>(num_hydrogens, static_cast<int>(tetra.size())));
        }
    } else if (neighbors.size() == 1) {
        const Eigen::Vector3d away = -neighbor_dirs[0];
        if (total_connections == 2) {
            dirs.push_back(away);
        } else if (total_connections == 3) {
            dirs = distribute_around_axis(away, std::acos(0.5), 2, any_perpendicular(away));
        } else {
            dirs = distribute_around_axis(away, std::acos(1.0 / 3.0), 3, any_perpendicular(away));
        }
    } else if (neighbors.size() == 2) {
        const Eigen::Vector3d away = normalized_or(-(neighbor_dirs[0] + neighbor_dirs[1]), any_perpendicular(neighbor_dirs[0]));
        const Eigen::Vector3d plane_normal = normalized_or(neighbor_dirs[0].cross(neighbor_dirs[1]), any_perpendicular(away));
        if (total_connections == 3) {
            dirs.push_back(away);
        } else {
            dirs = distribute_around_axis(away, std::acos(1.0 / std::sqrt(3.0)), 2, plane_normal);
        }
    } else {
        Eigen::Vector3d away = Eigen::Vector3d::Zero();
        for (const auto& dir : neighbor_dirs) {
            away -= dir;
        }
        away = normalized_or(away, neighbor_dirs[0].cross(neighbor_dirs[1]));
        dirs.push_back(away);
    }

    for (const Eigen::Vector3d& dir : dirs) {
        positions.push_back(center + bond_length * dir.normalized());
    }
    return positions;
}

void add_hydrogens(MolecularSystem& mol, int atom_index) {
    if (mol.atom(atom_index).Z == 1) {
        return;
    }
    const int count = missing_hydrogens(mol, atom_index);
    const std::vector<Eigen::Vector3d> positions = compute_hydrogen_positions(mol, atom_index, count);
    for (const Eigen::Vector3d& pos : positions) {
        const int h_index = mol.add_atom({1, pos, "", 0});
        mol.add_bond(atom_index, h_index, BondOrder::Single);
    }
}

void add_hydrogens(MolecularSystem& mol) {
    const int original_atom_count = mol.num_atoms();
    for (int i = 0; i < original_atom_count; ++i) {
        if (mol.atom(i).Z != 1) {
            add_hydrogens(mol, i);
        }
    }
}

void remove_hydrogens(MolecularSystem& mol) {
    for (int i = mol.num_atoms() - 1; i >= 0; --i) {
        if (mol.atom(i).Z == 1) {
            mol.remove_atom(i);
        }
    }
}

}  // namespace sbox::chem
