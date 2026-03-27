#include "core/molecular_system.h"

#include "core/covalent_radii.h"
#include "core/elements.h"

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace sbox::chem {
namespace {

void validate_atom_index(int i, int size) {
    if (i < 0 || i >= size) {
        throw std::out_of_range("Atom index out of range");
    }
}

void validate_bond_index(int i, int size) {
    if (i < 0 || i >= size) {
        throw std::out_of_range("Bond index out of range");
    }
}

}  // namespace

int MolecularSystem::num_atoms() const {
    return static_cast<int>(atoms_.size());
}

int MolecularSystem::num_bonds() const {
    return static_cast<int>(bonds_.size());
}

const Atom& MolecularSystem::atom(int i) const {
    validate_atom_index(i, num_atoms());
    return atoms_[static_cast<std::size_t>(i)];
}

Atom& MolecularSystem::atom(int i) {
    validate_atom_index(i, num_atoms());
    return atoms_[static_cast<std::size_t>(i)];
}

const Bond& MolecularSystem::bond(int i) const {
    validate_bond_index(i, num_bonds());
    return bonds_[static_cast<std::size_t>(i)];
}

int MolecularSystem::add_atom(Atom a) {
    atoms_.push_back(std::move(a));
    return num_atoms() - 1;
}

void MolecularSystem::remove_atom(int i) {
    validate_atom_index(i, num_atoms());
    atoms_.erase(atoms_.begin() + i);

    std::vector<Bond> new_bonds;
    new_bonds.reserve(bonds_.size());
    for (const Bond& bond : bonds_) {
        if (bond.atom_i == i || bond.atom_j == i) {
            continue;
        }

        Bond updated = bond;
        if (updated.atom_i > i) {
            --updated.atom_i;
        }
        if (updated.atom_j > i) {
            --updated.atom_j;
        }
        new_bonds.push_back(updated);
    }

    bonds_ = std::move(new_bonds);
}

int MolecularSystem::add_bond(int i, int j, BondOrder order) {
    validate_atom_index(i, num_atoms());
    validate_atom_index(j, num_atoms());
    if (i == j) {
        throw std::invalid_argument("Bond endpoints must be distinct atoms");
    }

    bonds_.push_back(Bond{i, j, order});
    return num_bonds() - 1;
}

void MolecularSystem::remove_bond(int i) {
    validate_bond_index(i, num_bonds());
    bonds_.erase(bonds_.begin() + i);
}

bool MolecularSystem::has_bond(int i, int j) const {
    validate_atom_index(i, num_atoms());
    validate_atom_index(j, num_atoms());

    for (const Bond& bond : bonds_) {
        if ((bond.atom_i == i && bond.atom_j == j) ||
            (bond.atom_i == j && bond.atom_j == i)) {
            return true;
        }
    }
    return false;
}

void MolecularSystem::perceive_bonds(double tolerance) {
    bonds_.clear();

    for (int i = 0; i < num_atoms(); ++i) {
        for (int j = i + 1; j < num_atoms(); ++j) {
            const double cutoff = tolerance * (covalent_radius(atom(i).Z) + covalent_radius(atom(j).Z));
            if (distance(i, j) < cutoff) {
                bonds_.push_back(Bond{i, j, BondOrder::Single});
            }
        }
    }
}

std::vector<int> MolecularSystem::neighbors(int atom_i) const {
    validate_atom_index(atom_i, num_atoms());

    std::vector<int> result;
    for (const Bond& bond : bonds_) {
        if (bond.atom_i == atom_i) {
            result.push_back(bond.atom_j);
        } else if (bond.atom_j == atom_i) {
            result.push_back(bond.atom_i);
        }
    }
    return result;
}

int MolecularSystem::coordination_number(int atom_i) const {
    return static_cast<int>(neighbors(atom_i).size());
}

double MolecularSystem::distance(int i, int j) const {
    return (atom(i).position - atom(j).position).norm();
}

double MolecularSystem::angle(int i, int j, int k) const {
    const Eigen::Vector3d v1 = atom(i).position - atom(j).position;
    const Eigen::Vector3d v2 = atom(k).position - atom(j).position;
    const double norm1 = v1.norm();
    const double norm2 = v2.norm();
    if (norm1 <= 0.0 || norm2 <= 0.0) {
        throw std::invalid_argument("Angle is undefined for zero-length vectors");
    }

    const double cos_theta = std::clamp(v1.dot(v2) / (norm1 * norm2), -1.0, 1.0);
    return std::acos(cos_theta);
}

double MolecularSystem::dihedral(int i, int j, int k, int l) const {
    const Eigen::Vector3d p0 = atom(i).position;
    const Eigen::Vector3d p1 = atom(j).position;
    const Eigen::Vector3d p2 = atom(k).position;
    const Eigen::Vector3d p3 = atom(l).position;

    const Eigen::Vector3d b0 = -(p1 - p0);
    const Eigen::Vector3d b1 = p2 - p1;
    const Eigen::Vector3d b2 = p3 - p2;

    const double b1_norm = b1.norm();
    if (b1_norm <= 0.0) {
        throw std::invalid_argument("Dihedral is undefined for zero-length central bond");
    }

    const Eigen::Vector3d b1_unit = b1 / b1_norm;
    const Eigen::Vector3d v = b0 - b0.dot(b1_unit) * b1_unit;
    const Eigen::Vector3d w = b2 - b2.dot(b1_unit) * b1_unit;

    const double v_norm = v.norm();
    const double w_norm = w.norm();
    if (v_norm <= 0.0 || w_norm <= 0.0) {
        throw std::invalid_argument("Dihedral is undefined for colinear atoms");
    }

    const Eigen::Vector3d v_unit = v / v_norm;
    const Eigen::Vector3d w_unit = w / w_norm;
    const double x = v_unit.dot(w_unit);
    const double y = b1_unit.cross(v_unit).dot(w_unit);
    return std::atan2(y, x);
}

Eigen::Vector3d MolecularSystem::center_of_mass() const {
    if (atoms_.empty()) {
        return Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d weighted_sum = Eigen::Vector3d::Zero();
    double total_mass = 0.0;
    for (const Atom& atom : atoms_) {
        const double mass = sbox::elements::get_element(atom.Z).atomic_mass;
        weighted_sum += mass * atom.position;
        total_mass += mass;
    }

    if (total_mass <= 0.0) {
        return Eigen::Vector3d::Zero();
    }
    return weighted_sum / total_mass;
}

void MolecularSystem::center() {
    const Eigen::Vector3d com = center_of_mass();
    for (Atom& atom : atoms_) {
        atom.position -= com;
    }
}

void MolecularSystem::clear() {
    atoms_.clear();
    bonds_.clear();
}

const std::string& MolecularSystem::name() const {
    return name_;
}

void MolecularSystem::set_name(const std::string& name) {
    name_ = name;
}

int MolecularSystem::charge() const {
    return charge_;
}

void MolecularSystem::set_charge(int charge) {
    charge_ = charge;
}

int MolecularSystem::multiplicity() const {
    return multiplicity_;
}

void MolecularSystem::set_multiplicity(int multiplicity) {
    multiplicity_ = multiplicity;
}

const std::vector<Atom>& MolecularSystem::atoms() const {
    return atoms_;
}

const std::vector<Bond>& MolecularSystem::bonds() const {
    return bonds_;
}

}  // namespace sbox::chem
