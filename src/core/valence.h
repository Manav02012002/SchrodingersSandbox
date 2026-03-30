#pragma once

#include "core/molecular_system.h"

#include <Eigen/Core>

#include <vector>

namespace sbox::chem {

int default_valence(int Z);
int current_valence(const MolecularSystem& mol, int atom_index);
int missing_hydrogens(const MolecularSystem& mol, int atom_index);

std::vector<Eigen::Vector3d> distribute_around_axis(
    const Eigen::Vector3d& axis,
    double angle_from_axis,
    int count,
    const Eigen::Vector3d& reference_dir = Eigen::Vector3d::Zero());

std::vector<Eigen::Vector3d> compute_hydrogen_positions(
    const MolecularSystem& mol,
    int atom_index,
    int num_hydrogens);

void add_hydrogens(MolecularSystem& mol);
void add_hydrogens(MolecularSystem& mol, int atom_index);
void remove_hydrogens(MolecularSystem& mol);

}  // namespace sbox::chem
