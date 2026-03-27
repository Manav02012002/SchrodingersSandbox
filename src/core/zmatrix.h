#pragma once

#include "core/molecular_system.h"

#include <vector>

namespace sbox::chem {

struct ZMatrixEntry {
    int Z;
    int bond_ref = -1;
    double bond_length = 0.0;
    int angle_ref = -1;
    double bond_angle = 0.0;
    int dihedral_ref = -1;
    double dihedral_angle = 0.0;
};

MolecularSystem zmatrix_to_cartesian(const std::vector<ZMatrixEntry>& zmat);

}  // namespace sbox::chem
