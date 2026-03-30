#pragma once

#include "chem/coordination.h"
#include "core/basis_set.h"
#include "core/molecular_system.h"

#include <string>
#include <vector>

namespace sbox::analysis {

struct DOrbitalEnergies {
    double dxy = 0.0;
    double dxz = 0.0;
    double dyz = 0.0;
    double dz2 = 0.0;
    double dx2y2 = 0.0;
    int mo_dxy = -1;
    int mo_dxz = -1;
    int mo_dyz = -1;
    int mo_dz2 = -1;
    int mo_dx2y2 = -1;
    bool ordering_warning = false;

    struct OrbitalGroup {
        std::string label;
        std::vector<std::string> orbitals;
        double average_energy = 0.0;
        int degeneracy = 0;
    };
    std::vector<OrbitalGroup> groups;

    double delta_oct() const;
    double delta_tet() const;
    double mean_energy() const;
};

DOrbitalEnergies extract_d_orbitals(const sbox::basis::MOData& mo_data,
                                    const sbox::chem::MolecularSystem& mol,
                                    int metal_atom_index);

void identify_splitting(DOrbitalEnergies& d_orbs, sbox::chem::CoordinationGeometry geom);

int d_electron_count(int metal_Z, int oxidation_state);
double octahedral_cfse_dq(int d_electrons, bool high_spin);

}  // namespace sbox::analysis
