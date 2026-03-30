#pragma once

#include "core/basis_set.h"
#include "core/molecular_system.h"

#include <string>
#include <utility>
#include <vector>

namespace sbox::analysis {

struct AtomContribution {
    int atom_index = -1;
    int Z = 0;
    std::string element;
    double total_weight = 0.0;
    std::vector<std::pair<std::string, double>> ao_contributions;
};

struct OrbitalComposition {
    int mo_index = -1;
    double energy_eV = 0.0;
    double occupation = 0.0;
    std::vector<AtomContribution> atom_contributions;
    std::string summary;
};

OrbitalComposition analyze_orbital_composition(
    const sbox::basis::MOData& mo_data,
    const sbox::chem::MolecularSystem& mol,
    int mo_index);

}  // namespace sbox::analysis
