#pragma once

#include "chem/coordination.h"
#include "core/molecular_system.h"

#include <string>
#include <vector>

namespace sbox::chem {

enum class LigandDenticity {
    Monodentate = 1,
    Bidentate = 2,
    Tridentate = 3,
    Tetradentate = 4,
    Hexadentate = 6,
};

struct LigandEntry {
    std::string name;
    std::string abbreviation;
    std::string formula;
    std::string category;
    LigandDenticity denticity = LigandDenticity::Monodentate;
    int charge = 0;
    MolecularSystem molecule;
    std::vector<int> donor_atoms;
    std::string field_strength;
    double typical_delta_oct_eV = 0.0;
};

class LigandLibrary {
public:
    LigandLibrary();

    const std::vector<LigandEntry>& all() const;
    std::vector<const LigandEntry*> by_category(const std::string& category) const;
    std::vector<const LigandEntry*> by_denticity(LigandDenticity d) const;
    const LigandEntry* find(const std::string& name_or_abbrev) const;
    std::vector<std::string> categories() const;

    LigandSpec to_ligand_spec(const LigandEntry& entry, int donor_index = 0) const;

private:
    std::vector<LigandEntry> ligands_;
    void build_library();
};

}  // namespace sbox::chem
