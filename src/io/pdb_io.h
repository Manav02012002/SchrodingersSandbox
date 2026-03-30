#pragma once

#include "core/molecular_system.h"

#include <Eigen/Core>

#include <string>
#include <utility>
#include <vector>

namespace sbox::io {

struct PDBAtom {
    int serial = 0;
    std::string name;
    std::string alt_loc;
    std::string residue_name;
    std::string chain_id;
    int residue_seq = 0;
    Eigen::Vector3d position = Eigen::Vector3d::Zero();  // Angstrom
    double occupancy = 0.0;
    double b_factor = 0.0;
    std::string element;
    int Z = 0;
};

struct PDBResidue {
    std::string name;
    int sequence = 0;
    std::string chain;
    std::vector<int> atom_indices;
};

struct PDBChain {
    std::string id;
    std::vector<int> residue_indices;
};

struct PDBData {
    std::string title;
    std::vector<PDBAtom> atoms;
    std::vector<PDBResidue> residues;
    std::vector<PDBChain> chains;
    std::vector<std::pair<int, int>> conect_bonds;

    sbox::chem::MolecularSystem to_molecular_system() const;
};

PDBData read_pdb(const std::string& filepath);
void write_pdb(const std::string& filepath, const PDBData& data);

}  // namespace sbox::io
