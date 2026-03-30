#pragma once

#include "core/molecular_system.h"

#include <Eigen/Core>

#include <string>
#include <vector>

namespace sbox::chem {

enum class CoordinationGeometry {
    Linear,
    TrigonalPlanar,
    TShaped,
    Tetrahedral,
    SquarePlanar,
    SeeSaw,
    TrigonalBipyramidal,
    SquarePyramidal,
    Octahedral,
    PentagonalBipyramidal,
    SquareAntiprismatic,
};

struct CoordinationTemplate {
    CoordinationGeometry geometry;
    std::string name;
    int coordination_number;
    std::vector<Eigen::Vector3d> directions;
};

struct LigandSpec {
    int donor_Z = 0;
    MolecularSystem ligand;
    double bond_length = 0.0;
    std::string name;
};

const CoordinationTemplate& get_template(CoordinationGeometry geom);
double default_ml_bond_length(int metal_Z, int ligand_Z);
Eigen::Matrix3d orient_ligand(const Eigen::Vector3d& current_direction,
                              const Eigen::Vector3d& target_direction);
MolecularSystem assemble_complex(int metal_Z,
                                 int metal_charge,
                                 CoordinationGeometry geometry,
                                 const std::vector<LigandSpec>& ligands);
CoordinationGeometry detect_coordination_geometry(const MolecularSystem& mol,
                                                 int metal_index,
                                                 double tolerance = 15.0);

}  // namespace sbox::chem
