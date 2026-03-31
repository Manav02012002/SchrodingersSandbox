#pragma once

#include "core/molecular_system.h"
#include "core/symmetry.h"

#include <Eigen/Core>
#include <imgui.h>

#include <string>
#include <vector>

namespace sbox::ui {

struct SymmetryElement {
    enum class Type { Cn, Sigma, Inversion, Sn };

    Type type = Type::Cn;
    Eigen::Vector3d axis = Eigen::Vector3d::UnitZ();
    Eigen::Vector3d position = Eigen::Vector3d::Zero();
    int order = 0;
    std::string label;
};

std::vector<SymmetryElement> extract_symmetry_elements(
    const sbox::chem::MolecularSystem& mol,
    sbox::chem::PointGroup pg,
    double tolerance = 0.1);

void draw_symmetry_overlays(
    const std::vector<SymmetryElement>& elements,
    const sbox::chem::MolecularSystem& mol,
    const Eigen::Matrix4f& view_matrix,
    const Eigen::Matrix4f& proj_matrix,
    const ImVec2& viewport_pos,
    const ImVec2& viewport_size,
    float molecule_radius);

}  // namespace sbox::ui
