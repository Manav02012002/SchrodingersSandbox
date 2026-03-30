#pragma once

#include "core/molecular_system.h"

#include <Eigen/Core>

#include <imgui.h>

#include <vector>

namespace sbox::ui {

void draw_charge_labels(
    const sbox::chem::MolecularSystem& mol,
    const std::vector<double>& charges,
    const Eigen::Matrix4f& view_matrix,
    const Eigen::Matrix4f& proj_matrix,
    const ImVec2& viewport_pos,
    const ImVec2& viewport_size);

void draw_bond_order_labels(
    const sbox::chem::MolecularSystem& mol,
    const Eigen::MatrixXd& mayer_bond_orders,
    const Eigen::Matrix4f& view_matrix,
    const Eigen::Matrix4f& proj_matrix,
    const ImVec2& viewport_pos,
    const ImVec2& viewport_size);

void draw_dipole_arrow(
    const sbox::chem::MolecularSystem& mol,
    const Eigen::Vector3d& dipole_debye,
    const Eigen::Matrix4f& view_matrix,
    const Eigen::Matrix4f& proj_matrix,
    const ImVec2& viewport_pos,
    const ImVec2& viewport_size);

}  // namespace sbox::ui
