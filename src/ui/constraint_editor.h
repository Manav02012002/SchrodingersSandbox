#pragma once

#include "core/molecular_system.h"
#include "editor/picking.h"
#include "ui/app_state.h"

#include <Eigen/Core>
#include <imgui.h>

namespace sbox::ui {

void draw_constraint_editor(AppState& state,
                            const sbox::chem::MolecularSystem& mol,
                            const sbox::editor::Selection& selection);

void draw_constraint_overlays(const std::vector<AppState::GeometricConstraint>& constraints,
                              const sbox::chem::MolecularSystem& mol,
                              const Eigen::Matrix4f& vp_matrix,
                              const ImVec2& viewport_pos,
                              const ImVec2& viewport_size);

}  // namespace sbox::ui
