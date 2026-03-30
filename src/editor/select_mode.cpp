#include "editor/select_mode.h"

#include "ui/context_menu.h"

#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>

namespace sbox::editor {

namespace {

float atom_render_radius_for_pick(int Z) {
    switch (Z) {
    case 1: return 0.6f;
    case 6: return 1.0f;
    case 7: return 0.95f;
    case 8: return 0.9f;
    case 16: return 1.2f;
    default: return 0.8f;
    }
}

ImVec2 world_to_screen(const Eigen::Vector3f& world_pos,
                       const Eigen::Matrix4f& vp_matrix,
                       const ImVec2& viewport_pos,
                       const ImVec2& viewport_size,
                       bool& visible) {
    const Eigen::Vector4f clip = vp_matrix * Eigen::Vector4f(world_pos.x(), world_pos.y(), world_pos.z(), 1.0f);
    visible = false;
    if (clip.w() <= 0.0f) {
        return {};
    }
    const Eigen::Vector3f ndc = clip.head<3>() / clip.w();
    if (ndc.x() < -1.0f || ndc.x() > 1.0f || ndc.y() < -1.0f || ndc.y() > 1.0f || ndc.z() < -1.0f || ndc.z() > 1.0f) {
        return {};
    }
    visible = true;
    return ImVec2(
        viewport_pos.x + (ndc.x() * 0.5f + 0.5f) * viewport_size.x,
        viewport_pos.y + (1.0f - (ndc.y() * 0.5f + 0.5f)) * viewport_size.y);
}

Eigen::Vector3f ray_plane_intersect(const Ray& ray, const Eigen::Vector3f& plane_normal, float plane_d) {
    const float denom = plane_normal.dot(ray.direction);
    if (std::abs(denom) <= 1.0e-6f) {
        return ray.origin;
    }
    const float t = -(plane_normal.dot(ray.origin) + plane_d) / denom;
    return ray.origin + t * ray.direction;
}

bool positions_changed(const std::vector<Eigen::Vector3d>& a, const std::vector<Eigen::Vector3d>& b) {
    if (a.size() != b.size()) {
        return true;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if ((a[i] - b[i]).norm() > 0.01) {
            return true;
        }
    }
    return false;
}

}  // namespace

void SelectMode::on_mouse_down(const Ray& ray, int button, bool shift,
                               sbox::chem::MolecularSystem& mol,
                               Selection& selection,
                               CommandStack& commands) {
    (void)commands;

    if (button == 0) {
        const PickResult hit = pick(ray, mol);
        dragging_ = true;
        drag_start_ray_ = ray;

        if (hit.type == PickResult::Type::None) {
            if (!shift) {
                selection.clear();
            }
            is_drag_move_ = false;
            drag_atom_indices_.clear();
            drag_original_positions_.clear();
            return;
        }

        if (shift) {
            if (hit.type == PickResult::Type::Atom) {
                selection.toggle_atom(hit.index);
            } else if (hit.type == PickResult::Type::Bond) {
                selection.toggle_bond(hit.index);
            }
            is_drag_move_ = false;
            drag_atom_indices_.clear();
            drag_original_positions_.clear();
            return;
        }

        const bool already_selected =
            (hit.type == PickResult::Type::Atom && selection.has_atom(hit.index)) ||
            (hit.type == PickResult::Type::Bond && selection.has_bond(hit.index));

        if (already_selected && hit.type == PickResult::Type::Atom) {
            is_drag_move_ = true;
            drag_atom_indices_ = selection.atoms;
            if (drag_atom_indices_.empty()) {
                drag_atom_indices_.push_back(hit.index);
            }
            drag_original_positions_.clear();
            drag_original_positions_.reserve(drag_atom_indices_.size());
            for (int atom_index : drag_atom_indices_) {
                drag_original_positions_.push_back(mol.atom(atom_index).position);
            }
            drag_plane_normal_ = ray.direction.normalized();
            drag_plane_d_ = -drag_plane_normal_.dot(mol.atom(hit.index).position.cast<float>());
            drag_start_ = ray_plane_intersect(ray, drag_plane_normal_, drag_plane_d_);
        } else {
            selection.clear();
            if (hit.type == PickResult::Type::Atom) {
                selection.atoms.push_back(hit.index);
            } else if (hit.type == PickResult::Type::Bond) {
                selection.bonds.push_back(hit.index);
            }
            is_drag_move_ = false;
            drag_atom_indices_.clear();
            drag_original_positions_.clear();
        }
        return;
    }

    if (button == 1) {
        if (context_menu_state_ != nullptr) {
            const PickResult hit = pick(ray, mol);
            context_menu_state_->show = true;
            context_menu_state_->position = ImGui::GetMousePos();
            context_menu_state_->clicked_atom = hit.type == PickResult::Type::Atom ? hit.index : -1;
            context_menu_state_->clicked_bond = hit.type == PickResult::Type::Bond ? hit.index : -1;
        }
    }
}

void SelectMode::on_mouse_up(const Ray& ray, int button,
                             sbox::chem::MolecularSystem& mol,
                             Selection& selection,
                             CommandStack& commands) {
    (void)selection;

    if (button != 0) {
        return;
    }

    if (is_drag_move_ && !drag_atom_indices_.empty()) {
        std::vector<Eigen::Vector3d> final_positions;
        final_positions.reserve(drag_atom_indices_.size());
        for (int atom_index : drag_atom_indices_) {
            final_positions.push_back(mol.atom(atom_index).position);
        }

        if (positions_changed(drag_original_positions_, final_positions)) {
            for (std::size_t i = 0; i < drag_atom_indices_.size(); ++i) {
                mol.atom(drag_atom_indices_[i]).position = drag_original_positions_[i];
            }
            commands.execute(
                std::make_unique<MoveAtomsCommand>(drag_atom_indices_, final_positions),
                mol);
        }
    }

    dragging_ = false;
    is_drag_move_ = false;
    drag_atom_indices_.clear();
    drag_original_positions_.clear();
    drag_start_ = ray.origin;
}

void SelectMode::on_mouse_move(const Ray& ray, float dx, float dy, bool dragging,
                               sbox::chem::MolecularSystem& mol,
                               Selection& selection,
                               CommandStack& commands) {
    (void)dx;
    (void)dy;
    (void)selection;
    (void)commands;

    if (!dragging || !dragging_ || !is_drag_move_ || drag_atom_indices_.empty()) {
        return;
    }

    const Eigen::Vector3f current = ray_plane_intersect(ray, drag_plane_normal_, drag_plane_d_);
    const Eigen::Vector3f displacement = current - drag_start_;
    for (std::size_t i = 0; i < drag_atom_indices_.size(); ++i) {
        mol.atom(drag_atom_indices_[i]).position = drag_original_positions_[i] + displacement.cast<double>();
    }
}

void SelectMode::on_key(int key, bool ctrl, bool shift,
                        sbox::chem::MolecularSystem& mol,
                        Selection& selection,
                        CommandStack& commands) {
    (void)shift;

    if ((key == GLFW_KEY_DELETE || key == GLFW_KEY_BACKSPACE)) {
        std::vector<int> atom_indices = selection.atoms;
        std::sort(atom_indices.begin(), atom_indices.end(), std::greater<int>());
        for (int atom_index : atom_indices) {
            commands.execute(std::make_unique<RemoveAtomCommand>(atom_index), mol);
        }

        if (atom_indices.empty()) {
            std::vector<int> bond_indices = selection.bonds;
            std::sort(bond_indices.begin(), bond_indices.end(), std::greater<int>());
            for (int bond_index : bond_indices) {
                commands.execute(std::make_unique<RemoveBondCommand>(bond_index), mol);
            }
        }

        selection.clear();
        return;
    }

    if (ctrl && key == GLFW_KEY_A) {
        selection.clear();
        selection.atoms.reserve(static_cast<std::size_t>(mol.num_atoms()));
        for (int i = 0; i < mol.num_atoms(); ++i) {
            selection.atoms.push_back(i);
        }
        return;
    }

    if (key == GLFW_KEY_ESCAPE) {
        selection.clear();
    }
}

void SelectMode::draw_overlay(ImDrawList* draw_list,
                              const sbox::chem::MolecularSystem& mol,
                              const Selection& selection,
                              const Eigen::Matrix4f& vp_matrix,
                              const ImVec2& viewport_pos,
                              const ImVec2& viewport_size) {
    const ImU32 accent = IM_COL32(45, 185, 185, 255);

    for (int atom_index : selection.atoms) {
        if (atom_index < 0 || atom_index >= mol.num_atoms()) {
            continue;
        }
        bool visible_center = false;
        const Eigen::Vector3f center = mol.atom(atom_index).position.cast<float>();
        const ImVec2 screen_center = world_to_screen(center, vp_matrix, viewport_pos, viewport_size, visible_center);
        if (!visible_center) {
            continue;
        }

        bool visible_radius = false;
        const ImVec2 screen_radius = world_to_screen(
            center + Eigen::Vector3f(atom_render_radius_for_pick(mol.atom(atom_index).Z), 0.0f, 0.0f),
            vp_matrix,
            viewport_pos,
            viewport_size,
            visible_radius);
        float radius = 12.0f;
        if (visible_radius) {
            radius = std::max(10.0f, std::sqrt((screen_radius.x - screen_center.x) * (screen_radius.x - screen_center.x) +
                                               (screen_radius.y - screen_center.y) * (screen_radius.y - screen_center.y)));
        }

        draw_list->AddCircle(screen_center, radius + 4.0f, accent, 0, 2.5f);
    }

    for (int bond_index : selection.bonds) {
        if (bond_index < 0 || bond_index >= mol.num_bonds()) {
            continue;
        }
        const sbox::chem::Bond& bond = mol.bond(bond_index);
        bool visible_a = false;
        bool visible_b = false;
        const ImVec2 a = world_to_screen(mol.atom(bond.atom_i).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible_a);
        const ImVec2 b = world_to_screen(mol.atom(bond.atom_j).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible_b);
        if (visible_a && visible_b) {
            draw_list->AddLine(a, b, accent, 4.0f);
        }
    }

}

}  // namespace sbox::editor
