#include "editor/fragment_mode.h"

#include "editor/command.h"

#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>

namespace sbox::editor {

namespace {

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
    if (std::abs(denom) < 1.0e-6f) {
        return ray.origin;
    }
    const float t = -(plane_normal.dot(ray.origin) + plane_d) / denom;
    return ray.origin + std::max(0.0f, t) * ray.direction;
}

}  // namespace

FragmentMode::FragmentMode(const FragmentLibrary* library)
    : library_(library) {}

void FragmentMode::set_fragment(const Fragment* frag) {
    selected_fragment_ = frag;
    show_preview_ = (frag != nullptr);
}

const Fragment* FragmentMode::selected_fragment() const {
    return selected_fragment_;
}

void FragmentMode::on_mouse_down(const Ray& ray, int button, bool shift,
                                 sbox::chem::MolecularSystem& mol,
                                 Selection& selection,
                                 CommandStack& commands) {
    (void)ray;
    (void)shift;
    selection.clear();

    if (button == 1) {
        selected_fragment_ = nullptr;
        show_preview_ = false;
        hover_atom_ = -1;
        return;
    }
    if (button != 0 || selected_fragment_ == nullptr || library_ == nullptr) {
        return;
    }

    const sbox::chem::MolecularSystem placed = library_->place(*selected_fragment_, preview_position_.cast<double>());
    const int first_atom_index = mol.num_atoms();
    commands.execute(std::make_unique<AddFragmentCommand>(placed, selected_fragment_->name), mol);

    if (hover_atom_ >= 0 && selected_fragment_->attachment_atom >= 0) {
        const int attachment_index = first_atom_index + selected_fragment_->attachment_atom;
        if (attachment_index >= 0 && attachment_index < mol.num_atoms() && !mol.has_bond(hover_atom_, attachment_index)) {
            commands.execute(std::make_unique<AddBondCommand>(
                hover_atom_, attachment_index, sbox::chem::BondOrder::Single), mol);
        }
    }
}

void FragmentMode::on_mouse_up(const Ray& ray, int button,
                               sbox::chem::MolecularSystem& mol,
                               Selection& selection,
                               CommandStack& commands) {
    (void)ray;
    (void)button;
    (void)mol;
    (void)selection;
    (void)commands;
}

void FragmentMode::on_mouse_move(const Ray& ray, float dx, float dy, bool dragging,
                                 sbox::chem::MolecularSystem& mol,
                                 Selection& selection,
                                 CommandStack& commands) {
    (void)dx;
    (void)dy;
    (void)dragging;
    (void)selection;
    (void)commands;

    if (selected_fragment_ == nullptr) {
        show_preview_ = false;
        hover_atom_ = -1;
        return;
    }

    const float plane_z = (mol.num_atoms() > 0) ? static_cast<float>(mol.center_of_mass().z()) : 0.0f;
    preview_position_ = ray_plane_intersect(ray, Eigen::Vector3f::UnitZ(), -plane_z);

    const PickResult atom_hit = pick_atom(ray, mol, 1.2f);
    hover_atom_ = -1;
    if (atom_hit.type == PickResult::Type::Atom) {
        hover_atom_ = atom_hit.index;
        preview_position_ = mol.atom(hover_atom_).position.cast<float>();
    }
    show_preview_ = true;
}

void FragmentMode::on_key(int key, bool ctrl, bool shift,
                          sbox::chem::MolecularSystem& mol,
                          Selection& selection,
                          CommandStack& commands) {
    (void)ctrl;
    (void)shift;
    (void)mol;
    (void)selection;
    (void)commands;

    if (key == GLFW_KEY_ESCAPE) {
        selected_fragment_ = nullptr;
        show_preview_ = false;
        hover_atom_ = -1;
    }
}

void FragmentMode::draw_overlay(ImDrawList* draw_list,
                                const sbox::chem::MolecularSystem& mol,
                                const Selection& selection,
                                const Eigen::Matrix4f& vp_matrix,
                                const ImVec2& viewport_pos,
                                const ImVec2& viewport_size) {
    (void)selection;
    if (selected_fragment_ == nullptr || library_ == nullptr || !show_preview_) {
        return;
    }

    const sbox::chem::MolecularSystem placed = library_->place(*selected_fragment_, preview_position_.cast<double>());
    const ImU32 atom_color = IM_COL32(120, 220, 255, 170);
    const ImU32 bond_color = IM_COL32(120, 220, 255, 120);
    const ImU32 attach_color = IM_COL32(60, 240, 180, 255);

    for (int i = 0; i < placed.num_bonds(); ++i) {
        const auto& bond = placed.bond(i);
        bool va = false;
        bool vb = false;
        const ImVec2 a = world_to_screen(placed.atom(bond.atom_i).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, va);
        const ImVec2 b = world_to_screen(placed.atom(bond.atom_j).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, vb);
        if (va && vb) {
            draw_list->AddLine(a, b, bond_color, 2.0f);
        }
    }

    for (int i = 0; i < placed.num_atoms(); ++i) {
        bool visible = false;
        const ImVec2 p = world_to_screen(placed.atom(i).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible);
        if (!visible) {
            continue;
        }
        const float radius = (i == selected_fragment_->attachment_atom) ? 8.0f : 6.0f;
        const ImU32 color = (i == selected_fragment_->attachment_atom) ? attach_color : atom_color;
        draw_list->AddCircleFilled(p, radius, color);
        draw_list->AddCircle(p, radius + 1.5f, IM_COL32(255, 255, 255, 120), 0, 1.5f);
    }

    if (hover_atom_ >= 0 && hover_atom_ < mol.num_atoms()) {
        bool visible = false;
        const ImVec2 p = world_to_screen(mol.atom(hover_atom_).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible);
        if (visible) {
            draw_list->AddCircle(p, 14.0f, attach_color, 0, 2.5f);
        }
    }
}

}  // namespace sbox::editor
