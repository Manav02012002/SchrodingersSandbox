#include "editor/draw_mode.h"

#include "core/covalent_radii.h"

#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>

namespace sbox::editor {

namespace {

Eigen::Vector3f ray_plane_intersect(const Ray& ray, const Eigen::Vector3f& plane_normal, float plane_d) {
    const float denom = plane_normal.dot(ray.direction);
    if (std::abs(denom) <= 1.0e-6f) {
        return ray.origin;
    }
    const float t = -(plane_normal.dot(ray.origin) + plane_d) / denom;
    return ray.origin + t * ray.direction;
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

float highlight_radius(int Z) {
    switch (Z) {
    case 1: return 12.0f;
    case 6: return 16.0f;
    case 7: return 15.0f;
    case 8: return 15.0f;
    default: return 13.0f;
    }
}

sbox::chem::BondOrder cycled_bond_order(sbox::chem::BondOrder current) {
    switch (current) {
    case sbox::chem::BondOrder::Single: return sbox::chem::BondOrder::Double;
    case sbox::chem::BondOrder::Double: return sbox::chem::BondOrder::Triple;
    default: return sbox::chem::BondOrder::Single;
    }
}

}  // namespace

void DrawMode::on_mouse_down(const Ray& ray, int button, bool shift,
                             sbox::chem::MolecularSystem& mol,
                             Selection& selection,
                             CommandStack& commands) {
    (void)shift;

    if (button == 1) {
        drawing_bond_ = false;
        bond_start_atom_ = -1;
        hovered_atom_ = -1;
        return;
    }
    if (button != 0) {
        return;
    }

    const PickResult hit = pick(ray, mol);
    draw_plane_normal_ = ray.direction.normalized();
    const Eigen::Vector3f plane_point =
        mol.num_atoms() > 0 ? Eigen::Vector3f(mol.center_of_mass().cast<float>()) : Eigen::Vector3f::Zero();
    draw_plane_d_ = -draw_plane_normal_.dot(plane_point);

    if (hit.type == PickResult::Type::Atom) {
        drawing_bond_ = true;
        bond_start_atom_ = hit.index;
        hovered_atom_ = hit.index;
        draw_plane_d_ = -draw_plane_normal_.dot(mol.atom(hit.index).position.cast<float>());
        bond_preview_end_ = mol.atom(hit.index).position.cast<float>();
        selection.clear();
        selection.atoms.push_back(hit.index);
        return;
    }

    const Eigen::Vector3f placement = ray_plane_intersect(ray, draw_plane_normal_, draw_plane_d_);
    commands.execute(
        std::make_unique<AddAtomCommand>(sbox::chem::Atom{current_element_, placement.cast<double>(), "", 0}),
        mol);
    selection.clear();
    selection.atoms.push_back(mol.num_atoms() - 1);
    drawing_bond_ = false;
    bond_start_atom_ = -1;
    hovered_atom_ = -1;
}

void DrawMode::on_mouse_up(const Ray& ray, int button,
                           sbox::chem::MolecularSystem& mol,
                           Selection& selection,
                           CommandStack& commands) {
    if (button != 0 || !drawing_bond_ || bond_start_atom_ < 0) {
        return;
    }

    const PickResult hit = pick(ray, mol);
    if (hit.type == PickResult::Type::Atom && hit.index != bond_start_atom_) {
        if (mol.has_bond(bond_start_atom_, hit.index)) {
            for (int i = 0; i < mol.num_bonds(); ++i) {
                const auto& bond = mol.bond(i);
                if ((bond.atom_i == bond_start_atom_ && bond.atom_j == hit.index) ||
                    (bond.atom_i == hit.index && bond.atom_j == bond_start_atom_)) {
                    commands.execute(
                        std::make_unique<ChangeBondOrderCommand>(i, cycled_bond_order(bond.order)),
                        mol);
                    break;
                }
            }
        } else {
            commands.execute(
                std::make_unique<AddBondCommand>(bond_start_atom_, hit.index, current_bond_order_),
                mol);
        }
        selection.clear();
        selection.atoms.push_back(hit.index);
    } else if (hit.type == PickResult::Type::None) {
        Eigen::Vector3f placement = ray_plane_intersect(ray, draw_plane_normal_, draw_plane_d_);
        const float target_length = static_cast<float>(
            sbox::chem::covalent_radius(mol.atom(bond_start_atom_).Z) + sbox::chem::covalent_radius(current_element_));
        Eigen::Vector3f direction = placement - mol.atom(bond_start_atom_).position.cast<float>();
        if (direction.norm() <= 1.0e-6f) {
            direction = Eigen::Vector3f::UnitX();
        } else {
            direction.normalize();
        }
        placement = mol.atom(bond_start_atom_).position.cast<float>() + direction * target_length;

        commands.execute(
            std::make_unique<AddAtomCommand>(sbox::chem::Atom{current_element_, placement.cast<double>(), "", 0}),
            mol);
        const int new_atom_index = mol.num_atoms() - 1;
        commands.execute(
            std::make_unique<AddBondCommand>(bond_start_atom_, new_atom_index, current_bond_order_),
            mol);
        selection.clear();
        selection.atoms.push_back(new_atom_index);
    }

    drawing_bond_ = false;
    bond_start_atom_ = -1;
    hovered_atom_ = -1;
}

void DrawMode::on_mouse_move(const Ray& ray, float dx, float dy, bool dragging,
                             sbox::chem::MolecularSystem& mol,
                             Selection& selection,
                             CommandStack& commands) {
    (void)dx;
    (void)dy;
    (void)dragging;
    (void)selection;
    (void)commands;

    const PickResult hover = pick(ray, mol);
    hovered_atom_ = (hover.type == PickResult::Type::Atom) ? hover.index : -1;

    if (drawing_bond_) {
        if (hovered_atom_ >= 0 && hovered_atom_ != bond_start_atom_) {
            bond_preview_end_ = mol.atom(hovered_atom_).position.cast<float>();
        } else {
            bond_preview_end_ = ray_plane_intersect(ray, draw_plane_normal_, draw_plane_d_);
        }
    } else {
        bond_preview_end_ = ray_plane_intersect(ray, draw_plane_normal_, draw_plane_d_);
    }
}

void DrawMode::on_key(int key, bool ctrl, bool shift,
                      sbox::chem::MolecularSystem& mol,
                      Selection& selection,
                      CommandStack& commands) {
    (void)ctrl;
    (void)shift;
    (void)mol;
    (void)selection;
    (void)commands;

    switch (key) {
    case GLFW_KEY_1: current_element_ = 1; break;
    case GLFW_KEY_6: current_element_ = 6; break;
    case GLFW_KEY_7: current_element_ = 7; break;
    case GLFW_KEY_8: current_element_ = 8; break;
    case GLFW_KEY_9: current_element_ = 9; break;
    case GLFW_KEY_ESCAPE:
        drawing_bond_ = false;
        bond_start_atom_ = -1;
        hovered_atom_ = -1;
        break;
    default:
        break;
    }
}

void DrawMode::draw_overlay(ImDrawList* draw_list,
                            const sbox::chem::MolecularSystem& mol,
                            const Selection& selection,
                            const Eigen::Matrix4f& vp_matrix,
                            const ImVec2& viewport_pos,
                            const ImVec2& viewport_size) {
    (void)selection;

    const ImU32 accent = IM_COL32(45, 185, 185, 255);
    const ImU32 subtle = IM_COL32(120, 180, 180, 180);

    if (drawing_bond_ && bond_start_atom_ >= 0 && bond_start_atom_ < mol.num_atoms()) {
        bool visible_start = false;
        bool visible_end = false;
        const ImVec2 a = world_to_screen(mol.atom(bond_start_atom_).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible_start);
        const ImVec2 b = world_to_screen(bond_preview_end_, vp_matrix, viewport_pos, viewport_size, visible_end);
        if (visible_start && visible_end) {
            const ImVec2 delta(b.x - a.x, b.y - a.y);
            const float length = std::sqrt(delta.x * delta.x + delta.y * delta.y);
            if (length > 1.0f) {
                const ImVec2 dir(delta.x / length, delta.y / length);
                const float dash = 8.0f;
                const float gap = 5.0f;
                float t = 0.0f;
                while (t < length) {
                    const float t1 = std::min(t + dash, length);
                    draw_list->AddLine(
                        ImVec2(a.x + dir.x * t, a.y + dir.y * t),
                        ImVec2(a.x + dir.x * t1, a.y + dir.y * t1),
                        accent,
                        2.0f);
                    t += dash + gap;
                }
            }
        }
        if (hovered_atom_ >= 0 && hovered_atom_ < mol.num_atoms()) {
            bool visible_hover = false;
            const ImVec2 p = world_to_screen(mol.atom(hovered_atom_).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible_hover);
            if (visible_hover) {
                draw_list->AddCircle(p, highlight_radius(mol.atom(hovered_atom_).Z), accent, 0, 2.5f);
            }
        }
        return;
    }

    if (hovered_atom_ >= 0 && hovered_atom_ < mol.num_atoms()) {
        bool visible = false;
        const ImVec2 p = world_to_screen(mol.atom(hovered_atom_).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible);
        if (visible) {
            draw_list->AddCircle(p, highlight_radius(mol.atom(hovered_atom_).Z), subtle, 0, 1.5f);
        }
    }
}

}  // namespace sbox::editor
