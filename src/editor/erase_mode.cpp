#include "editor/erase_mode.h"

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

float highlight_radius(int z) {
    switch (z) {
    case 1: return 12.0f;
    case 6: return 16.0f;
    case 7:
    case 8: return 15.0f;
    default: return 13.0f;
    }
}

}  // namespace

void EraseMode::on_mouse_down(const Ray& ray, int button, bool shift,
                              sbox::chem::MolecularSystem& mol,
                              Selection& selection,
                              CommandStack& commands) {
    (void)shift;
    if (button != 0) {
        return;
    }

    const PickResult hit = pick(ray, mol);
    if (hit.type == PickResult::Type::Atom) {
        commands.execute(std::make_unique<RemoveAtomCommand>(hit.index), mol);
        selection.clear();
    } else if (hit.type == PickResult::Type::Bond) {
        commands.execute(std::make_unique<RemoveBondCommand>(hit.index), mol);
        selection.clear();
    }
}

void EraseMode::on_mouse_up(const Ray& ray, int button,
                            sbox::chem::MolecularSystem& mol,
                            Selection& selection,
                            CommandStack& commands) {
    (void)ray;
    (void)button;
    (void)mol;
    (void)selection;
    (void)commands;
}

void EraseMode::on_mouse_move(const Ray& ray, float dx, float dy, bool dragging,
                              sbox::chem::MolecularSystem& mol,
                              Selection& selection,
                              CommandStack& commands) {
    (void)dx;
    (void)dy;
    (void)dragging;
    (void)selection;
    (void)commands;

    const PickResult hit = pick(ray, mol);
    hover_atom_ = (hit.type == PickResult::Type::Atom) ? hit.index : -1;
    hover_bond_ = (hit.type == PickResult::Type::Bond) ? hit.index : -1;
}

void EraseMode::on_key(int key, bool ctrl, bool shift,
                       sbox::chem::MolecularSystem& mol,
                       Selection& selection,
                       CommandStack& commands) {
    (void)key;
    (void)ctrl;
    (void)shift;
    (void)mol;
    (void)selection;
    (void)commands;
}

void EraseMode::draw_overlay(ImDrawList* draw_list,
                             const sbox::chem::MolecularSystem& mol,
                             const Selection& selection,
                             const Eigen::Matrix4f& vp_matrix,
                             const ImVec2& viewport_pos,
                             const ImVec2& viewport_size) {
    (void)selection;
    const ImU32 danger = IM_COL32(235, 70, 70, 255);

    if (hover_bond_ >= 0 && hover_bond_ < mol.num_bonds()) {
        const auto& bond = mol.bond(hover_bond_);
        bool visible_a = false;
        bool visible_b = false;
        const ImVec2 a = world_to_screen(mol.atom(bond.atom_i).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible_a);
        const ImVec2 b = world_to_screen(mol.atom(bond.atom_j).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible_b);
        if (visible_a && visible_b) {
            draw_list->AddLine(a, b, danger, 4.0f);
        }
    }

    if (hover_atom_ >= 0 && hover_atom_ < mol.num_atoms()) {
        bool visible = false;
        const ImVec2 p = world_to_screen(mol.atom(hover_atom_).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible);
        if (visible) {
            const float r = highlight_radius(mol.atom(hover_atom_).Z);
            draw_list->AddCircle(p, r, danger, 0, 2.5f);
            draw_list->AddLine(ImVec2(p.x - r * 0.5f, p.y - r * 0.5f), ImVec2(p.x + r * 0.5f, p.y + r * 0.5f), danger, 2.0f);
            draw_list->AddLine(ImVec2(p.x - r * 0.5f, p.y + r * 0.5f), ImVec2(p.x + r * 0.5f, p.y - r * 0.5f), danger, 2.0f);
        }
    }
}

}  // namespace sbox::editor
