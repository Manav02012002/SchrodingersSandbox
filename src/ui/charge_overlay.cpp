#include "ui/charge_overlay.h"

#include <imgui.h>

#include <Eigen/Core>

#include <cmath>
#include <cstdio>

namespace sbox::ui {

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

    const float x = viewport_pos.x + (ndc.x() * 0.5f + 0.5f) * viewport_size.x;
    const float y = viewport_pos.y + (1.0f - (ndc.y() * 0.5f + 0.5f)) * viewport_size.y;
    visible = x >= viewport_pos.x && x <= viewport_pos.x + viewport_size.x &&
              y >= viewport_pos.y && y <= viewport_pos.y + viewport_size.y;
    return ImVec2(x, y);
}

ImU32 charge_color(double charge) {
    if (charge > 0.01) {
        return IM_COL32(235, 90, 90, 255);
    }
    if (charge < -0.01) {
        return IM_COL32(90, 150, 255, 255);
    }
    return IM_COL32(180, 180, 180, 255);
}

ImU32 bond_order_color(double order) {
    if (std::abs(order - 1.5) < 0.25) {
        return IM_COL32(235, 165, 65, 255);
    }
    if (std::abs(order - 3.0) < 0.4) {
        return IM_COL32(90, 150, 255, 255);
    }
    if (std::abs(order - 2.0) < 0.35) {
        return IM_COL32(90, 220, 130, 255);
    }
    return IM_COL32(190, 190, 190, 255);
}

void draw_centered_text(ImDrawList* draw_list, const ImVec2& anchor, ImU32 color, const char* text) {
    const ImVec2 size = ImGui::CalcTextSize(text);
    draw_list->AddText(ImVec2(anchor.x - 0.5f * size.x, anchor.y - 0.5f * size.y), color, text);
}

}  // namespace

void draw_charge_labels(const sbox::chem::MolecularSystem& mol,
                        const std::vector<double>& charges,
                        const Eigen::Matrix4f& view_matrix,
                        const Eigen::Matrix4f& proj_matrix,
                        const ImVec2& viewport_pos,
                        const ImVec2& viewport_size) {
    if (charges.size() < static_cast<std::size_t>(mol.num_atoms())) {
        return;
    }

    ImDrawList* draw_list = ImGui::GetForegroundDrawList();
    const Eigen::Matrix4f vp = proj_matrix * view_matrix;

    char buffer[32];
    for (int i = 0; i < mol.num_atoms(); ++i) {
        bool visible = false;
        const Eigen::Vector3d pos = mol.atom(i).position + Eigen::Vector3d(0.0, 0.35, 0.0);
        const ImVec2 screen = world_to_screen(pos.cast<float>(), vp, viewport_pos, viewport_size, visible);
        if (!visible) {
            continue;
        }
        std::snprintf(buffer, sizeof(buffer), "%+.2f", charges[static_cast<std::size_t>(i)]);
        draw_centered_text(draw_list, screen, charge_color(charges[static_cast<std::size_t>(i)]), buffer);
    }
}

void draw_bond_order_labels(const sbox::chem::MolecularSystem& mol,
                            const Eigen::MatrixXd& mayer_bond_orders,
                            const Eigen::Matrix4f& view_matrix,
                            const Eigen::Matrix4f& proj_matrix,
                            const ImVec2& viewport_pos,
                            const ImVec2& viewport_size) {
    if (mayer_bond_orders.rows() < mol.num_atoms() || mayer_bond_orders.cols() < mol.num_atoms()) {
        return;
    }

    ImDrawList* draw_list = ImGui::GetForegroundDrawList();
    const Eigen::Matrix4f vp = proj_matrix * view_matrix;

    char buffer[32];
    for (const sbox::chem::Bond& bond : mol.bonds()) {
        const Eigen::Vector3d midpoint = 0.5 * (mol.atom(bond.atom_i).position + mol.atom(bond.atom_j).position);
        bool visible = false;
        const ImVec2 screen = world_to_screen(midpoint.cast<float>(), vp, viewport_pos, viewport_size, visible);
        if (!visible) {
            continue;
        }
        const double order = mayer_bond_orders(bond.atom_i, bond.atom_j);
        std::snprintf(buffer, sizeof(buffer), "%.1f", order);
        draw_centered_text(draw_list, screen, bond_order_color(order), buffer);
    }
}

void draw_dipole_arrow(const sbox::chem::MolecularSystem& mol,
                       const Eigen::Vector3d& dipole_debye,
                       const Eigen::Matrix4f& view_matrix,
                       const Eigen::Matrix4f& proj_matrix,
                       const ImVec2& viewport_pos,
                       const ImVec2& viewport_size) {
    if (mol.num_atoms() == 0 || dipole_debye.norm() < 1.0e-6) {
        return;
    }

    double total_Z = 0.0;
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (const sbox::chem::Atom& atom : mol.atoms()) {
        total_Z += static_cast<double>(atom.Z);
        center += static_cast<double>(atom.Z) * atom.position;
    }
    if (total_Z <= 0.0) {
        return;
    }
    center /= total_Z;

    const Eigen::Vector3d end = center + dipole_debye * 2.0;

    bool visible0 = false;
    bool visible1 = false;
    const Eigen::Matrix4f vp = proj_matrix * view_matrix;
    const ImVec2 start_screen = world_to_screen(center.cast<float>(), vp, viewport_pos, viewport_size, visible0);
    const ImVec2 end_screen = world_to_screen(end.cast<float>(), vp, viewport_pos, viewport_size, visible1);
    if (!visible0 || !visible1) {
        return;
    }

    ImDrawList* draw_list = ImGui::GetForegroundDrawList();
    const ImU32 color = IM_COL32(240, 220, 70, 255);
    draw_list->AddLine(start_screen, end_screen, color, 3.0f);

    const ImVec2 dir(end_screen.x - start_screen.x, end_screen.y - start_screen.y);
    const float len = std::sqrt(dir.x * dir.x + dir.y * dir.y);
    if (len > 1.0f) {
        const ImVec2 unit(dir.x / len, dir.y / len);
        const ImVec2 perp(-unit.y, unit.x);
        const ImVec2 head_base(end_screen.x - unit.x * 12.0f, end_screen.y - unit.y * 12.0f);
        draw_list->AddTriangleFilled(end_screen,
                                     ImVec2(head_base.x + perp.x * 6.0f, head_base.y + perp.y * 6.0f),
                                     ImVec2(head_base.x - perp.x * 6.0f, head_base.y - perp.y * 6.0f),
                                     color);
    }

    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "|u| = %.2f D", dipole_debye.norm());
    const ImVec2 label_pos(0.5f * (start_screen.x + end_screen.x), 0.5f * (start_screen.y + end_screen.y) - 14.0f);
    draw_centered_text(draw_list, label_pos, color, buffer);
}

}  // namespace sbox::ui
