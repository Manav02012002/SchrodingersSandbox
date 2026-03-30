#include "ui/constraint_editor.h"

#include "core/elements.h"

#include <imgui.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>

namespace sbox::ui {

namespace {

constexpr double kBohrToAngstrom = 0.529177;
constexpr double kRadToDeg = 57.29577951308232;

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

std::string atom_label(const sbox::chem::MolecularSystem& mol, int atom_index) {
    return std::string(sbox::elements::get_element(mol.atom(atom_index).Z).symbol) + std::to_string(atom_index + 1);
}

std::string constraint_type_label(AppState::GeometricConstraint::Type type) {
    switch (type) {
    case AppState::GeometricConstraint::Type::FreezeAtom: return "Freeze";
    case AppState::GeometricConstraint::Type::FixDistance: return "Distance";
    case AppState::GeometricConstraint::Type::FixAngle: return "Angle";
    case AppState::GeometricConstraint::Type::FixDihedral: return "Dihedral";
    }
    return "Unknown";
}

std::string constraint_atoms_label(const sbox::chem::MolecularSystem& mol, const AppState::GeometricConstraint& c) {
    std::string text;
    for (std::size_t i = 0; i < c.atom_indices.size(); ++i) {
        if (i > 0) {
            text += "-";
        }
        text += atom_label(mol, c.atom_indices[i]);
    }
    return text;
}

std::string constraint_value_label(const AppState::GeometricConstraint& c) {
    char buffer[64];
    switch (c.type) {
    case AppState::GeometricConstraint::Type::FreezeAtom:
        return "-";
    case AppState::GeometricConstraint::Type::FixDistance:
        std::snprintf(buffer, sizeof(buffer), "%.3f A", c.value * kBohrToAngstrom);
        return buffer;
    case AppState::GeometricConstraint::Type::FixAngle:
    case AppState::GeometricConstraint::Type::FixDihedral:
        std::snprintf(buffer, sizeof(buffer), "%.2f deg", c.value * kRadToDeg);
        return buffer;
    }
    return {};
}

void add_constraint(std::vector<AppState::GeometricConstraint>& constraints,
                    AppState::GeometricConstraint::Type type,
                    std::vector<int> atom_indices,
                    double value = 0.0) {
    constraints.push_back({type, std::move(atom_indices), value, true});
}

void draw_dashed_line(ImDrawList* draw, const ImVec2& a, const ImVec2& b, ImU32 color, float thickness = 2.0f) {
    const ImVec2 delta(b.x - a.x, b.y - a.y);
    const float length = std::sqrt(delta.x * delta.x + delta.y * delta.y);
    if (length < 1.0f) {
        return;
    }
    const ImVec2 dir(delta.x / length, delta.y / length);
    constexpr float dash = 6.0f;
    for (float t = 0.0f; t < length; t += dash * 2.0f) {
        const float t1 = std::min(length, t + dash);
        draw->AddLine(ImVec2(a.x + dir.x * t, a.y + dir.y * t),
                      ImVec2(a.x + dir.x * t1, a.y + dir.y * t1),
                      color, thickness);
    }
}

}  // namespace

void draw_constraint_editor(AppState& state,
                            const sbox::chem::MolecularSystem& mol,
                            const sbox::editor::Selection& selection) {
    if (!ImGui::Begin("Constraints")) {
        ImGui::End();
        return;
    }

    static float target_distance_ang = 0.0f;
    static float target_angle_deg = 0.0f;
    static float target_dihedral_deg = 0.0f;

    ImGui::TextUnformatted("Add Constraints from Selection");
    const auto& atoms = selection.atoms;
    if (atoms.size() == 1) {
        if (ImGui::Button("Freeze Atom")) {
            add_constraint(state.constraints, AppState::GeometricConstraint::Type::FreezeAtom, {atoms[0]});
        }
    } else if (atoms.size() == 2) {
        const double current_distance_ang = mol.distance(atoms[0], atoms[1]) * kBohrToAngstrom;
        target_distance_ang = target_distance_ang <= 0.0f ? static_cast<float>(current_distance_ang) : target_distance_ang;
        ImGui::Text("d(%s-%s) = %.3f A",
                    atom_label(mol, atoms[0]).c_str(), atom_label(mol, atoms[1]).c_str(), current_distance_ang);
        ImGui::InputFloat("Target distance (A)", &target_distance_ang);
        if (ImGui::Button("Fix Distance")) {
            add_constraint(state.constraints, AppState::GeometricConstraint::Type::FixDistance,
                           {atoms[0], atoms[1]}, target_distance_ang / static_cast<float>(kBohrToAngstrom));
        }
        ImGui::SameLine();
        if (ImGui::Button("Freeze Both Atoms")) {
            add_constraint(state.constraints, AppState::GeometricConstraint::Type::FreezeAtom, {atoms[0]});
            add_constraint(state.constraints, AppState::GeometricConstraint::Type::FreezeAtom, {atoms[1]});
        }
    } else if (atoms.size() == 3) {
        const double current_angle_deg = mol.angle(atoms[0], atoms[1], atoms[2]) * kRadToDeg;
        target_angle_deg = target_angle_deg == 0.0f ? static_cast<float>(current_angle_deg) : target_angle_deg;
        ImGui::Text("Angle = %.2f deg", current_angle_deg);
        ImGui::InputFloat("Target angle (deg)", &target_angle_deg);
        if (ImGui::Button("Fix Angle")) {
            add_constraint(state.constraints, AppState::GeometricConstraint::Type::FixAngle,
                           {atoms[0], atoms[1], atoms[2]}, target_angle_deg / static_cast<float>(kRadToDeg));
        }
    } else if (atoms.size() == 4) {
        const double current_dihedral_deg = mol.dihedral(atoms[0], atoms[1], atoms[2], atoms[3]) * kRadToDeg;
        target_dihedral_deg = target_dihedral_deg == 0.0f ? static_cast<float>(current_dihedral_deg) : target_dihedral_deg;
        ImGui::Text("Dihedral = %.2f deg", current_dihedral_deg);
        ImGui::InputFloat("Target dihedral (deg)", &target_dihedral_deg);
        if (ImGui::Button("Fix Dihedral")) {
            add_constraint(state.constraints, AppState::GeometricConstraint::Type::FixDihedral,
                           {atoms[0], atoms[1], atoms[2], atoms[3]}, target_dihedral_deg / static_cast<float>(kRadToDeg));
        }
    } else {
        ImGui::TextDisabled("Select 1-4 atoms to create a constraint.");
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Active Constraints");
    if (ImGui::BeginTable("ConstraintsTable", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("Active", ImGuiTableColumnFlags_WidthFixed, 56.0f);
        ImGui::TableSetupColumn("Type");
        ImGui::TableSetupColumn("Atoms");
        ImGui::TableSetupColumn("Value");
        ImGui::TableSetupColumn("Actions", ImGuiTableColumnFlags_WidthFixed, 60.0f);
        ImGui::TableHeadersRow();

        int remove_index = -1;
        for (int i = 0; i < static_cast<int>(state.constraints.size()); ++i) {
            auto& c = state.constraints[static_cast<std::size_t>(i)];
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Checkbox(("##active_" + std::to_string(i)).c_str(), &c.active);
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%s", constraint_type_label(c.type).c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%s", constraint_atoms_label(mol, c).c_str());
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("%s", constraint_value_label(c).c_str());
            ImGui::TableSetColumnIndex(4);
            if (ImGui::Button(("x##constraint_" + std::to_string(i)).c_str())) {
                remove_index = i;
            }
        }
        if (remove_index >= 0) {
            state.constraints.erase(state.constraints.begin() + remove_index);
        }
        ImGui::EndTable();
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Quick Actions");
    if (ImGui::Button("Freeze All Hydrogens")) {
        for (int i = 0; i < mol.num_atoms(); ++i) {
            if (mol.atom(i).Z == 1) {
                add_constraint(state.constraints, AppState::GeometricConstraint::Type::FreezeAtom, {i});
            }
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Freeze Heavy Atoms")) {
        for (int i = 0; i < mol.num_atoms(); ++i) {
            if (mol.atom(i).Z != 1) {
                add_constraint(state.constraints, AppState::GeometricConstraint::Type::FreezeAtom, {i});
            }
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Clear All")) {
        state.constraints.clear();
    }

    ImGui::End();
}

void draw_constraint_overlays(const std::vector<AppState::GeometricConstraint>& constraints,
                              const sbox::chem::MolecularSystem& mol,
                              const Eigen::Matrix4f& vp_matrix,
                              const ImVec2& viewport_pos,
                              const ImVec2& viewport_size) {
    ImDrawList* draw = ImGui::GetForegroundDrawList();
    char buffer[64];
    for (const auto& c : constraints) {
        if (!c.active || c.atom_indices.empty()) {
            continue;
        }
        if (c.type == AppState::GeometricConstraint::Type::FreezeAtom && c.atom_indices[0] < mol.num_atoms()) {
            bool visible = false;
            const ImVec2 pos = world_to_screen(mol.atom(c.atom_indices[0]).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible);
            if (visible) {
                draw->AddCircle(pos, 10.0f, IM_COL32(230, 80, 80, 255), 0, 2.0f);
                draw->AddText(ImVec2(pos.x - 4.0f, pos.y - 8.0f), IM_COL32(230, 80, 80, 255), "L");
            }
        } else if (c.type == AppState::GeometricConstraint::Type::FixDistance && c.atom_indices.size() >= 2) {
            bool v0 = false;
            bool v1 = false;
            const ImVec2 p0 = world_to_screen(mol.atom(c.atom_indices[0]).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, v0);
            const ImVec2 p1 = world_to_screen(mol.atom(c.atom_indices[1]).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, v1);
            if (v0 && v1) {
                draw_dashed_line(draw, p0, p1, IM_COL32(240, 210, 70, 255), 2.0f);
                const ImVec2 mid((p0.x + p1.x) * 0.5f, (p0.y + p1.y) * 0.5f);
                std::snprintf(buffer, sizeof(buffer), "%.2f A", c.value * kBohrToAngstrom);
                draw->AddText(ImVec2(mid.x + 4.0f, mid.y - 14.0f), IM_COL32(240, 210, 70, 255), buffer);
            }
        } else if (c.type == AppState::GeometricConstraint::Type::FixAngle && c.atom_indices.size() >= 3) {
            bool v0 = false;
            bool v1 = false;
            bool v2 = false;
            const ImVec2 p0 = world_to_screen(mol.atom(c.atom_indices[0]).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, v0);
            const ImVec2 p1 = world_to_screen(mol.atom(c.atom_indices[1]).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, v1);
            const ImVec2 p2 = world_to_screen(mol.atom(c.atom_indices[2]).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, v2);
            if (v0 && v1 && v2) {
                draw->AddLine(p1, p0, IM_COL32(240, 210, 70, 255), 2.0f);
                draw->AddLine(p1, p2, IM_COL32(240, 210, 70, 255), 2.0f);
                const float radius = 22.0f;
                const ImVec2 v01(p0.x - p1.x, p0.y - p1.y);
                const ImVec2 v21(p2.x - p1.x, p2.y - p1.y);
                const float a0 = std::atan2(v01.y, v01.x);
                const float a1 = std::atan2(v21.y, v21.x);
                const int segments = 20;
                for (int s = 0; s < segments; ++s) {
                    const float t0 = static_cast<float>(s) / segments;
                    const float t1 = static_cast<float>(s + 1) / segments;
                    const float ang0 = a0 + (a1 - a0) * t0;
                    const float ang1 = a0 + (a1 - a0) * t1;
                    draw->AddLine(ImVec2(p1.x + std::cos(ang0) * radius, p1.y + std::sin(ang0) * radius),
                                  ImVec2(p1.x + std::cos(ang1) * radius, p1.y + std::sin(ang1) * radius),
                                  IM_COL32(240, 210, 70, 255), 2.0f);
                }
                std::snprintf(buffer, sizeof(buffer), "%.1f deg", c.value * kRadToDeg);
                draw->AddText(ImVec2(p1.x + radius + 4.0f, p1.y - 8.0f), IM_COL32(240, 210, 70, 255), buffer);
            }
        } else if (c.type == AppState::GeometricConstraint::Type::FixDihedral && c.atom_indices.size() >= 4) {
            std::array<ImVec2, 4> pts{};
            std::array<bool, 4> vis{};
            bool all_visible = true;
            for (int i = 0; i < 4; ++i) {
                pts[i] = world_to_screen(mol.atom(c.atom_indices[i]).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, vis[i]);
                all_visible = all_visible && vis[i];
            }
            if (all_visible) {
                draw_dashed_line(draw, pts[0], pts[1], IM_COL32(240, 210, 70, 255), 2.0f);
                draw_dashed_line(draw, pts[1], pts[2], IM_COL32(240, 210, 70, 255), 2.0f);
                draw_dashed_line(draw, pts[2], pts[3], IM_COL32(240, 210, 70, 255), 2.0f);
                const ImVec2 mid((pts[1].x + pts[2].x) * 0.5f, (pts[1].y + pts[2].y) * 0.5f);
                std::snprintf(buffer, sizeof(buffer), "%.1f deg", c.value * kRadToDeg);
                draw->AddText(ImVec2(mid.x + 4.0f, mid.y - 14.0f), IM_COL32(240, 210, 70, 255), buffer);
            }
        }
    }
}

}  // namespace sbox::ui
