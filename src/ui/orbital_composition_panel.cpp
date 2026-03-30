#include "ui/orbital_composition_panel.h"

#include "analysis/orbital_composition.h"
#include "ui/plot_utils.h"

#include <imgui.h>

#include <algorithm>
#include <cstdio>
#include <string>

namespace sbox::ui {

namespace {

ImVec4 cpk_color_imvec4(int Z) {
    switch (Z) {
    case 1: return ImVec4(1.0f, 1.0f, 1.0f, 1.0f);
    case 6: return ImVec4(0.56f, 0.56f, 0.56f, 1.0f);
    case 7: return ImVec4(0.19f, 0.31f, 0.97f, 1.0f);
    case 8: return ImVec4(1.0f, 0.05f, 0.05f, 1.0f);
    case 9: return ImVec4(0.56f, 0.88f, 0.31f, 1.0f);
    case 15: return ImVec4(1.0f, 0.50f, 0.0f, 1.0f);
    case 16: return ImVec4(1.0f, 1.0f, 0.19f, 1.0f);
    case 17: return ImVec4(0.12f, 0.94f, 0.12f, 1.0f);
    case 35: return ImVec4(0.65f, 0.16f, 0.16f, 1.0f);
    default: return ImVec4(0.35f, 0.55f, 0.75f, 1.0f);
    }
}

std::string frontier_label(int mo_index, int homo) {
    if (homo < 0) {
        return "MO";
    }
    if (mo_index == homo) {
        return "HOMO";
    }
    if (mo_index == homo + 1) {
        return "LUMO";
    }
    if (mo_index < homo) {
        return "HOMO-" + std::to_string(homo - mo_index);
    }
    return "LUMO+" + std::to_string(mo_index - homo - 1);
}

}  // namespace

void draw_orbital_composition_panel(AppState& state,
                                    const sbox::backend::JobResult& result,
                                    const sbox::chem::MolecularSystem& mol) {
    if (!result.has_mo_data || result.mo_data.energies.size() == 0) {
        return;
    }

    if (!ImGui::Begin("Orbital Composition")) {
        ImGui::End();
        return;
    }

    if (ImGui::Button("<- Back to Dashboard")) {
        state.property_view = PropertyView::Dashboard;
    }
    ImGui::Separator();

    const int total_mo = result.mo_data.energies.size();
    const int homo = result.homo_index();
    const int lumo = result.lumo_index();
    if (state.selected_mo < 0) {
        state.selected_mo = homo >= 0 ? homo : 0;
    }
    state.selected_mo = std::clamp(state.selected_mo, 0, total_mo - 1);

    static int cached_mo = -1;
    static std::size_t cached_atom_count = 0;
    static sbox::analysis::OrbitalComposition cached;
    if (cached_mo != state.selected_mo || cached_atom_count != static_cast<std::size_t>(mol.num_atoms())) {
        cached = sbox::analysis::analyze_orbital_composition(result.mo_data, mol, state.selected_mo);
        cached_mo = state.selected_mo;
        cached_atom_count = static_cast<std::size_t>(mol.num_atoms());
    }

    ImGui::Text("MO #%d (%s)", state.selected_mo, frontier_label(state.selected_mo, homo).c_str());
    ImGui::Text("Energy: %.3f eV, Occupation: %.1f", cached.energy_eV, cached.occupation);
    ImGui::Separator();

    const float bar_width = ImGui::GetContentRegionAvail().x;
    const float bar_height = 24.0f;
    const ImVec2 bar_start = ImGui::GetCursorScreenPos();
    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    float x_offset = 0.0f;
    for (const auto& contrib : cached.atom_contributions) {
        const float segment_width = bar_width * static_cast<float>(contrib.total_weight);
        if (segment_width < 1.0f) {
            continue;
        }
        const ImVec4 color = cpk_color_imvec4(contrib.Z);
        draw_list->AddRectFilled(ImVec2(bar_start.x + x_offset, bar_start.y),
                                 ImVec2(bar_start.x + x_offset + segment_width, bar_start.y + bar_height),
                                 ImGui::GetColorU32(color),
                                 3.0f);
        if (segment_width > 40.0f) {
            char label[32];
            std::snprintf(label, sizeof(label), "%s %.0f%%", contrib.element.c_str(), contrib.total_weight * 100.0);
            draw_list->AddText(ImVec2(bar_start.x + x_offset + 4.0f, bar_start.y + 4.0f),
                               IM_COL32(255, 255, 255, 220),
                               label);
        }
        x_offset += segment_width;
    }
    draw_list->AddRect(bar_start,
                       ImVec2(bar_start.x + bar_width, bar_start.y + bar_height),
                       IM_COL32(220, 224, 232, 120),
                       3.0f);
    ImGui::Dummy(ImVec2(bar_width, bar_height + 4.0f));

    if (ImGui::BeginTable("OrbitalCompositionTable", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("Atom", ImGuiTableColumnFlags_WidthFixed, 48.0f);
        ImGui::TableSetupColumn("Element", ImGuiTableColumnFlags_WidthFixed, 64.0f);
        ImGui::TableSetupColumn("Weight", ImGuiTableColumnFlags_WidthFixed, 72.0f);
        ImGui::TableSetupColumn("AO Breakdown", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();

        const int shown = std::min<int>(8, static_cast<int>(cached.atom_contributions.size()));
        for (int i = 0; i < shown; ++i) {
            const auto& contrib = cached.atom_contributions[static_cast<std::size_t>(i)];
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%d", contrib.atom_index + 1);
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%s", contrib.element.c_str());
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%.1f%%", contrib.total_weight * 100.0);
            ImGui::TableSetColumnIndex(3);
            std::string breakdown;
            for (std::size_t j = 0; j < contrib.ao_contributions.size(); ++j) {
                if (j > 0) {
                    breakdown += "  ";
                }
                breakdown += contrib.ao_contributions[j].first + "(" +
                             std::to_string(static_cast<int>(std::round(contrib.ao_contributions[j].second * 100.0))) + "%)";
            }
            ImGui::TextWrapped("%s", breakdown.c_str());
        }
        ImGui::EndTable();
    }

    ImGui::Separator();
    ImGui::TextWrapped("%s", cached.summary.c_str());

    ImGui::Separator();
    if (ImGui::Button("<- Prev MO") && state.selected_mo > 0) {
        --state.selected_mo;
    }
    ImGui::SameLine();
    if (ImGui::Button("-> Next MO") && state.selected_mo + 1 < total_mo) {
        ++state.selected_mo;
    }
    ImGui::SameLine();
    if (ImGui::Button("HOMO") && homo >= 0) {
        state.selected_mo = homo;
    }
    ImGui::SameLine();
    if (ImGui::Button("LUMO") && lumo >= 0) {
        state.selected_mo = lumo;
    }

    ImGui::End();
}

}  // namespace sbox::ui
