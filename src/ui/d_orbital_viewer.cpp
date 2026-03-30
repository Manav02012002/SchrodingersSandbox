#include "ui/d_orbital_viewer.h"

#include "analysis/orbital_composition.h"
#include "core/elements.h"

#include <imgui.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

namespace sbox::ui {

namespace {

struct DOrbitalSlot {
    const char* key = "";
    const char* label = "";
    double energy = 0.0;
    int mo_index = -1;
};

bool is_transition_metal(int Z) {
    return (Z >= 21 && Z <= 30) || (Z >= 39 && Z <= 48) || (Z >= 57 && Z <= 80) || (Z >= 89 && Z <= 112);
}

sbox::chem::CoordinationGeometry inferred_geometry(const sbox::chem::MolecularSystem& mol, int metal_index) {
    const int cn = mol.coordination_number(metal_index);
    switch (cn) {
    case 4:
        return sbox::chem::CoordinationGeometry::SquarePlanar;
    case 6:
        return sbox::chem::CoordinationGeometry::Octahedral;
    default:
        return sbox::chem::CoordinationGeometry::Octahedral;
    }
}

std::vector<int> fill_electrons(const std::vector<double>& energies, int electrons, bool high_spin) {
    std::vector<int> occ(energies.size(), 0);
    std::vector<int> order(energies.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return energies[a] < energies[b]; });

    if (high_spin) {
        for (int idx : order) {
            if (electrons <= 0) {
                break;
            }
            occ[idx] += 1;
            --electrons;
        }
    }
    for (int idx : order) {
        while (electrons > 0 && occ[idx] < 2) {
            ++occ[idx];
            --electrons;
        }
    }
    return occ;
}

int count_unpaired(const std::vector<int>& occupancy) {
    int total = 0;
    for (int occ : occupancy) {
        if (occ == 1) {
            ++total;
        }
    }
    return total;
}

std::string electron_configuration_text(const sbox::analysis::DOrbitalEnergies& d_orbs,
                                        const std::vector<DOrbitalSlot>& orbitals,
                                        const std::vector<int>& occupancy) {
    if (!d_orbs.groups.empty()) {
        std::string text;
        for (std::size_t i = 0; i < d_orbs.groups.size(); ++i) {
            const auto& group = d_orbs.groups[i];
            int electrons = 0;
            for (const std::string& label : group.orbitals) {
                for (std::size_t j = 0; j < orbitals.size(); ++j) {
                    if (label == orbitals[j].key) {
                        electrons += occupancy[j];
                    }
                }
            }
            if (!text.empty()) {
                text += " ";
            }
            text += "(" + group.label + ")^" + std::to_string(electrons);
        }
        return text;
    }

    std::string text;
    for (std::size_t i = 0; i < orbitals.size(); ++i) {
        if (!text.empty()) {
            text += " ";
        }
        text += orbitals[i].label;
        text += "^";
        text += std::to_string(occupancy[i]);
    }
    return text;
}

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
    case 26: return ImVec4(0.88f, 0.40f, 0.20f, 1.0f);
    case 27: return ImVec4(0.94f, 0.56f, 0.62f, 1.0f);
    case 28: return ImVec4(0.31f, 0.82f, 0.31f, 1.0f);
    case 29: return ImVec4(0.78f, 0.50f, 0.20f, 1.0f);
    default: return ImVec4(0.35f, 0.55f, 0.75f, 1.0f);
    }
}

void draw_lobe(ImDrawList* draw, const ImVec2& center, const ImVec2& offset, ImU32 color) {
    draw->AddCircleFilled(ImVec2(center.x + offset.x, center.y + offset.y), 6.0f, color, 16);
}

void draw_d_shape(ImDrawList* draw, const ImVec2& min, const ImVec2& size, const char* label) {
    const ImVec2 center(min.x + size.x * 0.5f, min.y + size.y * 0.5f);
    const ImU32 pos = IM_COL32(70, 140, 245, 220);
    const ImU32 neg = IM_COL32(220, 80, 90, 220);
    const ImU32 frame = IM_COL32(130, 138, 152, 90);
    draw->AddRect(min, ImVec2(min.x + size.x, min.y + size.y), frame, 4.0f);

    if (std::string(label) == "dxy") {
        draw_lobe(draw, center, ImVec2(-10.0f, -10.0f), pos);
        draw_lobe(draw, center, ImVec2(10.0f, 10.0f), pos);
        draw_lobe(draw, center, ImVec2(10.0f, -10.0f), neg);
        draw_lobe(draw, center, ImVec2(-10.0f, 10.0f), neg);
    } else if (std::string(label) == "dxz") {
        draw_lobe(draw, center, ImVec2(-12.0f, 0.0f), pos);
        draw_lobe(draw, center, ImVec2(12.0f, 0.0f), pos);
        draw_lobe(draw, center, ImVec2(0.0f, -12.0f), neg);
        draw_lobe(draw, center, ImVec2(0.0f, 12.0f), neg);
        draw->AddLine(ImVec2(center.x - 10.0f, center.y), ImVec2(center.x + 10.0f, center.y), IM_COL32(200, 200, 200, 90), 1.0f);
    } else if (std::string(label) == "dyz") {
        draw_lobe(draw, center, ImVec2(-12.0f, 0.0f), neg);
        draw_lobe(draw, center, ImVec2(12.0f, 0.0f), neg);
        draw_lobe(draw, center, ImVec2(0.0f, -12.0f), pos);
        draw_lobe(draw, center, ImVec2(0.0f, 12.0f), pos);
        draw->AddLine(ImVec2(center.x, center.y - 10.0f), ImVec2(center.x, center.y + 10.0f), IM_COL32(200, 200, 200, 90), 1.0f);
    } else if (std::string(label) == "dz2") {
        draw->AddCircleFilled(ImVec2(center.x, center.y - 10.0f), 7.0f, pos, 16);
        draw->AddCircleFilled(ImVec2(center.x, center.y + 10.0f), 7.0f, neg, 16);
        draw->AddCircle(center, 9.0f, IM_COL32(230, 220, 140, 190), 20, 2.0f);
    } else {
        draw_lobe(draw, center, ImVec2(-12.0f, 0.0f), pos);
        draw_lobe(draw, center, ImVec2(12.0f, 0.0f), pos);
        draw_lobe(draw, center, ImVec2(0.0f, -12.0f), neg);
        draw_lobe(draw, center, ImVec2(0.0f, 12.0f), neg);
    }
}

std::vector<int> transition_metals(const sbox::chem::MolecularSystem& mol) {
    std::vector<int> metals;
    for (int i = 0; i < mol.num_atoms(); ++i) {
        if (is_transition_metal(mol.atom(i).Z)) {
            metals.push_back(i);
        }
    }
    return metals;
}

std::vector<DOrbitalSlot> slots_from_d_orbs(const sbox::analysis::DOrbitalEnergies& d_orbs) {
    return {
        {"dxy", "dxy", d_orbs.dxy, d_orbs.mo_dxy},
        {"dxz", "dxz", d_orbs.dxz, d_orbs.mo_dxz},
        {"dyz", "dyz", d_orbs.dyz, d_orbs.mo_dyz},
        {"dz2", "dz2", d_orbs.dz2, d_orbs.mo_dz2},
        {"dx2y2", "dx2-y2", d_orbs.dx2y2, d_orbs.mo_dx2y2},
    };
}

std::string free_ion_note(const sbox::analysis::DOrbitalEnergies& d_orbs) {
    char buffer[128];
    std::snprintf(buffer, sizeof(buffer), "Free ion reference: all five d levels degenerate at %.2f eV.", d_orbs.mean_energy());
    return buffer;
}

}  // namespace

void draw_d_orbital_viewer(AppState& state,
                           const sbox::backend::JobResult& result,
                           const sbox::chem::MolecularSystem& mol,
                           const sbox::analysis::DOrbitalEnergies& d_orbs) {
    if (!result.has_mo_data || mol.num_atoms() == 0) {
        return;
    }

    const std::vector<int> metals = transition_metals(mol);
    if (metals.empty()) {
        return;
    }

    if (!ImGui::Begin("d-Orbital Viewer")) {
        ImGui::End();
        return;
    }

    state.selected_crystal_field_metal = std::clamp(state.selected_crystal_field_metal, 0, static_cast<int>(metals.size()) - 1);
    const int metal_index = metals[static_cast<std::size_t>(state.selected_crystal_field_metal)];

    sbox::analysis::DOrbitalEnergies active_d_orbs = d_orbs;
    if (state.selected_crystal_field_metal != 0 || active_d_orbs.mo_dxy < 0 || active_d_orbs.mo_dxz < 0 ||
        active_d_orbs.mo_dyz < 0 || active_d_orbs.mo_dz2 < 0 || active_d_orbs.mo_dx2y2 < 0) {
        active_d_orbs = sbox::analysis::extract_d_orbitals(result.mo_data, mol, metal_index);
        sbox::analysis::identify_splitting(active_d_orbs, inferred_geometry(mol, metal_index));
    }

    const int oxidation_state = mol.atom(metal_index).formal_charge != 0 ? mol.atom(metal_index).formal_charge : 2;
    const int d_electrons = sbox::analysis::d_electron_count(mol.atom(metal_index).Z, oxidation_state);
    std::vector<DOrbitalSlot> orbitals = slots_from_d_orbs(active_d_orbs);
    state.selected_d_orbital = std::clamp(state.selected_d_orbital, 0, static_cast<int>(orbitals.size()) - 1);

    std::vector<double> level_energies;
    level_energies.reserve(orbitals.size());
    for (const DOrbitalSlot& orbital : orbitals) {
        level_energies.push_back(orbital.energy);
    }
    const std::vector<int> occupancy = fill_electrons(level_energies, d_electrons, state.d_orbital_high_spin);
    const int unpaired = count_unpaired(occupancy);

    ImGui::Text("Metal center: %s%d+ (d%d)",
                sbox::elements::get_element(mol.atom(metal_index).Z).symbol,
                oxidation_state,
                d_electrons);
    ImGui::Separator();

    ImGui::TextUnformatted("Orbital Selector");
    const float available = ImGui::GetContentRegionAvail().x;
    const float button_width = (available - 4.0f * ImGui::GetStyle().ItemSpacing.x) / 5.0f;
    for (std::size_t i = 0; i < orbitals.size(); ++i) {
        ImGui::PushID(static_cast<int>(i));
        const bool selected = static_cast<int>(i) == state.selected_d_orbital;
        const bool occupied = occupancy[i] > 0;
        const ImVec4 fill = selected ? ImVec4(0.10f, 0.56f, 0.66f, 1.0f)
                                     : (occupied ? ImVec4(0.16f, 0.36f, 0.78f, 1.0f) : ImVec4(0.22f, 0.24f, 0.28f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_Button, fill);
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(fill.x + 0.08f, fill.y + 0.08f, fill.z + 0.08f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, fill);
        char button_label[64];
        std::snprintf(button_label, sizeof(button_label), "%s\n%.2f eV", orbitals[i].label, orbitals[i].energy);
        if (ImGui::Button(button_label, ImVec2(button_width, 42.0f))) {
            state.selected_d_orbital = static_cast<int>(i);
            if (orbitals[i].mo_index >= 0) {
                state.selected_mo = orbitals[i].mo_index;
            }
        }
        ImGui::PopStyleColor(3);
        if (i + 1 < orbitals.size()) {
            ImGui::SameLine();
        }
        ImGui::PopID();
    }

    const float icon_y = ImGui::GetCursorScreenPos().y + 6.0f;
    ImDrawList* draw = ImGui::GetWindowDrawList();
    float icon_x = ImGui::GetCursorScreenPos().x;
    for (const auto& orbital : orbitals) {
        draw_d_shape(draw, ImVec2(icon_x, icon_y), ImVec2(button_width, 44.0f), orbital.key);
        icon_x += button_width + ImGui::GetStyle().ItemSpacing.x;
    }
    ImGui::Dummy(ImVec2(available, 56.0f));

    ImGui::Separator();
    ImGui::TextUnformatted("Occupancy Display");

    const float diagram_width = ImGui::GetContentRegionAvail().x;
    const float diagram_height = 180.0f;
    const ImVec2 origin = ImGui::GetCursorScreenPos();
    draw->AddRect(origin, ImVec2(origin.x + diagram_width, origin.y + diagram_height), IM_COL32(92, 100, 112, 100), 4.0f);

    const double e_min = *std::min_element(level_energies.begin(), level_energies.end()) - 0.3;
    const double e_max = *std::max_element(level_energies.begin(), level_energies.end()) + 0.3;
    auto y_from_e = [&](double energy) {
        const double t = (energy - e_min) / std::max(1.0e-6, e_max - e_min);
        return static_cast<float>(origin.y + diagram_height - 24.0f - t * (diagram_height - 48.0f));
    };

    for (std::size_t i = 0; i < orbitals.size(); ++i) {
        const float y = y_from_e(orbitals[i].energy);
        const float x0 = origin.x + 72.0f;
        const float x1 = origin.x + diagram_width - 88.0f;
        const ImU32 col = occupancy[i] > 0 ? IM_COL32(80, 150, 245, 255) : IM_COL32(160, 166, 176, 220);
        draw->AddLine(ImVec2(x0, y), ImVec2(x1, y), col, static_cast<int>(i) == state.selected_d_orbital ? 3.5f : 2.0f);
        draw->AddText(ImVec2(origin.x + 12.0f, y - 8.0f), IM_COL32(230, 234, 242, 220), orbitals[i].label);
        for (int e = 0; e < occupancy[i]; ++e) {
            const char* arrow = e == 0 ? "^" : "v";
            const ImU32 arrow_col = e == 0 ? IM_COL32(248, 210, 96, 255) : IM_COL32(235, 135, 182, 255);
            draw->AddText(ImVec2(x0 + 24.0f + 14.0f * static_cast<float>(e), y - 22.0f), arrow_col, arrow);
        }
    }
    if (!active_d_orbs.groups.empty() && active_d_orbs.groups.size() >= 2) {
        const float x = origin.x + diagram_width - 42.0f;
        const float y0 = y_from_e(active_d_orbs.groups.front().average_energy);
        const float y1 = y_from_e(active_d_orbs.groups.back().average_energy);
        draw->AddLine(ImVec2(x, y0), ImVec2(x, y1), IM_COL32(228, 220, 130, 220), 2.0f);
        char gap[64];
        std::snprintf(gap, sizeof(gap), "D = %.2f eV", std::abs(y1 - y0) > 0.0f ? std::abs(active_d_orbs.groups.back().average_energy - active_d_orbs.groups.front().average_energy) : 0.0);
        draw->AddText(ImVec2(x - 28.0f, 0.5f * (y0 + y1) - 8.0f), IM_COL32(238, 226, 150, 230), gap);
    }
    if (state.d_orbital_compare_free_ion) {
        const float y = y_from_e(active_d_orbs.mean_energy());
        draw->AddLine(ImVec2(origin.x + 60.0f, y), ImVec2(origin.x + diagram_width - 10.0f, y), IM_COL32(180, 180, 180, 110), 1.0f);
    }
    ImGui::Dummy(ImVec2(diagram_width, diagram_height + 4.0f));

    ImGui::Text("Configuration: %s", electron_configuration_text(active_d_orbs, orbitals, occupancy).c_str());
    ImGui::Text("Unpaired electrons: %d", unpaired);
    ImGui::Text("Spin state: %s", state.d_orbital_high_spin ? "high-spin" : "low-spin");

    ImGui::Separator();
    ImGui::TextUnformatted("d-Orbital Mixing");

    static int cached_metal_index = -1;
    static std::array<int, 5> cached_mo_indices = {{-2, -2, -2, -2, -2}};
    static std::array<sbox::analysis::OrbitalComposition, 5> cached_compositions;

    auto refresh_compositions = [&]() {
        for (std::size_t i = 0; i < orbitals.size(); ++i) {
            cached_mo_indices[i] = orbitals[i].mo_index;
            if (orbitals[i].mo_index >= 0) {
                cached_compositions[i] = sbox::analysis::analyze_orbital_composition(result.mo_data, mol, orbitals[i].mo_index);
            } else {
                cached_compositions[i] = {};
            }
        }
        cached_metal_index = metal_index;
    };

    bool needs_refresh = cached_metal_index != metal_index;
    for (std::size_t i = 0; i < orbitals.size(); ++i) {
        if (cached_mo_indices[i] != orbitals[i].mo_index) {
            needs_refresh = true;
        }
    }
    if (needs_refresh) {
        refresh_compositions();
    }

    for (std::size_t i = 0; i < orbitals.size(); ++i) {
        if (orbitals[i].mo_index < 0) {
            ImGui::TextDisabled("%s: no mapped MO", orbitals[i].label);
            continue;
        }
        const auto& comp = cached_compositions[i];
        ImGui::PushID(static_cast<int>(i) + 100);
        ImGui::Text("%s MO (#%d): %s", orbitals[i].label, orbitals[i].mo_index, comp.summary.c_str());

        const float bar_width = ImGui::GetContentRegionAvail().x;
        const float bar_height = 20.0f;
        const ImVec2 bar_start = ImGui::GetCursorScreenPos();
        float x_offset = 0.0f;
        for (const auto& contrib : comp.atom_contributions) {
            const float w = bar_width * static_cast<float>(contrib.total_weight);
            if (w < 1.0f) {
                continue;
            }
            const ImVec4 color = cpk_color_imvec4(contrib.Z);
            draw->AddRectFilled(ImVec2(bar_start.x + x_offset, bar_start.y),
                                ImVec2(bar_start.x + x_offset + w, bar_start.y + bar_height),
                                ImGui::GetColorU32(color),
                                2.0f);
            if (w > 44.0f) {
                char label[32];
                std::snprintf(label, sizeof(label), "%s %.0f%%", contrib.element.c_str(), contrib.total_weight * 100.0);
                draw->AddText(ImVec2(bar_start.x + x_offset + 4.0f, bar_start.y + 2.0f), IM_COL32(255, 255, 255, 220), label);
            }
            x_offset += w;
        }
        draw->AddRect(bar_start, ImVec2(bar_start.x + bar_width, bar_start.y + bar_height), IM_COL32(225, 228, 235, 120), 2.0f);
        ImGui::Dummy(ImVec2(bar_width, bar_height + 4.0f));
        ImGui::PopID();
    }

    ImGui::Separator();
    if (ImGui::Button("<- Previous Metal") && state.selected_crystal_field_metal > 0) {
        --state.selected_crystal_field_metal;
    }
    ImGui::SameLine();
    if (ImGui::Button("Next Metal ->") && state.selected_crystal_field_metal + 1 < static_cast<int>(metals.size())) {
        ++state.selected_crystal_field_metal;
    }
    ImGui::SameLine();
    if (ImGui::Button("Recalculate")) {
        cached_metal_index = -1;
    }
    ImGui::SameLine();
    ImGui::Checkbox("Compare with Free Ion", &state.d_orbital_compare_free_ion);

    bool low_spin = !state.d_orbital_high_spin;
    if (ImGui::Checkbox("High Spin", &state.d_orbital_high_spin)) {
        low_spin = !state.d_orbital_high_spin;
    }
    ImGui::SameLine();
    if (ImGui::Checkbox("Low Spin", &low_spin)) {
        state.d_orbital_high_spin = !low_spin;
    }

    if (state.d_orbital_compare_free_ion) {
        ImGui::TextWrapped("%s", free_ion_note(active_d_orbs).c_str());
    }

    ImGui::End();
}

}  // namespace sbox::ui
