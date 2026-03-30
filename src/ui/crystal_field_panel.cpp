#include "ui/crystal_field_panel.h"

#include "analysis/crystal_field.h"
#include "chem/coordination.h"
#include "core/elements.h"

#include <imgui.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace sbox::ui {

namespace {

constexpr double kEvToCm1 = 8065.54;

bool is_transition_metal(int Z) {
    return (Z >= 21 && Z <= 30) || (Z >= 39 && Z <= 48) || (Z >= 57 && Z <= 80) || (Z >= 89 && Z <= 112);
}

sbox::chem::CoordinationGeometry auto_geometry(const sbox::chem::MolecularSystem& mol, int metal_index) {
    const int cn = mol.coordination_number(metal_index);
    switch (cn) {
    case 4: return sbox::chem::CoordinationGeometry::Octahedral;  // fallback later overridden by UI if needed
    case 5: return sbox::chem::CoordinationGeometry::SquarePyramidal;
    case 6: return sbox::chem::CoordinationGeometry::Octahedral;
    default: return sbox::chem::CoordinationGeometry::Octahedral;
    }
}

const char* geometry_name(int mode) {
    switch (mode) {
    case 0: return "Auto-detect";
    case 1: return "Octahedral";
    case 2: return "Tetrahedral";
    case 3: return "Square Planar";
    default: return "Auto-detect";
    }
}

sbox::chem::CoordinationGeometry geometry_from_mode(int mode, const sbox::chem::MolecularSystem& mol, int metal_index) {
    switch (mode) {
    case 1: return sbox::chem::CoordinationGeometry::Octahedral;
    case 2: return sbox::chem::CoordinationGeometry::Tetrahedral;
    case 3: return sbox::chem::CoordinationGeometry::SquarePlanar;
    default: return auto_geometry(mol, metal_index);
    }
}

std::string orbital_label_for_index(int index, int homo) {
    if (index == homo) return "HOMO";
    if (index == homo + 1) return "LUMO";
    if (index < homo) return "HOMO-" + std::to_string(homo - index);
    return "LUMO+" + std::to_string(index - homo - 1);
}

std::string ligand_name_for_atom(const sbox::chem::MolecularSystem& mol, int atom_index) {
    return sbox::elements::get_element(mol.atom(atom_index).Z).symbol;
}

struct LevelInfo {
    std::string label;
    double energy = 0.0;
};

int count_unpaired(const std::vector<int>& occupancy) {
    int n = 0;
    for (int occ : occupancy) {
        if (occ == 1) {
            ++n;
        }
    }
    return n;
}

std::vector<int> fill_electrons(const std::vector<double>& energies, int electrons, bool high_spin) {
    std::vector<int> occ(energies.size(), 0);
    std::vector<int> order(energies.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return energies[a] < energies[b]; });

    if (high_spin) {
        for (int idx : order) {
            if (electrons <= 0) break;
            occ[idx] += 1;
            --electrons;
        }
    }
    for (int idx : order) {
        while (electrons > 0 && occ[idx] < 2) {
            occ[idx] += 1;
            --electrons;
        }
    }
    return occ;
}

}  // namespace

void draw_crystal_field_panel(AppState& state,
                              const sbox::backend::JobResult& result,
                              const sbox::chem::MolecularSystem& mol) {
    (void)state;
    if (!result.has_mo_data || mol.num_atoms() == 0) {
        return;
    }

    std::vector<int> metals;
    for (int i = 0; i < mol.num_atoms(); ++i) {
        if (is_transition_metal(mol.atom(i).Z)) {
            metals.push_back(i);
        }
    }
    if (metals.empty()) {
        return;
    }

    if (!ImGui::Begin("Crystal Field Diagram")) {
        ImGui::End();
        return;
    }

    static int selected_metal_slot = 0;
    static int geometry_mode = 0;
    static bool high_spin = true;
    static std::map<std::string, double> delta_cache;

    selected_metal_slot = std::clamp(selected_metal_slot, 0, static_cast<int>(metals.size()) - 1);
    if (metals.size() > 1) {
        std::vector<const char*> labels;
        std::vector<std::string> storage;
        for (int atom_index : metals) {
            storage.push_back(std::string(sbox::elements::get_element(mol.atom(atom_index).Z).symbol) + std::to_string(atom_index + 1));
        }
        labels.reserve(storage.size());
        for (const std::string& s : storage) labels.push_back(s.c_str());
        ImGui::Combo("Metal", &selected_metal_slot, labels.data(), static_cast<int>(labels.size()));
    }
    const int metal_index = metals[static_cast<std::size_t>(selected_metal_slot)];
    const int metal_Z = mol.atom(metal_index).Z;
    const int oxidation_state = mol.atom(metal_index).formal_charge != 0 ? mol.atom(metal_index).formal_charge : std::max(0, result.optimized_geometry.charge() == 0 ? 2 : result.optimized_geometry.charge());
    const int d_count = sbox::analysis::d_electron_count(metal_Z, oxidation_state);

    ImGui::Text("Metal: %s%d+ (d%d)",
                sbox::elements::get_element(metal_Z).symbol,
                oxidation_state,
                d_count);

    const char* modes[] = {"Auto-detect", "Octahedral", "Tetrahedral", "Square Planar"};
    ImGui::Combo("Geometry", &geometry_mode, modes, IM_ARRAYSIZE(modes));
    bool low_spin = !high_spin;
    if (ImGui::Checkbox("High Spin", &high_spin)) {
        low_spin = !high_spin;
    }
    ImGui::SameLine();
    if (ImGui::Checkbox("Low Spin", &low_spin)) {
        high_spin = !low_spin;
    }

    sbox::analysis::DOrbitalEnergies d_orbs = sbox::analysis::extract_d_orbitals(result.mo_data, mol, metal_index);
    const sbox::chem::CoordinationGeometry geom = geometry_from_mode(geometry_mode, mol, metal_index);
    sbox::analysis::identify_splitting(d_orbs, geom);

    const std::vector<LevelInfo> levels = {
        {"dxy", d_orbs.dxy},
        {"dxz", d_orbs.dxz},
        {"dyz", d_orbs.dyz},
        {"dz2", d_orbs.dz2},
        {"dx2-y2", d_orbs.dx2y2},
    };

    std::vector<double> level_energies;
    level_energies.reserve(levels.size());
    for (const auto& level : levels) {
        level_energies.push_back(level.energy);
    }
    const std::vector<int> occupancy = fill_electrons(level_energies, d_count, high_spin);
    const int unpaired = count_unpaired(occupancy);

    const float plot_width = ImGui::GetContentRegionAvail().x;
    const float plot_height = 260.0f;
    const ImVec2 origin = ImGui::GetCursorScreenPos();
    ImDrawList* draw = ImGui::GetWindowDrawList();
    draw->AddRect(origin, ImVec2(origin.x + plot_width, origin.y + plot_height), IM_COL32(90, 96, 110, 110), 4.0f);

    const double e_min = std::min({d_orbs.dxy, d_orbs.dxz, d_orbs.dyz, d_orbs.dz2, d_orbs.dx2y2}) - 0.3;
    const double e_max = std::max({d_orbs.dxy, d_orbs.dxz, d_orbs.dyz, d_orbs.dz2, d_orbs.dx2y2}) + 0.3;
    auto y_from_e = [&](double e) {
        const double t = (e - e_min) / std::max(1.0e-6, e_max - e_min);
        return static_cast<float>(origin.y + plot_height - 28.0f - t * (plot_height - 56.0f));
    };

    for (std::size_t i = 0; i < levels.size(); ++i) {
        const float y = y_from_e(levels[i].energy);
        const float x0 = origin.x + 56.0f;
        const float x1 = origin.x + plot_width - 70.0f;
        const ImU32 line_col = occupancy[i] > 0 ? IM_COL32(70, 150, 245, 255) : IM_COL32(160, 165, 176, 220);
        draw->AddLine(ImVec2(x0, y), ImVec2(x1, y), line_col, 2.5f);
        draw->AddText(ImVec2(x1 + 8.0f, y - 8.0f), IM_COL32(230, 234, 242, 220), levels[i].label.c_str());
        for (int e = 0; e < occupancy[i]; ++e) {
            const float arrow_x = x0 + 20.0f + 14.0f * static_cast<float>(e);
            draw->AddText(ImVec2(arrow_x, y - 22.0f), e == 0 ? IM_COL32(240, 200, 80, 255) : IM_COL32(230, 120, 180, 255), e == 0 ? "^" : "v");
        }
    }

    if (!d_orbs.groups.empty() && d_orbs.groups.size() >= 2) {
        const double lower = d_orbs.groups.front().average_energy;
        const double upper = d_orbs.groups.back().average_energy;
        const float x = origin.x + 24.0f;
        const float y0 = y_from_e(lower);
        const float y1 = y_from_e(upper);
        draw->AddLine(ImVec2(x, y0), ImVec2(x, y1), IM_COL32(220, 220, 130, 220), 2.0f);
        draw->AddText(ImVec2(x + 8.0f, 0.5f * (y0 + y1) - 8.0f),
                      IM_COL32(240, 230, 120, 230),
                      ("Delta = " + std::to_string(std::round(std::abs(upper - lower) * 100.0) / 100.0) + " eV").c_str());
    }
    ImGui::Dummy(ImVec2(plot_width, plot_height + 8.0f));

    const double delta_ev = geom == sbox::chem::CoordinationGeometry::Tetrahedral ? d_orbs.delta_tet() : d_orbs.delta_oct();
    const double wavelength_nm = std::abs(delta_ev) > 1.0e-6 ? 1240.0 / std::abs(delta_ev) : 0.0;
    const double mu_b = std::sqrt(static_cast<double>(unpaired) * (static_cast<double>(unpaired) + 2.0));
    const double cfse_dq = geom == sbox::chem::CoordinationGeometry::Octahedral
                               ? sbox::analysis::octahedral_cfse_dq(d_count, high_spin)
                               : 0.0;

    ImGui::Text("Crystal field splitting: Delta = %.2f eV (%.0f cm^-1)", delta_ev, delta_ev * kEvToCm1);
    if (wavelength_nm > 0.0) {
        ImGui::Text("Predicted absorption: ~%.0f nm", wavelength_nm);
    }
    ImGui::Text("Spin state: %s", high_spin ? "high-spin" : "low-spin");
    ImGui::Text("Unpaired electrons: %d", unpaired);
    ImGui::Text("Expected magnetic moment: %.2f muB", mu_b);
    if (geom == sbox::chem::CoordinationGeometry::Octahedral) {
        ImGui::Text("CFSE: %.2f Dq", cfse_dq);
    }
    if (d_orbs.ordering_warning) {
        ImGui::TextColored(ImVec4(0.95f, 0.45f, 0.30f, 1.0f), "Computed ordering disagrees with the idealized crystal field pattern.");
    }

    static const std::vector<std::string> series = {"I-", "Br-", "Cl-", "F-", "OH-", "H2O", "NH3", "en", "NO2-", "CN-", "CO"};
    std::vector<std::string> present;
    for (int n : mol.neighbors(metal_index)) {
        present.push_back(ligand_name_for_atom(mol, n));
    }
    std::string series_text = "Weak <- ";
    for (std::size_t i = 0; i < series.size(); ++i) {
        bool highlight = std::find(present.begin(), present.end(), series[i]) != present.end();
        if (highlight) {
            series_text += "[" + series[i] + "]";
        } else {
            series_text += series[i];
        }
        if (i + 1 < series.size()) {
            series_text += " < ";
        }
    }
    series_text += " -> Strong";
    ImGui::Separator();
    ImGui::TextWrapped("Ligand field strength: %s", series_text.c_str());

    const std::string cache_key = std::string(sbox::elements::get_element(metal_Z).symbol) + ":" + (present.empty() ? "unknown" : present.front());
    delta_cache[cache_key] = delta_ev;
    if (!delta_cache.empty()) {
        ImGui::TextUnformatted("Comparison cache:");
        for (const auto& [key, value] : delta_cache) {
            ImGui::BulletText("%s: Delta = %.2f eV", key.c_str(), value);
        }
    }

    ImGui::End();
}

}  // namespace sbox::ui
