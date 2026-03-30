#include "ui/complex_builder.h"

#include "analysis/crystal_field.h"
#include "chem/coordination.h"
#include "core/elements.h"
#include "editor/command.h"
#include "ui/app_state.h"

#include <imgui.h>

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace sbox::ui {

namespace {

using sbox::chem::CoordinationGeometry;
using sbox::chem::LigandDenticity;
using sbox::chem::LigandEntry;

constexpr double kPi = 3.14159265358979323846;

struct GeometryChoice {
    CoordinationGeometry geometry;
    const char* label;
};

const std::array<GeometryChoice, 11> kGeometryChoices = {{
    {CoordinationGeometry::Linear, "Linear (2)"},
    {CoordinationGeometry::TrigonalPlanar, "Trigonal Planar (3)"},
    {CoordinationGeometry::TShaped, "T-Shaped (3)"},
    {CoordinationGeometry::Tetrahedral, "Tetrahedral (4)"},
    {CoordinationGeometry::SquarePlanar, "Square Planar (4)"},
    {CoordinationGeometry::SeeSaw, "See-Saw (4)"},
    {CoordinationGeometry::TrigonalBipyramidal, "Trigonal Bipyramidal (5)"},
    {CoordinationGeometry::SquarePyramidal, "Square Pyramidal (5)"},
    {CoordinationGeometry::Octahedral, "Octahedral (6)"},
    {CoordinationGeometry::PentagonalBipyramidal, "Pentagonal Bipyramidal (7)"},
    {CoordinationGeometry::SquareAntiprismatic, "Square Antiprismatic (8)"},
}};

struct MetalChoice {
    int Z;
    std::vector<int> oxidation_states;
};

const std::array<MetalChoice, 30> kTransitionMetals = {{
    {21, {3}},
    {22, {2, 3, 4}},
    {23, {2, 3, 5}},
    {24, {2, 3, 6}},
    {25, {2, 4, 7}},
    {26, {2, 3}},
    {27, {2, 3}},
    {28, {2, 3}},
    {29, {1, 2}},
    {30, {2}},
    {39, {3}},
    {40, {4}},
    {41, {3, 5}},
    {42, {4, 6}},
    {43, {4, 7}},
    {44, {2, 3}},
    {45, {1, 3}},
    {46, {2, 4}},
    {47, {1}},
    {48, {2}},
    {72, {4}},
    {73, {5}},
    {74, {4, 6}},
    {75, {4, 7}},
    {76, {4}},
    {77, {3, 4}},
    {78, {2, 4}},
    {79, {1, 3}},
    {80, {2}},
    {0, {}},
}};

int geometry_choice_index(CoordinationGeometry geometry) {
    for (std::size_t i = 0; i < kGeometryChoices.size(); ++i) {
        if (kGeometryChoices[i].geometry == geometry) {
            return static_cast<int>(i);
        }
    }
    return 0;
}

const std::vector<int>& common_oxidation_states(int Z) {
    static const std::vector<int> fallback = {2};
    for (const MetalChoice& metal : kTransitionMetals) {
        if (metal.Z == Z) {
            return metal.oxidation_states.empty() ? fallback : metal.oxidation_states;
        }
    }
    return fallback;
}

void ensure_site_count(AppState::ComplexBuilderState& builder) {
    const int count = sbox::chem::get_template(builder.geometry).coordination_number;
    if (builder.site_ligands.size() != static_cast<std::size_t>(count)) {
        builder.site_ligands.assign(static_cast<std::size_t>(count), "water");
    }
}

const LigandEntry* ligand_for_name(const sbox::chem::LigandLibrary& library, const std::string& name) {
    const LigandEntry* ligand = library.find(name);
    if (ligand != nullptr) {
        return ligand;
    }
    return library.find("water");
}

int denticity_value(LigandDenticity denticity) {
    return static_cast<int>(denticity);
}

void clear_overlapping_sites(std::vector<std::string>& site_ligands,
                             const sbox::chem::LigandLibrary& library,
                             int site_index) {
    for (int leader = 0; leader < static_cast<int>(site_ligands.size()); ++leader) {
        const LigandEntry* ligand = ligand_for_name(library, site_ligands[static_cast<std::size_t>(leader)]);
        if (ligand == nullptr) {
            continue;
        }
        const int width = denticity_value(ligand->denticity);
        const int end = leader + width - 1;
        if (site_index >= leader && site_index <= end) {
            for (int i = leader; i <= std::min(end, static_cast<int>(site_ligands.size()) - 1); ++i) {
                site_ligands[static_cast<std::size_t>(i)] = "water";
            }
        }
    }
}

bool site_is_occupied_by_previous(const AppState::ComplexBuilderState& builder,
                                  const sbox::chem::LigandLibrary& library,
                                  int site_index,
                                  int* leader_out = nullptr,
                                  const LigandEntry** ligand_out = nullptr) {
    for (int leader = 0; leader < site_index; ++leader) {
        const LigandEntry* ligand = ligand_for_name(library, builder.site_ligands[static_cast<std::size_t>(leader)]);
        if (ligand == nullptr) {
            continue;
        }
        const int width = denticity_value(ligand->denticity);
        if (site_index < leader + width) {
            if (leader_out != nullptr) {
                *leader_out = leader;
            }
            if (ligand_out != nullptr) {
                *ligand_out = ligand;
            }
            return true;
        }
    }
    return false;
}

void assign_ligand(AppState::ComplexBuilderState& builder,
                   const sbox::chem::LigandLibrary& library,
                   int site_index,
                   const LigandEntry& ligand) {
    clear_overlapping_sites(builder.site_ligands, library, site_index);
    const int width = denticity_value(ligand.denticity);
    for (int i = 0; i < width && site_index + i < static_cast<int>(builder.site_ligands.size()); ++i) {
        builder.site_ligands[static_cast<std::size_t>(site_index + i)] = ligand.name;
    }
}

void assign_uniform_ligand(AppState::ComplexBuilderState& builder,
                           const LigandEntry& ligand) {
    std::fill(builder.site_ligands.begin(), builder.site_ligands.end(), std::string("water"));
    const int width = denticity_value(ligand.denticity);
    for (int site = 0; site < static_cast<int>(builder.site_ligands.size()); site += width) {
        if (site + width > static_cast<int>(builder.site_ligands.size())) {
            break;
        }
        for (int i = 0; i < width; ++i) {
            builder.site_ligands[static_cast<std::size_t>(site + i)] = ligand.name;
        }
    }
}

bool selection_tiles_geometry(const AppState::ComplexBuilderState& builder,
                              const sbox::chem::LigandLibrary& library) {
    int site = 0;
    while (site < static_cast<int>(builder.site_ligands.size())) {
        if (site_is_occupied_by_previous(builder, library, site)) {
            ++site;
            continue;
        }
        const LigandEntry* ligand = ligand_for_name(library, builder.site_ligands[static_cast<std::size_t>(site)]);
        if (ligand == nullptr) {
            return false;
        }
        const int width = denticity_value(ligand->denticity);
        if (site + width > static_cast<int>(builder.site_ligands.size())) {
            return false;
        }
        for (int i = 0; i < width; ++i) {
            if (builder.site_ligands[static_cast<std::size_t>(site + i)] != ligand->name) {
                return false;
            }
        }
        site += width;
    }
    return true;
}

std::string oxidation_string(int oxidation_state) {
    if (oxidation_state > 0) {
        return "+" + std::to_string(oxidation_state);
    }
    return std::to_string(oxidation_state);
}

std::string superscript_charge(int charge) {
    if (charge == 0) {
        return "";
    }
    if (charge == 1) {
        return "+";
    }
    if (charge == -1) {
        return "-";
    }
    return charge > 0 ? std::to_string(charge) + "+" : std::to_string(-charge) + "-";
}

std::string summarize_complex(const AppState::ComplexBuilderState& builder,
                              const sbox::chem::LigandLibrary& library) {
    std::map<std::string, int> counts;
    int total_charge = builder.oxidation_state;
    for (int site = 0; site < static_cast<int>(builder.site_ligands.size()); ++site) {
        if (site_is_occupied_by_previous(builder, library, site)) {
            continue;
        }
        const LigandEntry* ligand = ligand_for_name(library, builder.site_ligands[static_cast<std::size_t>(site)]);
        if (ligand == nullptr) {
            continue;
        }
        counts[ligand->formula.empty() ? ligand->abbreviation : ligand->formula] += 1;
        total_charge += ligand->charge;
    }

    std::ostringstream out;
    out << "[" << sbox::elements::get_element(builder.metal_Z).symbol;
    for (const auto& [label, count] : counts) {
        out << "(" << label << ")";
        if (count > 1) {
            out << count;
        }
    }
    out << "]";
    out << superscript_charge(total_charge);
    return out.str();
}

double expected_delta(const AppState::ComplexBuilderState& builder,
                      const sbox::chem::LigandLibrary& library) {
    double sum = 0.0;
    int count = 0;
    for (int site = 0; site < static_cast<int>(builder.site_ligands.size()); ++site) {
        if (site_is_occupied_by_previous(builder, library, site)) {
            continue;
        }
        const LigandEntry* ligand = ligand_for_name(library, builder.site_ligands[static_cast<std::size_t>(site)]);
        if (ligand != nullptr && ligand->typical_delta_oct_eV > 0.0) {
            sum += ligand->typical_delta_oct_eV;
            ++count;
        }
    }
    return count > 0 ? sum / static_cast<double>(count) : 0.0;
}

std::string expected_spin_text(const AppState::ComplexBuilderState& builder,
                               const sbox::chem::LigandLibrary& library) {
    const double delta = expected_delta(builder, library);
    if (builder.geometry == CoordinationGeometry::Octahedral) {
        return delta >= 1.6 ? "low-spin octahedral" : "high-spin octahedral";
    }
    if (builder.geometry == CoordinationGeometry::SquarePlanar) {
        return "square-planar ligand field";
    }
    if (builder.geometry == CoordinationGeometry::Tetrahedral) {
        return "tetrahedral ligand field";
    }
    return "coordination complex";
}

void draw_metal_selector(AppState::ComplexBuilderState& builder) {
    const ImDrawList* unused = nullptr;
    (void)unused;
    static const std::array<std::array<int, 10>, 3> rows = {{
        {{21, 22, 23, 24, 25, 26, 27, 28, 29, 30}},
        {{39, 40, 41, 42, 43, 44, 45, 46, 47, 48}},
        {{72, 73, 74, 75, 76, 77, 78, 79, 80, 0}},
    }};

    for (std::size_t r = 0; r < rows.size(); ++r) {
        for (std::size_t c = 0; c < rows[r].size(); ++c) {
            const int Z = rows[r][c];
            if (Z == 0) {
                ImGui::Dummy(ImVec2(32.0f, 24.0f));
            } else {
                const bool active = builder.metal_Z == Z;
                if (active) {
                    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.14f, 0.46f, 0.56f, 1.0f));
                    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.18f, 0.55f, 0.66f, 1.0f));
                }
                if (ImGui::Button(sbox::elements::get_element(Z).symbol, ImVec2(32.0f, 24.0f))) {
                    builder.metal_Z = Z;
                    const auto& ox = common_oxidation_states(Z);
                    if (std::find(ox.begin(), ox.end(), builder.oxidation_state) == ox.end()) {
                        builder.oxidation_state = ox.front();
                    }
                }
                if (active) {
                    ImGui::PopStyleColor(2);
                }
            }
            if (c + 1 < rows[r].size()) {
                ImGui::SameLine();
            }
        }
    }
}

void draw_geometry_schematic(CoordinationGeometry geometry) {
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    const ImVec2 size(160.0f, 120.0f);
    ImDrawList* draw = ImGui::GetWindowDrawList();
    draw->AddRect(p0, ImVec2(p0.x + size.x, p0.y + size.y), IM_COL32(80, 88, 98, 180), 4.0f);
    const ImVec2 center(p0.x + size.x * 0.5f, p0.y + size.y * 0.5f);
    draw->AddCircleFilled(center, 6.0f, IM_COL32(240, 210, 110, 255));

    auto add_site = [&](float x, float y) {
        const ImVec2 pt(p0.x + x, p0.y + y);
        draw->AddLine(center, pt, IM_COL32(110, 120, 136, 180), 1.5f);
        draw->AddCircleFilled(pt, 5.0f, IM_COL32(86, 176, 196, 220));
    };

    switch (geometry) {
    case CoordinationGeometry::Linear:
        add_site(size.x * 0.5f, 18.0f);
        add_site(size.x * 0.5f, size.y - 18.0f);
        break;
    case CoordinationGeometry::TrigonalPlanar:
        add_site(size.x * 0.5f, 18.0f);
        add_site(30.0f, size.y - 24.0f);
        add_site(size.x - 30.0f, size.y - 24.0f);
        break;
    case CoordinationGeometry::TShaped:
        add_site(26.0f, size.y * 0.5f);
        add_site(size.x - 26.0f, size.y * 0.5f);
        add_site(size.x * 0.5f, 20.0f);
        break;
    case CoordinationGeometry::Tetrahedral:
        add_site(size.x * 0.5f, 18.0f);
        add_site(32.0f, size.y - 28.0f);
        add_site(size.x - 32.0f, size.y - 28.0f);
        add_site(size.x * 0.5f, size.y * 0.5f + 8.0f);
        break;
    case CoordinationGeometry::SquarePlanar:
        add_site(size.x * 0.5f, 18.0f);
        add_site(24.0f, size.y * 0.5f);
        add_site(size.x - 24.0f, size.y * 0.5f);
        add_site(size.x * 0.5f, size.y - 18.0f);
        break;
    case CoordinationGeometry::SeeSaw:
        add_site(size.x * 0.5f, 18.0f);
        add_site(size.x * 0.5f, size.y - 18.0f);
        add_site(26.0f, size.y * 0.5f + 10.0f);
        add_site(size.x - 28.0f, 28.0f);
        break;
    case CoordinationGeometry::TrigonalBipyramidal:
        add_site(size.x * 0.5f, 14.0f);
        add_site(size.x * 0.5f, size.y - 14.0f);
        add_site(26.0f, size.y - 28.0f);
        add_site(size.x - 26.0f, size.y - 28.0f);
        add_site(size.x * 0.5f, size.y * 0.5f + 18.0f);
        break;
    case CoordinationGeometry::SquarePyramidal:
        add_site(size.x * 0.5f, 14.0f);
        add_site(24.0f, size.y * 0.5f);
        add_site(size.x - 24.0f, size.y * 0.5f);
        add_site(size.x * 0.5f, size.y - 18.0f);
        add_site(size.x * 0.5f, size.y * 0.5f + 20.0f);
        break;
    case CoordinationGeometry::Octahedral:
        add_site(size.x * 0.5f, 14.0f);
        add_site(size.x * 0.5f, size.y - 14.0f);
        add_site(24.0f, size.y * 0.5f);
        add_site(size.x - 24.0f, size.y * 0.5f);
        add_site(42.0f, 30.0f);
        add_site(size.x - 42.0f, size.y - 30.0f);
        break;
    case CoordinationGeometry::PentagonalBipyramidal:
        add_site(size.x * 0.5f, 12.0f);
        add_site(size.x * 0.5f, size.y - 12.0f);
        for (int i = 0; i < 5; ++i) {
            const float angle = static_cast<float>(i) * 2.0f * static_cast<float>(kPi) / 5.0f - static_cast<float>(kPi) * 0.5f;
            add_site(center.x - p0.x + std::cos(angle) * 42.0f, center.y - p0.y + std::sin(angle) * 32.0f);
        }
        break;
    case CoordinationGeometry::SquareAntiprismatic:
        add_site(42.0f, 30.0f);
        add_site(size.x - 42.0f, 30.0f);
        add_site(42.0f, size.y - 30.0f);
        add_site(size.x - 42.0f, size.y - 30.0f);
        add_site(size.x * 0.5f, 18.0f);
        add_site(size.x - 24.0f, size.y * 0.5f);
        add_site(size.x * 0.5f, size.y - 18.0f);
        add_site(24.0f, size.y * 0.5f);
        break;
    }

    ImGui::Dummy(size);
}

void draw_component_bars(const LigandEntry& ligand) {
    ImVec2 p = ImGui::GetCursorScreenPos();
    const ImVec2 size(ImGui::GetContentRegionAvail().x, 10.0f);
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const ImU32 fill = ligand.charge < 0 ? IM_COL32(70, 120, 220, 220) : IM_COL32(210, 90, 90, 220);
    draw->AddRectFilled(p, ImVec2(p.x + size.x, p.y + size.y), IM_COL32(42, 46, 56, 180), 3.0f);
    const float fraction = std::min(1.0f, 0.25f + 0.15f * static_cast<float>(denticity_value(ligand.denticity)));
    draw->AddRectFilled(p, ImVec2(p.x + size.x * fraction, p.y + size.y), fill, 3.0f);
    ImGui::Dummy(ImVec2(size.x, size.y + 4.0f));
}

}  // namespace

void draw_complex_builder(AppState& state,
                          sbox::chem::MolecularSystem& mol,
                          sbox::editor::CommandStack& commands,
                          const sbox::chem::LigandLibrary& ligand_library) {
    AppState::ComplexBuilderState& builder = state.complex_builder;
    ensure_site_count(builder);

    if (!ImGui::Begin("Complex Builder")) {
        ImGui::End();
        return;
    }

    ImGui::TextUnformatted("Step 1: Select Metal");
    draw_metal_selector(builder);

    const auto& metal = sbox::elements::get_element(builder.metal_Z);
    const auto& oxidation_states = common_oxidation_states(builder.metal_Z);
    ImGui::Separator();
    ImGui::Text("Selected: %s (%s)", metal.symbol, metal.name);
    ImGui::Text("Common oxidation states: ");
    ImGui::SameLine();
    for (std::size_t i = 0; i < oxidation_states.size(); ++i) {
        if (ImGui::SmallButton((std::string(i == 0 ? "" : " ") + oxidation_string(oxidation_states[i])).c_str())) {
            builder.oxidation_state = oxidation_states[i];
        }
        if (i + 1 < oxidation_states.size()) {
            ImGui::SameLine();
        }
    }
    ImGui::InputInt("Oxidation State", &builder.oxidation_state);
    const int d_count = sbox::analysis::d_electron_count(builder.metal_Z, builder.oxidation_state);
    ImGui::Text("d-electron count: d%d (%s%s)", d_count, metal.symbol, oxidation_string(builder.oxidation_state).c_str());

    ImGui::Separator();
    ImGui::TextUnformatted("Step 2: Select Geometry");
    int geometry_index = geometry_choice_index(builder.geometry);
    if (ImGui::Combo("Geometry", &geometry_index, [](void* data, int idx, const char** out_text) {
            const auto* choices = static_cast<const std::array<GeometryChoice, 11>*>(data);
            *out_text = (*choices)[static_cast<std::size_t>(idx)].label;
            return true;
        }, const_cast<std::array<GeometryChoice, 11>*>(&kGeometryChoices), static_cast<int>(kGeometryChoices.size()))) {
        builder.geometry = kGeometryChoices[static_cast<std::size_t>(geometry_index)].geometry;
        ensure_site_count(builder);
    }
    draw_geometry_schematic(builder.geometry);

    ImGui::Separator();
    ImGui::TextUnformatted("Step 3: Assign Ligands");
    if (ImGui::Checkbox("All same ligand", &builder.all_same)) {
        if (builder.all_same) {
            const LigandEntry* ligand = ligand_library.find(builder.uniform_ligand);
            if (ligand != nullptr) {
                assign_uniform_ligand(builder, *ligand);
            }
        }
    }

    if (builder.all_same) {
        const LigandEntry* selected = ligand_library.find(builder.uniform_ligand);
        const char* preview = selected != nullptr ? selected->abbreviation.c_str() : "Select ligand";
        if (ImGui::BeginCombo("Uniform Ligand", preview)) {
            for (const std::string& category : ligand_library.categories()) {
                ImGui::TextDisabled("%s", category.c_str());
                for (const LigandEntry* ligand : ligand_library.by_category(category)) {
                    const bool chosen = ligand->name == builder.uniform_ligand;
                    const std::string label = ligand->abbreviation + " (" + ligand->name + ")";
                    if (ImGui::Selectable(label.c_str(), chosen)) {
                        builder.uniform_ligand = ligand->name;
                        assign_uniform_ligand(builder, *ligand);
                    }
                }
                ImGui::Separator();
            }
            ImGui::EndCombo();
        }
    }

    ImGui::Columns(4, "complex_builder_sites", false);
    ImGui::TextUnformatted("Site");
    ImGui::NextColumn();
    ImGui::TextUnformatted("Ligand");
    ImGui::NextColumn();
    ImGui::TextUnformatted("Info");
    ImGui::NextColumn();
    ImGui::TextUnformatted("Field");
    ImGui::NextColumn();
    ImGui::Separator();

    for (int site = 0; site < static_cast<int>(builder.site_ligands.size()); ++site) {
        int leader = -1;
        const LigandEntry* occupying = nullptr;
        const bool occupied = site_is_occupied_by_previous(builder, ligand_library, site, &leader, &occupying);

        ImGui::Text("Site %d", site + 1);
        ImGui::NextColumn();

        if (occupied && occupying != nullptr) {
            ImGui::TextDisabled("<- %s (%d-dentate)", occupying->abbreviation.c_str(), denticity_value(occupying->denticity));
            ImGui::NextColumn();
            ImGui::TextDisabled("occupied by site %d", leader + 1);
            ImGui::NextColumn();
            ImGui::TextDisabled("%s", occupying->field_strength.c_str());
            ImGui::NextColumn();
            continue;
        }

        const LigandEntry* selected = ligand_for_name(ligand_library, builder.site_ligands[static_cast<std::size_t>(site)]);
        const char* preview = selected != nullptr ? selected->abbreviation.c_str() : "Select ligand";
        ImGui::BeginDisabled(builder.all_same);
        if (ImGui::BeginCombo(("##site_" + std::to_string(site)).c_str(), preview)) {
            for (const std::string& category : ligand_library.categories()) {
                ImGui::TextDisabled("%s", category.c_str());
                for (const LigandEntry* ligand : ligand_library.by_category(category)) {
                    const bool chosen = selected != nullptr && ligand->name == selected->name;
                    const std::string label = ligand->abbreviation + " (" + ligand->name + ")";
                    if (ImGui::Selectable(label.c_str(), chosen)) {
                        assign_ligand(builder, ligand_library, site, *ligand);
                    }
                }
                ImGui::Separator();
            }
            ImGui::EndCombo();
        }
        ImGui::EndDisabled();

        ImGui::NextColumn();
        if (selected != nullptr && !selected->donor_atoms.empty()) {
            const int donor_index = selected->donor_atoms.front();
            const int donor_Z = selected->molecule.atom(donor_index).Z;
            ImGui::Text("%s donor, charge %d", sbox::elements::get_element(donor_Z).symbol, selected->charge);
        } else {
            ImGui::TextUnformatted("n/a");
        }
        ImGui::NextColumn();
        if (selected != nullptr) {
            ImGui::Text("%s", selected->field_strength.c_str());
            draw_component_bars(*selected);
        } else {
            ImGui::TextUnformatted("n/a");
        }
        ImGui::NextColumn();
    }
    ImGui::Columns(1);

    ImGui::Separator();
    ImGui::TextUnformatted("Step 4: Build");
    const std::string summary = summarize_complex(builder, ligand_library);
    int total_charge = builder.oxidation_state;
    for (int site = 0; site < static_cast<int>(builder.site_ligands.size()); ++site) {
        if (site_is_occupied_by_previous(builder, ligand_library, site)) {
            continue;
        }
        const LigandEntry* ligand = ligand_for_name(ligand_library, builder.site_ligands[static_cast<std::size_t>(site)]);
        if (ligand != nullptr) {
            total_charge += ligand->charge;
        }
    }
    const double delta_ev = expected_delta(builder, ligand_library);
    ImGui::Text("Complex: %s", summary.c_str());
    ImGui::Text("Total charge: %s", oxidation_string(total_charge).c_str());
    ImGui::Text("d-electrons: %d", d_count);
    ImGui::Text("Expected: %s (Delta_oct ~= %.2f eV)", expected_spin_text(builder, ligand_library).c_str(), delta_ev);
    ImGui::Checkbox("Auto-optimize with xTB", &builder.auto_optimize);

    const bool valid_tiling = selection_tiles_geometry(builder, ligand_library);
    if (!valid_tiling) {
        ImGui::TextColored(ImVec4(0.92f, 0.36f, 0.30f, 1.0f), "Selected denticities do not fill the coordination geometry cleanly.");
    }

    if (ImGui::Button("Build Complex") && valid_tiling) {
        std::vector<sbox::chem::LigandSpec> specs;
        specs.reserve(builder.site_ligands.size());
        for (int site = 0; site < static_cast<int>(builder.site_ligands.size()); ++site) {
            const LigandEntry* ligand = ligand_for_name(ligand_library, builder.site_ligands[static_cast<std::size_t>(site)]);
            if (ligand == nullptr) {
                continue;
            }
            specs.push_back(ligand_library.to_ligand_spec(*ligand));
        }

        sbox::chem::MolecularSystem complex = sbox::chem::assemble_complex(
            builder.metal_Z,
            builder.oxidation_state,
            builder.geometry,
            specs
        );
        complex.set_name(summary);
        complex.set_charge(total_charge);
        commands.execute(std::make_unique<sbox::editor::AddFragmentCommand>(complex, summary), mol);
        state.molecule_loaded = true;
        state.view_mode = ViewMode::MolecularOrbital;
        builder.built = true;

        if (builder.auto_optimize) {
            state.computation.method = sbox::backend::Method::GFN2_XTB;
            state.computation.charge = total_charge;
            state.computation.multiplicity = std::max(1, d_count % 2 == 0 ? 1 : 2);
            state.computation.optimize = true;
            state.computation.job_completed = false;
            state.computation.last_error.clear();
            state.computation.run_requested = true;
        }
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Step 5: Post-build options");
    ImGui::BeginDisabled(!builder.built);
    if (ImGui::Button("Optimize Geometry (xTB)")) {
        state.computation.method = sbox::backend::Method::GFN2_XTB;
        state.computation.charge = total_charge;
        state.computation.multiplicity = std::max(1, d_count % 2 == 0 ? 1 : 2);
        state.computation.optimize = true;
        state.computation.job_completed = false;
        state.computation.last_error.clear();
        state.computation.run_requested = true;
    }
    if (ImGui::Button("Run Crystal Field Analysis")) {
        state.computation.method = sbox::backend::Method::HF;
        state.computation.basis = sbox::backend::BasisSetType::STO_3G;
        state.computation.charge = total_charge;
        state.computation.multiplicity = std::max(1, d_count % 2 == 0 ? 1 : 2);
        state.computation.optimize = false;
        state.computation.job_completed = false;
        state.computation.last_error.clear();
        state.computation.run_requested = true;
        state.property_view = PropertyView::Dashboard;
    }
    if (ImGui::Button("Modify Ligands")) {
        builder.built = false;
    }
    ImGui::EndDisabled();

    ImGui::End();
}

}  // namespace sbox::ui
