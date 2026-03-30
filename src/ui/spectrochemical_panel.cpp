#include "ui/spectrochemical_panel.h"

#include "analysis/crystal_field.h"
#include "chem/coordination.h"
#include "core/elements.h"

#include <imgui.h>
#include <implot.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace sbox::ui {

namespace {

struct MetalOption {
    const char* label;
    int Z;
    int oxidation_state;
};

struct SeriesLigand {
    const char* label;
    const char* key;
};

constexpr std::array<MetalOption, 8> kMetals = {{
    {"Fe2+", 26, 2},
    {"Fe3+", 26, 3},
    {"Co2+", 27, 2},
    {"Co3+", 27, 3},
    {"Ni2+", 28, 2},
    {"Cu2+", 29, 2},
    {"Cr3+", 24, 3},
    {"Mn2+", 25, 2},
}};

constexpr std::array<SeriesLigand, 13> kSeriesLigands = {{
    {"I-", "i"}, {"Br-", "br"}, {"Cl-", "cl"}, {"F-", "f"}, {"OH-", "oh"}, {"H2O", "water"},
    {"NCS-", "scn"}, {"py", "py"}, {"NH3", "nh3"}, {"en", "en"}, {"NO2-", "no2"}, {"CN-", "cn"}, {"CO", "co"},
}};

constexpr std::array<const char*, 13> kExpectedOrder = {
    "i","br","cl","f","oh","water","nh3","en","no2","cn","co","py","scn"
};

int geometry_mode(sbox::chem::CoordinationGeometry geometry) {
    return geometry == sbox::chem::CoordinationGeometry::Tetrahedral ? 1 : 0;
}

sbox::chem::CoordinationGeometry geometry_from_mode(int mode) {
    return mode == 1 ? sbox::chem::CoordinationGeometry::Tetrahedral : sbox::chem::CoordinationGeometry::Octahedral;
}

sbox::backend::BasisSetType method_basis(sbox::backend::Method method) {
    switch (method) {
    case sbox::backend::Method::HF: return sbox::backend::BasisSetType::STO_3G;
    case sbox::backend::Method::DFT_B3LYP: return sbox::backend::BasisSetType::B6_31Gd;
    default: return sbox::backend::BasisSetType::STO_3G;
    }
}

double delta_for_geometry(const sbox::analysis::DOrbitalEnergies& d, sbox::chem::CoordinationGeometry g) {
    return g == sbox::chem::CoordinationGeometry::Tetrahedral ? d.delta_tet() : d.delta_oct();
}

const char* wavelength_color(double nm) {
    if (nm < 380.0) return "UV / colorless";
    if (nm < 450.0) return "violet";
    if (nm < 485.0) return "blue";
    if (nm < 500.0) return "cyan";
    if (nm < 565.0) return "green";
    if (nm < 590.0) return "yellow";
    if (nm < 625.0) return "orange";
    if (nm <= 780.0) return "red";
    return "IR / colorless";
}

int expected_rank(const std::string& key) {
    for (int i = 0; i < static_cast<int>(kExpectedOrder.size()); ++i) {
        if (key == kExpectedOrder[static_cast<std::size_t>(i)]) {
            return i;
        }
    }
    return 1000;
}

std::string result_label(const AppState::SpectrochemicalState::LigandResult& result) {
    return result.ligand_name + " (" + std::to_string(std::round(result.delta_eV * 100.0) / 100.0) + ")";
}

sbox::analysis::DOrbitalEnergies idealized_d_orbitals(double delta_eV, sbox::chem::CoordinationGeometry geometry) {
    sbox::analysis::DOrbitalEnergies d;
    if (geometry == sbox::chem::CoordinationGeometry::Tetrahedral) {
        d.dz2 = -0.6 * delta_eV;
        d.dx2y2 = -0.6 * delta_eV;
        d.dxy = 0.4 * delta_eV;
        d.dxz = 0.4 * delta_eV;
        d.dyz = 0.4 * delta_eV;
    } else {
        d.dxy = -0.4 * delta_eV;
        d.dxz = -0.4 * delta_eV;
        d.dyz = -0.4 * delta_eV;
        d.dz2 = 0.6 * delta_eV;
        d.dx2y2 = 0.6 * delta_eV;
    }
    sbox::analysis::identify_splitting(d, geometry);
    return d;
}

}  // namespace

void draw_spectrochemical_panel(AppState& state,
                                sbox::backend::BackendManager& backend,
                                const sbox::chem::LigandLibrary& ligand_library) {
    if (!ImGui::Begin("Spectrochemical Series")) {
        ImGui::End();
        return;
    }

    AppState::SpectrochemicalState& spec = state.spectrochemical;
    if (spec.ligand_selected.size() != kSeriesLigands.size()) {
        spec.ligand_selected.assign(kSeriesLigands.size(), true);
    }

    int selected_metal = 0;
    for (int i = 0; i < static_cast<int>(kMetals.size()); ++i) {
        if (kMetals[static_cast<std::size_t>(i)].Z == spec.metal_Z &&
            kMetals[static_cast<std::size_t>(i)].oxidation_state == spec.oxidation_state) {
            selected_metal = i;
            break;
        }
    }

    if (ImGui::BeginCombo("Metal", kMetals[static_cast<std::size_t>(selected_metal)].label)) {
        for (int i = 0; i < static_cast<int>(kMetals.size()); ++i) {
            const bool selected = i == selected_metal;
            if (ImGui::Selectable(kMetals[static_cast<std::size_t>(i)].label, selected)) {
                spec.metal_Z = kMetals[static_cast<std::size_t>(i)].Z;
                spec.oxidation_state = kMetals[static_cast<std::size_t>(i)].oxidation_state;
            }
            if (selected) {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }

    int geom_mode = geometry_mode(spec.geometry);
    const char* geometry_labels[] = {"Octahedral", "Tetrahedral"};
    if (ImGui::Combo("Geometry", &geom_mode, geometry_labels, 2)) {
        spec.geometry = geometry_from_mode(geom_mode);
    }

    int method_mode = spec.method == sbox::backend::Method::HF ? 1 : (spec.method == sbox::backend::Method::DFT_B3LYP ? 2 : 0);
    const char* method_labels[] = {"GFN2-xTB", "HF/STO-3G", "DFT/B3LYP/6-31G*"};
    if (ImGui::Combo("Method", &method_mode, method_labels, 3)) {
        spec.method = method_mode == 1 ? sbox::backend::Method::HF
                     : method_mode == 2 ? sbox::backend::Method::DFT_B3LYP
                                        : sbox::backend::Method::GFN2_XTB;
    }

    ImGui::TextUnformatted("Ligands to test:");
    for (int i = 0; i < static_cast<int>(kSeriesLigands.size()); ++i) {
        bool checked = spec.ligand_selected[static_cast<std::size_t>(i)];
        if (ImGui::Checkbox(kSeriesLigands[static_cast<std::size_t>(i)].label, &checked)) {
            spec.ligand_selected[static_cast<std::size_t>(i)] = checked;
        }
        if ((i % 4) != 3 && i + 1 < static_cast<int>(kSeriesLigands.size())) {
            ImGui::SameLine();
        }
    }

    if (ImGui::Button("Run Series")) {
        spec.results.clear();
        spec.series_running = false;
        spec.series_complete = false;
        const int coordination_number = sbox::chem::get_template(spec.geometry).coordination_number;
        for (int i = 0; i < static_cast<int>(kSeriesLigands.size()); ++i) {
            if (!spec.ligand_selected[static_cast<std::size_t>(i)]) {
                continue;
            }
            const sbox::chem::LigandEntry* ligand = ligand_library.find(kSeriesLigands[static_cast<std::size_t>(i)].key);
            if (ligand == nullptr) {
                continue;
            }
            const int denticity = static_cast<int>(ligand->denticity);
            if (coordination_number % denticity != 0) {
                continue;
            }
            std::vector<sbox::chem::LigandSpec> specs;
            const int ligand_count = coordination_number / denticity;
            specs.reserve(static_cast<std::size_t>(ligand_count));
            for (int n = 0; n < ligand_count; ++n) {
                specs.push_back(ligand_library.to_ligand_spec(*ligand));
            }
            sbox::chem::MolecularSystem complex = sbox::chem::assemble_complex(
                spec.metal_Z, spec.oxidation_state, spec.geometry, specs);
            complex.set_name(std::string(sbox::elements::get_element(spec.metal_Z).symbol) + "-" + ligand->abbreviation);

            sbox::backend::JobSpec job;
            job.geometry = complex;
            job.charge = spec.oxidation_state + ligand_count * ligand->charge;
            job.multiplicity = 1;
            job.method = spec.method;
            job.basis = method_basis(spec.method);
            job.properties = {
                sbox::backend::PropertyRequest::MoldenFile,
            };

            const int job_id = backend.submit(job);
            spec.results.push_back({ligand->abbreviation, job_id, false, 0.0, 0.0, {}, complex});
        }
        spec.series_running = !spec.results.empty();
    }

    bool all_complete = !spec.results.empty();
    for (auto& result : spec.results) {
        if (!result.complete) {
            const sbox::backend::JobResult* backend_result = backend.result(result.job_id);
            if (backend_result != nullptr && backend_result->converged() &&
                backend_result->has_mo_data && backend_result->mo_data.coefficients.cols() > 0) {
                result.energy = backend_result->total_energy;
                result.d_orbitals = sbox::analysis::extract_d_orbitals(backend_result->mo_data, result.geometry, 0);
                sbox::analysis::identify_splitting(result.d_orbitals, spec.geometry);
                result.delta_eV = delta_for_geometry(result.d_orbitals, spec.geometry);
                result.complete = true;
            } else if (backend_result != nullptr && backend_result->converged()) {
                result.energy = backend_result->total_energy;
                const sbox::chem::LigandEntry* ligand = ligand_library.find(result.ligand_name);
                if (ligand != nullptr && ligand->typical_delta_oct_eV > 0.0) {
                    result.delta_eV = ligand->typical_delta_oct_eV;
                    result.d_orbitals = idealized_d_orbitals(result.delta_eV, spec.geometry);
                    result.complete = true;
                } else if (backend_result->homo_lumo_gap_eV() > 0.0) {
                    result.delta_eV = backend_result->homo_lumo_gap_eV();
                    result.d_orbitals = idealized_d_orbitals(result.delta_eV, spec.geometry);
                    result.complete = true;
                } else {
                    result.complete = true;
                }
            } else if (backend_result == nullptr || backend.is_running(result.job_id)) {
                all_complete = false;
            } else {
                result.complete = true;
            }
        }
        if (!result.complete) {
            all_complete = false;
        }
    }
    spec.series_complete = spec.series_running && all_complete;
    if (spec.series_complete) {
        spec.series_running = false;
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Progress");
    if (ImGui::BeginTable("##spectro_jobs", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Ligand");
        ImGui::TableSetupColumn("Status");
        ImGui::TableSetupColumn("Energy");
        ImGui::TableSetupColumn("Delta (eV)");
        ImGui::TableHeadersRow();
        for (const auto& result : spec.results) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(result.ligand_name.c_str());
            ImGui::TableSetColumnIndex(1);
            if (result.complete) {
                ImGui::TextUnformatted("done");
            } else if (backend.is_running(result.job_id)) {
                ImGui::TextUnformatted("running");
            } else {
                ImGui::TextUnformatted("queued");
            }
            ImGui::TableSetColumnIndex(2);
            if (result.complete) ImGui::Text("%.6f", result.energy);
            ImGui::TableSetColumnIndex(3);
            if (result.complete) ImGui::Text("%.3f", result.delta_eV);
        }
        ImGui::EndTable();
    }

    if (spec.series_complete) {
        std::vector<AppState::SpectrochemicalState::LigandResult> completed = spec.results;
        std::sort(completed.begin(), completed.end(), [](const auto& a, const auto& b) { return a.delta_eV < b.delta_eV; });

        ImGui::Separator();
        std::string computed = "Computed: ";
        for (std::size_t i = 0; i < completed.size(); ++i) {
            computed += result_label(completed[i]);
            if (i + 1 < completed.size()) computed += " < ";
        }
        ImGui::TextWrapped("%s", computed.c_str());
        ImGui::TextWrapped("Expected: I- < Br- < Cl- < F- < OH- < H2O < NH3 < en < NO2- < CN- < CO");

        int correct_pairs = 0;
        int total_pairs = 0;
        for (std::size_t i = 1; i < completed.size(); ++i) {
            const bool correct = expected_rank(completed[i - 1].ligand_name) <= expected_rank(completed[i].ligand_name);
            ImGui::BulletText("%s %s %s", correct ? "[ok]" : "[x]", completed[i - 1].ligand_name.c_str(), completed[i].ligand_name.c_str());
            ++total_pairs;
            if (correct) ++correct_pairs;
        }
        ImGui::Text("Ordering: %d/%d pairs correct (%.0f%%)", correct_pairs, total_pairs,
                    total_pairs > 0 ? 100.0 * static_cast<double>(correct_pairs) / static_cast<double>(total_pairs) : 0.0);

        std::vector<double> xs(completed.size());
        std::vector<double> ys(completed.size());
        std::vector<double> expected(completed.size());
        for (std::size_t i = 0; i < completed.size(); ++i) {
            xs[i] = static_cast<double>(i);
            ys[i] = completed[i].delta_eV;
            const sbox::chem::LigandEntry* ligand = ligand_library.find(completed[i].ligand_name);
            expected[i] = ligand != nullptr ? ligand->typical_delta_oct_eV : 0.0;
        }

        if (ImPlot::BeginPlot("##spectro_bar", ImVec2(-1.0f, 240.0f))) {
            ImPlot::SetupAxes(nullptr, "Delta (eV)", ImPlotAxisFlags_NoTickLabels, 0);
            ImPlot::PlotBars("Computed", xs.data(), ys.data(), static_cast<int>(ys.size()), 0.6);
            ImPlot::PlotLine("Expected", xs.data(), expected.data(), static_cast<int>(expected.size()));
            for (int i = 0; i < static_cast<int>(completed.size()); ++i) {
                ImPlot::Annotation(xs[static_cast<std::size_t>(i)], ys[static_cast<std::size_t>(i)], ImVec4(1,1,1,1), ImVec2(0,0), true, "%s", completed[static_cast<std::size_t>(i)].ligand_name.c_str());
            }
            ImPlot::EndPlot();
        }

        if (ImPlot::BeginPlot("##spectro_scatter", ImVec2(-1.0f, 220.0f))) {
            ImPlot::SetupAxes("Expected Delta (eV)", "Computed Delta (eV)");
            ImPlot::PlotScatter("Delta", expected.data(), ys.data(), static_cast<int>(ys.size()));
            const double max_delta = std::max(*std::max_element(expected.begin(), expected.end()),
                                              *std::max_element(ys.begin(), ys.end()));
            const double ref_x[] = {0.0, max_delta};
            const double ref_y[] = {0.0, max_delta};
            ImPlot::PlotLine("y=x", ref_x, ref_y, 2);
            ImPlot::EndPlot();
        }

        if (!completed.empty()) {
            const auto target = std::find_if(completed.begin(), completed.end(), [](const auto& r) { return r.ligand_name == "water"; });
            const auto& chosen = target != completed.end() ? *target : completed.front();
            const double wavelength_nm = chosen.delta_eV > 1.0e-6 ? 1240.0 / chosen.delta_eV : 0.0;
            ImGui::Separator();
            ImGui::Text("Predicted absorption wavelength for [%s(%s)] : ~%.0f nm (Delta = %.2f eV)",
                        sbox::elements::get_element(spec.metal_Z).symbol,
                        chosen.ligand_name.c_str(),
                        wavelength_nm,
                        chosen.delta_eV);
            ImGui::Text("Color: %s", wavelength_color(wavelength_nm));
        }
    }

    ImGui::End();
}

}  // namespace sbox::ui
