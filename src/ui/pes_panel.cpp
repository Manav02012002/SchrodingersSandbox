#include "ui/pes_panel.h"

#include "backend/job_types.h"
#include "core/elements.h"
#include "io/xyz_io.h"
#include "ui/file_dialog.h"
#include "ui/plot_utils.h"

#include <json.hpp>

#include <implot.h>
#include <imgui.h>

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace sbox::ui {

namespace {

using json = nlohmann::json;

constexpr double kBohrToAngstrom = 0.529177;
constexpr double kHartreeToKcal = 627.509;
constexpr double kRadToDeg = 57.29577951308232;
using CoordType = sbox::backend::JobSpec::ScanSpec::CoordType;

int coord_type_atom_count(CoordType type) {
    switch (type) {
    case CoordType::Distance: return 2;
    case CoordType::Angle: return 3;
    case CoordType::Dihedral: return 4;
    }
    return 2;
}

std::string atom_label(const sbox::chem::MolecularSystem& mol, int atom_index) {
    return std::string(sbox::elements::get_element(mol.atom(atom_index).Z).symbol) + std::to_string(atom_index + 1);
}

std::string coord_atoms_label(const sbox::chem::MolecularSystem& mol, const std::vector<int>& atoms) {
    std::string label;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        if (i > 0) {
            label += "-";
        }
        label += atom_label(mol, atoms[i]);
    }
    return label;
}

double current_coordinate_value(const sbox::chem::MolecularSystem& mol,
                                CoordType type,
                                const std::vector<int>& atoms) {
    if (type == CoordType::Distance && atoms.size() >= 2) {
        return mol.distance(atoms[0], atoms[1]) * kBohrToAngstrom;
    }
    if (type == CoordType::Angle && atoms.size() >= 3) {
        return mol.angle(atoms[0], atoms[1], atoms[2]) * kRadToDeg;
    }
    if (type == CoordType::Dihedral && atoms.size() >= 4) {
        return mol.dihedral(atoms[0], atoms[1], atoms[2], atoms[3]) * kRadToDeg;
    }
    return 0.0;
}

std::string current_value_label(const sbox::chem::MolecularSystem& mol,
                                CoordType type,
                                const std::vector<int>& atoms) {
    char buffer[128];
    if (type == CoordType::Distance && atoms.size() >= 2) {
        std::snprintf(buffer,
                      sizeof(buffer),
                      "Current: d(%s) = %.3f A",
                      coord_atoms_label(mol, atoms).c_str(),
                      current_coordinate_value(mol, type, atoms));
        return buffer;
    }
    if (type == CoordType::Angle && atoms.size() >= 3) {
        std::snprintf(buffer,
                      sizeof(buffer),
                      "Current: angle(%s) = %.2f deg",
                      coord_atoms_label(mol, atoms).c_str(),
                      current_coordinate_value(mol, type, atoms));
        return buffer;
    }
    if (type == CoordType::Dihedral && atoms.size() >= 4) {
        std::snprintf(buffer,
                      sizeof(buffer),
                      "Current: dihedral(%s) = %.2f deg",
                      coord_atoms_label(mol, atoms).c_str(),
                      current_coordinate_value(mol, type, atoms));
        return buffer;
    }
    return {};
}

std::optional<json> load_partial_scan_json(const sbox::backend::BackendManager& backend, int job_id) {
    const std::string work_dir = backend.work_dir(job_id);
    if (work_dir.empty()) {
        return std::nullopt;
    }
    const std::filesystem::path result_path = std::filesystem::path(work_dir) / "result.json";
    if (!std::filesystem::exists(result_path)) {
        return std::nullopt;
    }
    std::ifstream input(result_path);
    if (!input) {
        return std::nullopt;
    }
    json parsed;
    input >> parsed;
    return parsed;
}

std::string default_csv_name(bool is_2d) {
    return is_2d ? "pes_scan_2d.csv" : "pes_scan_1d.csv";
}

std::string default_xyz_name() {
    return "pes_scan.xyz";
}

void write_scan_csv(const std::string& path, const AppState::PESScanState& pes) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Could not open PES CSV for writing");
    }
    if (pes.result.is_2d) {
        out << "coord1,coord2,energy_hartree,relative_energy_kcalmol\n";
        if (pes.result.energies.empty()) {
            return;
        }
        const double min_energy = *std::min_element(pes.result.energies.begin(), pes.result.energies.end());
        for (std::size_t i = 0; i < pes.result.energies.size(); ++i) {
            const double rel = (pes.result.energies[i] - min_energy) * kHartreeToKcal;
            out << pes.result.coord1_values[i] << ','
                << (i < pes.result.coord2_values.size() ? pes.result.coord2_values[i] : 0.0) << ','
                << pes.result.energies[i] << ','
                << rel << '\n';
        }
    } else {
        out << "coord1,energy_hartree,relative_energy_kcalmol\n";
        if (pes.result.energies.empty()) {
            return;
        }
        const double min_energy = *std::min_element(pes.result.energies.begin(), pes.result.energies.end());
        for (std::size_t i = 0; i < pes.result.energies.size(); ++i) {
            const double rel = (pes.result.energies[i] - min_energy) * kHartreeToKcal;
            out << pes.result.coord1_values[i] << ','
                << pes.result.energies[i] << ','
                << rel << '\n';
        }
    }
}

void write_scan_xyz(const std::string& path, const std::vector<sbox::chem::MolecularSystem>& geometries) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Could not open PES XYZ for writing");
    }
    for (std::size_t frame = 0; frame < geometries.size(); ++frame) {
        const auto& mol = geometries[frame];
        out << mol.num_atoms() << '\n';
        out << "Scan frame " << frame << '\n';
        for (const auto& atom : mol.atoms()) {
            out << sbox::elements::get_element(atom.Z).symbol << ' '
                << atom.position.x() * kBohrToAngstrom << ' '
                << atom.position.y() * kBohrToAngstrom << ' '
                << atom.position.z() * kBohrToAngstrom << '\n';
        }
    }
}

void draw_coordinate_editor(const char* label_prefix,
                            sbox::backend::JobSpec::ScanSpec::ScanCoordinate& coord,
                            bool& coord_set,
                            const sbox::chem::MolecularSystem& mol,
                            const sbox::editor::Selection& selection,
                            int default_steps) {
    ImGui::SeparatorText(label_prefix);

    int coord_type = static_cast<int>(coord.type);
    const char* types[] = {"Distance", "Angle", "Dihedral"};
    if (ImGui::Combo((std::string("Type##") + label_prefix).c_str(), &coord_type, types, IM_ARRAYSIZE(types))) {
        coord.type = static_cast<CoordType>(coord_type);
        coord_set = false;
    }

    const int required_atoms = coord_type_atom_count(coord.type);
    std::string selected_label = "Selected atoms: ";
    if (selection.atoms.empty()) {
        selected_label += "(none)";
    } else {
        selected_label += coord_atoms_label(mol, selection.atoms);
    }
    ImGui::TextWrapped("%s", selected_label.c_str());

    const bool selection_valid = static_cast<int>(selection.atoms.size()) == required_atoms;
    if (selection_valid) {
        ImGui::TextUnformatted(current_value_label(mol, coord.type, selection.atoms).c_str());
        if (ImGui::Button((std::string("Use Selection##") + label_prefix).c_str())) {
            coord.atom_indices = selection.atoms;
            coord_set = true;
            const double current = current_coordinate_value(mol, coord.type, coord.atom_indices);
            const double delta = coord.type == CoordType::Distance ? 0.5 : 30.0;
            coord.start = current - delta;
            coord.end = current + delta;
            coord.steps = std::max(default_steps, 2);
        }
    } else {
        ImGui::TextDisabled("Select %d atom%s in the viewport for this coordinate.",
                            required_atoms,
                            required_atoms == 1 ? "" : "s");
    }

    if (coord_set && static_cast<int>(coord.atom_indices.size()) == required_atoms) {
        ImGui::Text("Atoms: %s", coord_atoms_label(mol, coord.atom_indices).c_str());
    }

    float start_value = static_cast<float>(coord.start);
    float end_value = static_cast<float>(coord.end);
    int steps = std::max(coord.steps, 2);
    ImGui::InputFloat((std::string("Start##") + label_prefix).c_str(), &start_value);
    ImGui::InputFloat((std::string("End##") + label_prefix).c_str(), &end_value);
    ImGui::InputInt((std::string("Steps##") + label_prefix).c_str(), &steps);
    steps = std::max(steps, 2);
    coord.start = start_value;
    coord.end = end_value;
    coord.steps = steps;
    ImGui::SameLine();
    ImGui::TextDisabled("%s", coord.type == CoordType::Distance ? "A" : "deg");
}

std::string estimate_time_label(const AppState::PESScanState& pes, sbox::backend::Method method) {
    const int total_points = std::max(1, pes.coord1.steps) * (pes.is_2d ? std::max(1, pes.coord2.steps) : 1);
    const double seconds_per_point =
        sbox::backend::method_is_xtb(method) ? 2.0 : (method == sbox::backend::Method::HF ? 6.0 : 15.0);
    const double total_seconds = total_points * seconds_per_point;
    char buffer[128];
    std::snprintf(buffer,
                  sizeof(buffer),
                  "~%d single-point optimisations, estimated %.1f minutes.",
                  total_points,
                  total_seconds / 60.0);
    return buffer;
}

void draw_1d_results(AppState::PESScanState& pes) {
    const auto& energies = pes.result.energies;
    const auto& coords = pes.result.coord1_values;
    if (energies.empty() || coords.empty()) {
        return;
    }

    std::vector<float> xs(coords.begin(), coords.end());
    const double min_energy = *std::min_element(energies.begin(), energies.end());
    std::vector<float> rel(energies.size());
    for (std::size_t i = 0; i < energies.size(); ++i) {
        rel[i] = static_cast<float>((energies[i] - min_energy) * kHartreeToKcal);
    }

    const auto min_it = std::min_element(energies.begin(), energies.end());
    const auto max_it = std::max_element(energies.begin(), energies.end());
    const int min_idx = static_cast<int>(std::distance(energies.begin(), min_it));
    const int max_idx = static_cast<int>(std::distance(energies.begin(), max_it));

    if (ImPlot::BeginPlot("PES 1D Scan", ImVec2(-1.0f, 250.0f))) {
        ImPlot::SetupAxes("Coordinate", "Relative Energy (kcal/mol)");
        plot_line_styled("Profile", xs.data(), rel.data(), static_cast<int>(xs.size()), ImVec4(0.15f, 0.70f, 0.88f, 1.0f), 2.0f);
        ImPlot::PlotScatter("Minimum", &xs[static_cast<std::size_t>(min_idx)], &rel[static_cast<std::size_t>(min_idx)], 1,
                            {ImPlotProp_Marker, ImPlotMarker_Diamond, ImPlotProp_MarkerSize, 7.0f, ImPlotProp_LineColor, ImVec4(0.25f, 0.9f, 0.35f, 1.0f)});
        ImPlot::PlotScatter("Maximum", &xs[static_cast<std::size_t>(max_idx)], &rel[static_cast<std::size_t>(max_idx)], 1,
                            {ImPlotProp_Marker, ImPlotMarker_Square, ImPlotProp_MarkerSize, 7.0f, ImPlotProp_LineColor, ImVec4(0.95f, 0.35f, 0.25f, 1.0f)});
        if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            const ImPlotPoint mouse = ImPlot::GetPlotMousePos();
            double best = std::numeric_limits<double>::max();
            int best_idx = -1;
            for (std::size_t i = 0; i < coords.size(); ++i) {
                const double dx = coords[i] - mouse.x;
                const double dy = rel[i] - mouse.y;
                const double dist2 = dx * dx + dy * dy;
                if (dist2 < best) {
                    best = dist2;
                    best_idx = static_cast<int>(i);
                }
            }
            pes.viewed_point = best_idx;
        }
        ImPlot::EndPlot();
    }

    ImGui::Text("Minimum: %.3f, E = %.6f Hartree", coords[static_cast<std::size_t>(min_idx)], energies[static_cast<std::size_t>(min_idx)]);
    ImGui::Text("Maximum: %.3f, E = %.6f Hartree", coords[static_cast<std::size_t>(max_idx)], energies[static_cast<std::size_t>(max_idx)]);
    ImGui::Text("Barrier: %.2f kcal/mol", (*max_it - *min_it) * kHartreeToKcal);
}

void draw_2d_surface_wireframe(const AppState::PESScanState& pes,
                               const ImVec2& origin,
                               const ImVec2& size) {
    if (pes.result.energies.empty() || pes.result.steps_1 <= 1 || pes.result.steps_2 <= 1) {
        return;
    }
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const double min_e = *std::min_element(pes.result.energies.begin(), pes.result.energies.end());
    const double max_e = *std::max_element(pes.result.energies.begin(), pes.result.energies.end());
    auto project = [&](double x, double y, double e) {
        const float px = static_cast<float>((x - pes.coord1.start) / std::max(1e-9, pes.coord1.end - pes.coord1.start));
        const float py = static_cast<float>((y - pes.coord2.start) / std::max(1e-9, pes.coord2.end - pes.coord2.start));
        const float pz = static_cast<float>((e - min_e) / std::max(1e-9, max_e - min_e));
        return ImVec2(origin.x + 30.0f + (px - py) * size.x * 0.35f + size.x * 0.35f,
                      origin.y + size.y - (px + py) * size.y * 0.2f - pz * size.y * 0.45f);
    };
    for (int j = 0; j < pes.result.steps_2; ++j) {
        for (int i = 0; i < pes.result.steps_1; ++i) {
            const int idx = j * pes.result.steps_1 + i;
            if (idx >= static_cast<int>(pes.result.energies.size())) {
                continue;
            }
            const auto p = project(pes.result.coord1_values[static_cast<std::size_t>(idx)],
                                   pes.result.coord2_values[static_cast<std::size_t>(idx)],
                                   pes.result.energies[static_cast<std::size_t>(idx)]);
            const float t = static_cast<float>((pes.result.energies[static_cast<std::size_t>(idx)] - min_e) / std::max(1e-9, max_e - min_e));
            const ImU32 color = ImGui::GetColorU32(ImVec4(t, 0.2f, 1.0f - t, 0.9f));
            if (i + 1 < pes.result.steps_1) {
                const int idx_r = j * pes.result.steps_1 + (i + 1);
                if (idx_r < static_cast<int>(pes.result.energies.size())) {
                    const auto q = project(pes.result.coord1_values[static_cast<std::size_t>(idx_r)],
                                           pes.result.coord2_values[static_cast<std::size_t>(idx_r)],
                                           pes.result.energies[static_cast<std::size_t>(idx_r)]);
                    draw->AddLine(p, q, color, 1.5f);
                }
            }
            if (j + 1 < pes.result.steps_2) {
                const int idx_u = (j + 1) * pes.result.steps_1 + i;
                if (idx_u < static_cast<int>(pes.result.energies.size())) {
                    const auto q = project(pes.result.coord1_values[static_cast<std::size_t>(idx_u)],
                                           pes.result.coord2_values[static_cast<std::size_t>(idx_u)],
                                           pes.result.energies[static_cast<std::size_t>(idx_u)]);
                    draw->AddLine(p, q, color, 1.5f);
                }
            }
        }
    }
}

void draw_2d_results(AppState::PESScanState& pes) {
    if (pes.result.energies.empty() || pes.result.coord1_values.empty() || pes.result.coord2_values.empty()) {
        return;
    }

    std::vector<float> heatmap(pes.result.energies.size());
    const double min_energy = *std::min_element(pes.result.energies.begin(), pes.result.energies.end());
    const double max_energy = *std::max_element(pes.result.energies.begin(), pes.result.energies.end());
    for (std::size_t i = 0; i < pes.result.energies.size(); ++i) {
        heatmap[i] = static_cast<float>((pes.result.energies[i] - min_energy) * kHartreeToKcal);
    }

    if (ImPlot::BeginPlot("PES 2D Scan", ImVec2(-1.0f, 300.0f))) {
        ImPlot::SetupAxes("Coordinate 1", "Coordinate 2");
        ImPlot::PlotHeatmap("##pes2d",
                            heatmap.data(),
                            std::max(1, pes.result.steps_2),
                            std::max(1, pes.result.steps_1),
                            0.0,
                            (max_energy - min_energy) * kHartreeToKcal,
                            nullptr,
                            ImPlotPoint(pes.coord1.start, pes.coord2.start),
                            ImPlotPoint(pes.coord1.end, pes.coord2.end));
        if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            const ImPlotPoint mouse = ImPlot::GetPlotMousePos();
            double best = std::numeric_limits<double>::max();
            int best_idx = -1;
            for (std::size_t i = 0; i < pes.result.energies.size(); ++i) {
                const double dx = pes.result.coord1_values[i] - mouse.x;
                const double dy = pes.result.coord2_values[i] - mouse.y;
                const double dist2 = dx * dx + dy * dy;
                if (dist2 < best) {
                    best = dist2;
                    best_idx = static_cast<int>(i);
                }
            }
            pes.viewed_point = best_idx;
        }
        ImPlot::EndPlot();
    }

    const ImVec2 plot_origin = ImGui::GetCursorScreenPos();
    const ImVec2 plot_size(ImGui::GetContentRegionAvail().x, 180.0f);
    ImGui::InvisibleButton("##pes_surface_wireframe", plot_size);
    draw_2d_surface_wireframe(pes, plot_origin, plot_size);
}

}  // namespace

void draw_pes_panel(AppState& state,
                    sbox::backend::BackendManager& backend,
                    const sbox::chem::MolecularSystem& mol,
                    const sbox::editor::Selection& selection) {
    if (!ImGui::Begin("Potential Energy Surface")) {
        ImGui::End();
        return;
    }

    AppState::PESScanState& pes = state.pes;

    if (pes.scan_running && pes.scan_job_id >= 0) {
        if (const sbox::backend::JobResult* completed = backend.result(pes.scan_job_id)) {
            pes.result = completed->scan_result;
            pes.scan_running = false;
            pes.scan_complete = completed->has_scan;
            if (pes.viewed_point < 0 && !pes.result.geometries.empty()) {
                pes.viewed_point = 0;
            }
        }
    }

    ImGui::TextUnformatted("Scan Setup");
    if (ImGui::RadioButton("1D Scan", !pes.is_2d)) {
        pes.is_2d = false;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("2D Scan", pes.is_2d)) {
        pes.is_2d = true;
    }

    draw_coordinate_editor("Coordinate 1", pes.coord1, pes.coord1_set, mol, selection, 20);
    if (pes.is_2d) {
        draw_coordinate_editor("Coordinate 2", pes.coord2, pes.coord2_set, mol, selection, 10);
    }

    ImGui::SeparatorText("Method");
    ImGui::Text("Method: %s", sbox::backend::method_display_name(state.computation.method));
    if (sbox::backend::method_needs_basis(state.computation.method)) {
        ImGui::Text("Basis: %s", sbox::backend::basis_display_name(state.computation.basis));
    } else {
        ImGui::TextDisabled("Basis: not used for this method");
    }
    ImGui::TextDisabled("%s", estimate_time_label(pes, state.computation.method).c_str());

    const bool xtb_scan_unsupported = sbox::backend::method_is_xtb(state.computation.method);
    if (xtb_scan_unsupported) {
        ImGui::TextColored(ImVec4(0.95f, 0.70f, 0.25f, 1.0f), "PES scan currently uses the PySCF scan driver. Select HF or DFT.");
    }

    const bool can_run =
        mol.num_atoms() > 0 &&
        pes.coord1_set &&
        (!pes.is_2d || pes.coord2_set) &&
        !pes.scan_running &&
        !xtb_scan_unsupported;
    if (!can_run) {
        ImGui::BeginDisabled();
    }
    if (ImGui::Button("Run Scan")) {
        sbox::backend::JobSpec spec;
        spec.geometry = mol;
        spec.geometry.set_charge(state.computation.charge);
        spec.geometry.set_multiplicity(state.computation.multiplicity);
        spec.method = state.computation.method;
        spec.basis = state.computation.basis;
        spec.charge = state.computation.charge;
        spec.multiplicity = state.computation.multiplicity;
        spec.max_scf_cycles = 200;
        spec.scf_convergence = 1.0e-8;
        spec.run_pes_scan = true;
        spec.scan.is_2d = pes.is_2d;
        spec.scan.coord1 = pes.coord1;
        spec.scan.coord2 = pes.coord2;
        pes.scan_job_id = backend.submit(spec);
        pes.scan_running = true;
        pes.scan_complete = false;
        pes.result = {};
        pes.viewed_point = -1;
    }
    if (!can_run) {
        ImGui::EndDisabled();
    }

    if (pes.scan_running && pes.scan_job_id >= 0) {
        ImGui::SeparatorText("Progress");
        const auto progress = backend.get_progress(pes.scan_job_id);
        ImGui::Text("Point %d / %d", progress.step, progress.total);
        const float fraction = progress.total > 0 ? static_cast<float>(progress.step) / static_cast<float>(progress.total) : 0.0f;
        ImGui::ProgressBar(fraction, ImVec2(-1.0f, 0.0f));
        ImGui::Text("Current energy: %.6f Hartree", progress.energy);

        if (!pes.is_2d) {
            if (const auto partial = load_partial_scan_json(backend, pes.scan_job_id); partial.has_value()) {
                if (partial->contains("energies") && partial->contains("coordinate_values_1")) {
                    const auto energies = (*partial)["energies"].get<std::vector<double>>();
                    const auto coords = (*partial)["coordinate_values_1"].get<std::vector<double>>();
                    if (!energies.empty() && energies.size() == coords.size()) {
                        std::vector<float> xs(coords.begin(), coords.end());
                        std::vector<float> ys(energies.size());
                        for (std::size_t i = 0; i < energies.size(); ++i) {
                            ys[i] = static_cast<float>(energies[i]);
                        }
                        if (ImPlot::BeginPlot("Live 1D Scan", ImVec2(-1.0f, 180.0f))) {
                            ImPlot::SetupAxes("Coordinate", "Energy (Hartree)");
                            plot_line_styled("Energy", xs.data(), ys.data(), static_cast<int>(xs.size()), ImVec4(0.15f, 0.70f, 0.88f, 1.0f), 2.0f);
                            ImPlot::EndPlot();
                        }
                    }
                }
            }
        }
    }

    if (pes.scan_complete) {
        ImGui::SeparatorText(pes.result.is_2d ? "2D Scan Results" : "1D Scan Results");
        if (pes.result.is_2d) {
            draw_2d_results(pes);
        } else {
            draw_1d_results(pes);
        }

        if (!pes.result.geometries.empty()) {
            ImGui::SeparatorText("Trajectory Along Scan");
            const int max_point = static_cast<int>(pes.result.geometries.size()) - 1;
            if (pes.viewed_point < 0) {
                pes.viewed_point = 0;
            }
            ImGui::SliderInt("Point", &pes.viewed_point, 0, max_point);
            ImGui::Text("Showing point %d / %d", pes.viewed_point + 1, max_point + 1);
            if (!pes.result.energies.empty() && pes.viewed_point >= 0 && pes.viewed_point < static_cast<int>(pes.result.energies.size())) {
                ImGui::Text("Energy: %.6f Hartree", pes.result.energies[static_cast<std::size_t>(pes.viewed_point)]);
            }
        }

        ImGui::SeparatorText("Export");
        if (ImGui::Button("Export Energies (CSV)")) {
            const std::string path = save_file_dialog("Export PES CSV", "csv", default_csv_name(pes.result.is_2d).c_str());
            if (!path.empty()) {
                write_scan_csv(path, pes);
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Export Geometries (XYZ)")) {
            const std::string path = save_file_dialog("Export PES XYZ", "xyz", default_xyz_name().c_str());
            if (!path.empty()) {
                write_scan_xyz(path, pes.result.geometries);
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Export for Excel")) {
            const std::string path = save_file_dialog("Export PES Table", "csv", "pes_for_excel.csv");
            if (!path.empty()) {
                write_scan_csv(path, pes);
            }
        }
    }

    ImGui::End();
}

}  // namespace sbox::ui
