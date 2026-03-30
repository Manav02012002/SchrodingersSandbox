#include "ui/reaction_path_panel.h"

#include "backend/job_types.h"
#include "core/elements.h"
#include "io/pdb_io.h"
#include "io/sdf_io.h"
#include "io/xyz_io.h"
#include "ui/file_dialog.h"
#include "ui/plot_utils.h"

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
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace sbox::ui {

namespace {

constexpr double kBohrToAngstrom = 0.529177;
constexpr double kHartreeToKcal = 627.509;
constexpr double kHartreeToEv = 27.2114;
constexpr double kRadToDeg = 57.29577951308232;

using CoordType = sbox::backend::JobSpec::ScanSpec::CoordType;

struct CoordOption {
    std::string label;
    CoordType type = CoordType::Distance;
    std::vector<int> atoms;
};

std::string atom_label(const sbox::chem::MolecularSystem& mol, int atom_index) {
    return std::string(sbox::elements::get_element(mol.atom(atom_index).Z).symbol) + std::to_string(atom_index + 1);
}

std::string formula_string(const sbox::chem::MolecularSystem& mol) {
    std::map<int, int> counts;
    for (const auto& atom : mol.atoms()) {
        counts[atom.Z]++;
    }
    std::string formula;
    auto it_c = counts.find(6);
    if (it_c != counts.end()) {
        formula += "C";
        if (it_c->second > 1) {
            formula += std::to_string(it_c->second);
        }
        counts.erase(6);
    }
    auto it_h = counts.find(1);
    if (it_h != counts.end()) {
        formula += "H";
        if (it_h->second > 1) {
            formula += std::to_string(it_h->second);
        }
        counts.erase(1);
    }
    for (const auto& [z, count] : counts) {
        formula += sbox::elements::get_element(z).symbol;
        if (count > 1) {
            formula += std::to_string(count);
        }
    }
    return formula.empty() ? "Unknown" : formula;
}

std::optional<sbox::chem::MolecularSystem> load_molecule_any(const std::string& path) {
    const std::string ext = std::filesystem::path(path).extension().string();
    if (ext == ".xyz") {
        return sbox::io::read_xyz(path);
    }
    if (ext == ".sdf" || ext == ".mol") {
        return sbox::io::read_sdf(path);
    }
    if (ext == ".pdb" || ext == ".ent") {
        return sbox::io::read_pdb(path).to_molecular_system();
    }
    return std::nullopt;
}

bool validate_reaction_pair(const sbox::chem::MolecularSystem& reactant,
                            const sbox::chem::MolecularSystem& product,
                            std::vector<std::string>& warnings) {
    warnings.clear();
    if (reactant.num_atoms() != product.num_atoms()) {
        warnings.push_back("Atom count mismatch: reactant has " + std::to_string(reactant.num_atoms())
                           + " atoms, product has " + std::to_string(product.num_atoms()) + ".");
        return false;
    }
    bool ok = true;
    double mean_disp = 0.0;
    for (int i = 0; i < reactant.num_atoms(); ++i) {
        if (reactant.atom(i).Z != product.atom(i).Z) {
            warnings.push_back("Atom types differ at position " + std::to_string(i + 1) + ": reactant has "
                               + std::string(sbox::elements::get_element(reactant.atom(i).Z).symbol) + ", product has "
                               + std::string(sbox::elements::get_element(product.atom(i).Z).symbol) + ".");
            ok = false;
        }
        mean_disp += (reactant.atom(i).position - product.atom(i).position).norm();
    }
    if (reactant.num_atoms() > 0) {
        mean_disp /= reactant.num_atoms();
        if (mean_disp > 5.0) {
            warnings.push_back("Large displacement detected. Consider reordering atoms.");
        }
    }
    return ok;
}

std::vector<CoordOption> build_coordinate_options(const sbox::chem::MolecularSystem& mol) {
    std::vector<CoordOption> options;
    for (const auto& bond : mol.bonds()) {
        options.push_back({"d(" + atom_label(mol, bond.atom_i) + "-" + atom_label(mol, bond.atom_j) + ")",
                           CoordType::Distance,
                           {bond.atom_i, bond.atom_j}});
    }

    std::set<std::tuple<int, int, int>> seen_angles;
    for (int center = 0; center < mol.num_atoms(); ++center) {
        const auto neigh = mol.neighbors(center);
        for (std::size_t i = 0; i < neigh.size(); ++i) {
            for (std::size_t j = i + 1; j < neigh.size(); ++j) {
                auto tup = std::make_tuple(neigh[i], center, neigh[j]);
                if (seen_angles.insert(tup).second) {
                    options.push_back({"a(" + atom_label(mol, neigh[i]) + "-" + atom_label(mol, center) + "-" + atom_label(mol, neigh[j]) + ")",
                                       CoordType::Angle,
                                       {neigh[i], center, neigh[j]}});
                }
            }
        }
    }

    std::set<std::tuple<int, int, int, int>> seen_dihedrals;
    for (const auto& bond : mol.bonds()) {
        const auto left = mol.neighbors(bond.atom_i);
        const auto right = mol.neighbors(bond.atom_j);
        for (int a : left) {
            if (a == bond.atom_j) {
                continue;
            }
            for (int d : right) {
                if (d == bond.atom_i) {
                    continue;
                }
                auto tup = std::make_tuple(a, bond.atom_i, bond.atom_j, d);
                if (seen_dihedrals.insert(tup).second) {
                    options.push_back({"t(" + atom_label(mol, a) + "-" + atom_label(mol, bond.atom_i) + "-"
                                           + atom_label(mol, bond.atom_j) + "-" + atom_label(mol, d) + ")",
                                       CoordType::Dihedral,
                                       {a, bond.atom_i, bond.atom_j, d}});
                }
            }
        }
    }
    return options;
}

double coordinate_value(const sbox::chem::MolecularSystem& mol, const CoordOption& option) {
    switch (option.type) {
    case CoordType::Distance:
        return mol.distance(option.atoms[0], option.atoms[1]) * kBohrToAngstrom;
    case CoordType::Angle:
        return mol.angle(option.atoms[0], option.atoms[1], option.atoms[2]) * kRadToDeg;
    case CoordType::Dihedral:
        return mol.dihedral(option.atoms[0], option.atoms[1], option.atoms[2], option.atoms[3]) * kRadToDeg;
    }
    return 0.0;
}

void update_tracked_values(AppState::ReactionPathState& path) {
    for (auto& tracked : path.tracked) {
        tracked.values.clear();
        tracked.values.reserve(path.result.path_geometries.size());
        CoordOption option{tracked.label, tracked.type, tracked.atoms};
        for (const auto& geom : path.result.path_geometries) {
            tracked.values.push_back(coordinate_value(geom, option));
        }
    }
}

Eigen::Vector3d lerp(const Eigen::Vector3d& a, const Eigen::Vector3d& b, double t) {
    return a * (1.0 - t) + b * t;
}

sbox::chem::MolecularSystem interpolated_frame(const sbox::chem::MolecularSystem& a,
                                               const sbox::chem::MolecularSystem& b,
                                               double t) {
    sbox::chem::MolecularSystem out = a;
    for (int i = 0; i < out.num_atoms() && i < b.num_atoms(); ++i) {
        out.atom(i).position = lerp(a.atom(i).position, b.atom(i).position, t);
    }
    out.perceive_bonds();
    return out;
}

void apply_neb_view_frame(AppState::ReactionPathState& path,
                          sbox::chem::MolecularSystem& mol,
                          sbox::render::MolRenderer& renderer) {
    if (path.result.path_geometries.empty()) {
        return;
    }
    path.player.total_frames = static_cast<int>(path.result.path_geometries.size());
    path.player.current_frame = std::clamp(path.player.current_frame, 0, path.player.total_frames - 1);

    constexpr float frame_duration = 0.12f;
    float interp_t = 0.0f;
    if (path.smooth_interpolation && path.player.playing) {
        interp_t = std::clamp(path.player.accumulator / frame_duration, 0.0f, 1.0f);
    }

    const int current = path.player.current_frame;
    const int next = std::min(current + 1, path.player.total_frames - 1);
    if (path.smooth_interpolation && next != current && interp_t > 0.0f) {
        mol = interpolated_frame(path.result.path_geometries[static_cast<std::size_t>(current)],
                                 path.result.path_geometries[static_cast<std::size_t>(next)],
                                 interp_t);
    } else {
        mol = path.result.path_geometries[static_cast<std::size_t>(current)];
    }
    renderer.upload(mol);
    path.player.last_applied_frame = current;
}

void write_xyz_frames(const std::string& path, const std::vector<sbox::chem::MolecularSystem>& frames) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Could not open reaction-path XYZ file");
    }
    for (std::size_t i = 0; i < frames.size(); ++i) {
        const auto& mol = frames[i];
        out << mol.num_atoms() << '\n';
        out << "Image " << i << '\n';
        for (const auto& atom : mol.atoms()) {
            out << sbox::elements::get_element(atom.Z).symbol << ' '
                << atom.position.x() * kBohrToAngstrom << ' '
                << atom.position.y() * kBohrToAngstrom << ' '
                << atom.position.z() * kBohrToAngstrom << '\n';
        }
    }
}

void write_energy_csv(const std::string& path, const AppState::ReactionPathState& rp) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Could not open reaction-path CSV file");
    }
    out << "image,energy_hartree,relative_energy_kcalmol\n";
    if (rp.result.path_energies.empty()) {
        return;
    }
    const double ref = rp.result.path_energies.front();
    for (std::size_t i = 0; i < rp.result.path_energies.size(); ++i) {
        out << i << ','
            << rp.result.path_energies[i] << ','
            << (rp.result.path_energies[i] - ref) * kHartreeToKcal << '\n';
    }
}

void draw_spline_profile(const std::vector<float>& xs, const std::vector<float>& ys) {
    if (xs.size() < 2) {
        return;
    }
    std::vector<float> sx;
    std::vector<float> sy;
    for (std::size_t i = 0; i + 1 < xs.size(); ++i) {
        const float x0 = xs[i > 0 ? i - 1 : i];
        const float x1 = xs[i];
        const float x2 = xs[i + 1];
        const float x3 = xs[(i + 2 < xs.size()) ? i + 2 : i + 1];
        const float y0 = ys[i > 0 ? i - 1 : i];
        const float y1 = ys[i];
        const float y2 = ys[i + 1];
        const float y3 = ys[(i + 2 < ys.size()) ? i + 2 : i + 1];
        for (int s = 0; s < 12; ++s) {
            const float t = static_cast<float>(s) / 12.0f;
            const float tt = t * t;
            const float ttt = tt * t;
            const float x = 0.5f * ((2.0f * x1) + (-x0 + x2) * t + (2.0f * x0 - 5.0f * x1 + 4.0f * x2 - x3) * tt
                                    + (-x0 + 3.0f * x1 - 3.0f * x2 + x3) * ttt);
            const float y = 0.5f * ((2.0f * y1) + (-y0 + y2) * t + (2.0f * y0 - 5.0f * y1 + 4.0f * y2 - y3) * tt
                                    + (-y0 + 3.0f * y1 - 3.0f * y2 + y3) * ttt);
            sx.push_back(x);
            sy.push_back(y);
        }
    }
    sx.push_back(xs.back());
    sy.push_back(ys.back());
    plot_line_styled("Spline", sx.data(), sy.data(), static_cast<int>(sx.size()), ImVec4(0.20f, 0.75f, 0.90f, 1.0f), 2.0f);
}

}  // namespace

void draw_reaction_path_panel(AppState& state,
                              sbox::backend::BackendManager& backend,
                              sbox::chem::MolecularSystem& mol,
                              sbox::render::MolRenderer& renderer) {
    if (!ImGui::Begin("Reaction Path")) {
        ImGui::End();
        return;
    }

    auto& rp = state.reaction_path;

    if (rp.ts_frequency_running && rp.ts_frequency_job_id >= 0) {
        if (const auto* result = backend.result(rp.ts_frequency_job_id)) {
            rp.ts_frequency_result = *result;
            rp.ts_frequency_running = false;
        }
    }
    if (rp.neb_running && rp.neb_job_id >= 0) {
        if (const auto* result = backend.result(rp.neb_job_id)) {
            rp.result = result->neb_result;
            rp.neb_running = false;
            rp.neb_complete = result->has_neb;
            rp.player.current_frame = 0;
            rp.player.total_frames = static_cast<int>(rp.result.path_geometries.size());
            rp.player.last_applied_frame = -1;
            update_tracked_values(rp);
        }
    }

    ImGui::SeparatorText("Setup");
    if (ImGui::Button("Load Reactant from file...")) {
        const std::string path = open_file_dialog("Load Reactant Geometry", "xyz,sdf,mol,pdb,ent");
        if (!path.empty()) {
            if (auto loaded = load_molecule_any(path); loaded.has_value()) {
                rp.reactant = *loaded;
                rp.reactant_set = true;
            }
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Use current geometry##reactant")) {
        rp.reactant = mol;
        rp.reactant_set = mol.num_atoms() > 0;
    }
    if (rp.reactant_set) {
        ImGui::Text("Reactant loaded: %s (%d atoms)", formula_string(rp.reactant).c_str(), rp.reactant.num_atoms());
    }

    if (ImGui::Button("Load Product from file...")) {
        const std::string path = open_file_dialog("Load Product Geometry", "xyz,sdf,mol,pdb,ent");
        if (!path.empty()) {
            if (auto loaded = load_molecule_any(path); loaded.has_value()) {
                rp.product = *loaded;
                rp.product_set = true;
            }
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Build in editor##product")) {
        rp.product = mol;
        rp.product_set = mol.num_atoms() > 0;
    }
    if (rp.product_set) {
        ImGui::Text("Product loaded: %s (%d atoms)", formula_string(rp.product).c_str(), rp.product.num_atoms());
    }

    std::vector<std::string> warnings;
    const bool valid_pair = rp.reactant_set && rp.product_set && validate_reaction_pair(rp.reactant, rp.product, warnings);
    for (const auto& warning : warnings) {
        ImGui::TextColored(ImVec4(0.95f, 0.70f, 0.25f, 1.0f), "Warning: %s", warning.c_str());
    }

    ImGui::SliderInt("Number of images", &rp.num_images, 5, 21);
    ImGui::InputInt("Max NEB steps", &rp.max_neb_steps);
    rp.max_neb_steps = std::max(rp.max_neb_steps, 1);
    ImGui::Text("Method: %s", sbox::backend::method_display_name(state.computation.method));
    if (sbox::backend::method_needs_basis(state.computation.method)) {
        ImGui::Text("Basis: %s", sbox::backend::basis_display_name(state.computation.basis));
    }

    if (!valid_pair || rp.neb_running) {
        ImGui::BeginDisabled();
    }
    if (ImGui::Button("Run NEB")) {
        sbox::backend::JobSpec spec;
        spec.run_neb = true;
        spec.method = state.computation.method;
        spec.basis = state.computation.basis;
        spec.charge = state.computation.charge;
        spec.multiplicity = state.computation.multiplicity;
        spec.neb.reactant = rp.reactant;
        spec.neb.product = rp.product;
        spec.neb.num_images = rp.num_images;
        spec.neb.max_neb_steps = rp.max_neb_steps;
        rp.neb_job_id = backend.submit(spec);
        rp.neb_running = true;
        rp.neb_complete = false;
        rp.result = {};
        rp.player = {};
    }
    if (!valid_pair || rp.neb_running) {
        ImGui::EndDisabled();
    }

    if (rp.neb_running) {
        ImGui::SeparatorText("Progress");
        const auto progress = backend.get_progress(rp.neb_job_id);
        ImGui::Text("Image %d / %d", progress.step, rp.num_images);
        ImGui::ProgressBar(rp.num_images > 0 ? static_cast<float>(progress.step) / static_cast<float>(rp.num_images) : 0.0f, ImVec2(-1.0f, 0.0f));
        if (!progress.message.empty()) {
            ImGui::TextWrapped("%s", progress.message.c_str());
        }
    }

    if (rp.neb_complete && !rp.result.path_energies.empty()) {
        ImGui::SeparatorText("Energy Profile");
        std::vector<float> xs(rp.result.path_energies.size());
        std::vector<float> rel(rp.result.path_energies.size());
        const double ref = rp.result.path_energies.front();
        for (std::size_t i = 0; i < xs.size(); ++i) {
            xs[i] = static_cast<float>(i);
            rel[i] = static_cast<float>((rp.result.path_energies[i] - ref) * kHartreeToKcal);
        }

        if (ImPlot::BeginPlot("Reaction Energy Profile", ImVec2(-1.0f, 250.0f))) {
            ImPlot::SetupAxes("Reaction Coordinate", "Relative Energy (kcal/mol)");
            draw_spline_profile(xs, rel);
            ImPlot::PlotScatter("Images", xs.data(), rel.data(), static_cast<int>(xs.size()),
                                {ImPlotProp_Marker, ImPlotMarker_Circle, ImPlotProp_MarkerSize, 5.0f});

            const int ts = rp.result.ts_index;
            if (!xs.empty()) {
                const float react_y = rel.front();
                const float prod_y = rel.back();
                ImPlot::PlotScatter("Reactant", &xs.front(), &react_y, 1,
                                    {ImPlotProp_LineColor, ImVec4(0.25f, 0.85f, 0.35f, 1.0f), ImPlotProp_Marker, ImPlotMarker_Square});
                ImPlot::PlotScatter("Product", &xs.back(), &prod_y, 1,
                                    {ImPlotProp_LineColor, ImVec4(0.25f, 0.45f, 0.95f, 1.0f), ImPlotProp_Marker, ImPlotMarker_Square});
                if (ts >= 0 && ts < static_cast<int>(rel.size())) {
                    float tsx = xs[static_cast<std::size_t>(ts)];
                    float tsy = rel[static_cast<std::size_t>(ts)];
                    ImPlot::PlotScatter("TS", &tsx, &tsy, 1,
                                        {ImPlotProp_LineColor, ImVec4(0.95f, 0.25f, 0.25f, 1.0f), ImPlotProp_Marker, ImPlotMarker_Diamond, ImPlotProp_MarkerSize, 8.0f});
                }
            }

            if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                const ImPlotPoint mouse = ImPlot::GetPlotMousePos();
                double best = std::numeric_limits<double>::max();
                int best_idx = -1;
                for (std::size_t i = 0; i < xs.size(); ++i) {
                    const double dx = xs[i] - mouse.x;
                    const double dy = rel[i] - mouse.y;
                    const double dist2 = dx * dx + dy * dy;
                    if (dist2 < best) {
                        best = dist2;
                        best_idx = static_cast<int>(i);
                    }
                }
                if (best_idx >= 0) {
                    rp.player.current_frame = best_idx;
                    rp.player.last_applied_frame = -1;
                }
            }
            ImPlot::EndPlot();
        }

        if (rp.result.ts_index >= 0 && rp.result.ts_index < static_cast<int>(rp.result.path_energies.size())) {
            ImGui::Text("Forward barrier: %.2f kcal/mol", rp.result.forward_barrier * kHartreeToKcal);
            ImGui::Text("Reverse barrier: %.2f kcal/mol", rp.result.reverse_barrier * kHartreeToKcal);
            ImGui::Text("Reaction energy: %.2f kcal/mol",
                        (rp.result.path_energies.back() - rp.result.path_energies.front()) * kHartreeToKcal);
        }

        ImGui::SeparatorText("Path Playback");
        rp.player.total_frames = static_cast<int>(rp.result.path_geometries.size());
        if (rp.player.playing) {
            rp.player.accumulator += ImGui::GetIO().DeltaTime * rp.player.playback_speed;
            constexpr float frame_duration = 0.12f;
            while (rp.player.accumulator >= frame_duration) {
                rp.player.accumulator -= frame_duration;
                ++rp.player.current_frame;
                if (rp.player.current_frame >= rp.player.total_frames) {
                    if (rp.player.looping) {
                        rp.player.current_frame = 0;
                    } else {
                        rp.player.current_frame = rp.player.total_frames - 1;
                        rp.player.playing = false;
                        break;
                    }
                }
            }
        }
        if (ImGui::Button("|<##neb")) {
            rp.player.current_frame = 0;
            rp.player.last_applied_frame = -1;
        }
        ImGui::SameLine();
        if (ImGui::Button("<##neb")) {
            rp.player.current_frame = std::max(0, rp.player.current_frame - 1);
            rp.player.last_applied_frame = -1;
        }
        ImGui::SameLine();
        if (ImGui::Button(rp.player.playing ? "Pause##neb" : "Play##neb")) {
            rp.player.playing = !rp.player.playing;
        }
        ImGui::SameLine();
        if (ImGui::Button(">##neb")) {
            rp.player.current_frame = std::min(rp.player.total_frames - 1, rp.player.current_frame + 1);
            rp.player.last_applied_frame = -1;
        }
        ImGui::SameLine();
        if (ImGui::Button(">|##neb")) {
            rp.player.current_frame = rp.player.total_frames - 1;
            rp.player.last_applied_frame = -1;
        }
        ImGui::SameLine();
        ImGui::Checkbox("Loop##neb", &rp.player.looping);
        ImGui::SliderInt("Image", &rp.player.current_frame, 0, std::max(0, rp.player.total_frames - 1));
        ImGui::SliderFloat("Speed##neb", &rp.player.playback_speed, 0.1f, 5.0f, "%.1fx");
        ImGui::Checkbox("Smooth interpolation", &rp.smooth_interpolation);
        apply_neb_view_frame(rp, mol, renderer);

        ImGui::SeparatorText("Transition State Info");
        if (rp.result.ts_index >= 0) {
            const double reaction_kcal = (rp.result.path_energies.back() - rp.result.path_energies.front()) * kHartreeToKcal;
            ImGui::Text("Transition State: Image #%d", rp.result.ts_index + 1);
            ImGui::Text("TS Energy: %.6f Hartree", rp.result.ts_energy);
            ImGui::Text("Forward Barrier: %.2f kcal/mol (%.2f eV)", rp.result.forward_barrier * kHartreeToKcal, rp.result.forward_barrier * kHartreeToEv);
            ImGui::Text("Reverse Barrier: %.2f kcal/mol (%.2f eV)", rp.result.reverse_barrier * kHartreeToKcal, rp.result.reverse_barrier * kHartreeToEv);
            ImGui::Text("Reaction Energy: %.2f kcal/mol (%s)", std::abs(reaction_kcal), reaction_kcal <= 0.0 ? "exothermic" : "endothermic");

            if (ImGui::Button("View TS Geometry")) {
                rp.player.current_frame = rp.result.ts_index;
                rp.player.last_applied_frame = -1;
                apply_neb_view_frame(rp, mol, renderer);
            }
            ImGui::SameLine();
            if (!rp.ts_frequency_running) {
                if (ImGui::Button("Run Frequency at TS")) {
                    sbox::backend::JobSpec spec;
                    spec.geometry = rp.result.path_geometries[static_cast<std::size_t>(rp.result.ts_index)];
                    spec.method = state.computation.method;
                    spec.basis = state.computation.basis;
                    spec.charge = state.computation.charge;
                    spec.multiplicity = state.computation.multiplicity;
                    spec.properties = {
                        sbox::backend::PropertyRequest::Frequencies,
                        sbox::backend::PropertyRequest::MoldenFile,
                    };
                    rp.ts_frequency_job_id = backend.submit(spec);
                    rp.ts_frequency_running = true;
                    rp.ts_frequency_result = {};
                }
            } else {
                ImGui::TextDisabled("Frequency job running...");
            }
            ImGui::SameLine();
            ImGui::TextDisabled("IRC not yet implemented - coming in a future update.");

            if (!rp.ts_frequency_running && rp.ts_frequency_result.has_frequencies && !rp.ts_frequency_result.frequencies_cm1.empty()) {
                const auto min_it = std::min_element(rp.ts_frequency_result.frequencies_cm1.begin(),
                                                    rp.ts_frequency_result.frequencies_cm1.end());
                if (min_it != rp.ts_frequency_result.frequencies_cm1.end() && *min_it < 0.0) {
                    ImGui::Text("Imaginary frequency: %.1fi cm^-1", std::abs(*min_it));
                }
            }
        }

        ImGui::SeparatorText("Geometric Analysis Along Path");
        const auto options = build_coordinate_options(rp.reactant_set ? rp.reactant : mol);
        static std::array<int, 3> selected_indices = {{-1, -1, -1}};
        for (int slot = 0; slot < 3; ++slot) {
            std::string preview = selected_indices[slot] >= 0 && selected_indices[slot] < static_cast<int>(options.size())
                                      ? options[static_cast<std::size_t>(selected_indices[slot])].label
                                      : "None";
            if (ImGui::BeginCombo(("Track coordinate " + std::to_string(slot + 1)).c_str(), preview.c_str())) {
                if (ImGui::Selectable("None", selected_indices[slot] == -1)) {
                    selected_indices[slot] = -1;
                }
                for (int i = 0; i < static_cast<int>(options.size()); ++i) {
                    if (ImGui::Selectable(options[static_cast<std::size_t>(i)].label.c_str(), selected_indices[slot] == i)) {
                        selected_indices[slot] = i;
                    }
                }
                ImGui::EndCombo();
            }
        }
        rp.tracked.clear();
        for (int idx : selected_indices) {
            if (idx >= 0 && idx < static_cast<int>(options.size())) {
                AppState::ReactionPathState::TrackedCoordinate tracked;
                tracked.label = options[static_cast<std::size_t>(idx)].label;
                tracked.type = options[static_cast<std::size_t>(idx)].type;
                tracked.atoms = options[static_cast<std::size_t>(idx)].atoms;
                rp.tracked.push_back(std::move(tracked));
            }
        }
        update_tracked_values(rp);
        if (!rp.tracked.empty() && ImPlot::BeginPlot("Tracked Coordinates", ImVec2(-1.0f, 220.0f))) {
            ImPlot::SetupAxes("Image", "Value");
            std::vector<float> image_x(rp.result.path_geometries.size());
            std::iota(image_x.begin(), image_x.end(), 0.0f);
            const std::array<ImVec4, 3> colors = {
                ImVec4(0.18f, 0.72f, 0.95f, 1.0f),
                ImVec4(0.92f, 0.52f, 0.18f, 1.0f),
                ImVec4(0.38f, 0.82f, 0.42f, 1.0f),
            };
            for (std::size_t i = 0; i < rp.tracked.size(); ++i) {
                std::vector<float> values(rp.tracked[i].values.begin(), rp.tracked[i].values.end());
                plot_line_styled(rp.tracked[i].label.c_str(), image_x.data(), values.data(), static_cast<int>(values.size()), colors[i % colors.size()], 2.0f);
            }
            ImPlot::EndPlot();
        }

        ImGui::SeparatorText("Export");
        if (ImGui::Button("Export Path (XYZ)")) {
            const std::string path = save_file_dialog("Export NEB Path XYZ", "xyz", "neb_path.xyz");
            if (!path.empty()) {
                write_xyz_frames(path, rp.result.path_geometries);
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Export TS Geometry (XYZ)")) {
            const std::string path = save_file_dialog("Export TS XYZ", "xyz", "transition_state.xyz");
            if (!path.empty() && rp.result.ts_index >= 0) {
                sbox::io::write_xyz(path, rp.result.path_geometries[static_cast<std::size_t>(rp.result.ts_index)]);
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Export Energies (CSV)")) {
            const std::string path = save_file_dialog("Export NEB Energies", "csv", "neb_energies.csv");
            if (!path.empty()) {
                write_energy_csv(path, rp);
            }
        }
    }

    ImGui::End();
}

}  // namespace sbox::ui
