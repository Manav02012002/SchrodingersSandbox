#include "ui/optimization_panel.h"

#include "core/elements.h"
#include "io/xyz_io.h"
#include "ui/file_dialog.h"
#include "ui/plot_utils.h"

#include <implot.h>
#include <imgui.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace sbox::ui {

namespace {

constexpr double kHartreeToEv = 27.2114;
constexpr double kHartreeToKcal = 627.509;
constexpr double kBohrToAngstrom = 0.529177;
constexpr double kRadToDeg = 57.29577951308232;

struct BondChange {
    std::string label;
    double initial = 0.0;
    double final = 0.0;
    double delta = 0.0;
};

struct AngleChange {
    std::string label;
    double initial = 0.0;
    double final = 0.0;
    double delta = 0.0;
};

std::string atom_label(const sbox::chem::MolecularSystem& mol, int atom_index) {
    return std::string(sbox::elements::get_element(mol.atom(atom_index).Z).symbol) + std::to_string(atom_index + 1);
}

void write_xyz_trajectory(const std::string& path, const std::vector<sbox::chem::MolecularSystem>& frames) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Could not open trajectory file for writing: " + path);
    }
    for (std::size_t i = 0; i < frames.size(); ++i) {
        const auto& mol = frames[i];
        out << mol.num_atoms() << '\n';
        out << "Frame " << i << '\n';
        for (const auto& atom : mol.atoms()) {
            out << sbox::elements::get_element(atom.Z).symbol << ' '
                << atom.position.x() * kBohrToAngstrom << ' '
                << atom.position.y() * kBohrToAngstrom << ' '
                << atom.position.z() * kBohrToAngstrom << '\n';
        }
    }
}

std::vector<BondChange> largest_bond_changes(const sbox::chem::MolecularSystem& initial,
                                             const sbox::chem::MolecularSystem& final) {
    std::vector<BondChange> changes;
    for (const auto& bond : initial.bonds()) {
        if (bond.atom_i >= final.num_atoms() || bond.atom_j >= final.num_atoms()) {
            continue;
        }
        const double d0 = initial.distance(bond.atom_i, bond.atom_j) * kBohrToAngstrom;
        const double d1 = final.distance(bond.atom_i, bond.atom_j) * kBohrToAngstrom;
        changes.push_back({
            atom_label(initial, bond.atom_i) + "-" + atom_label(initial, bond.atom_j),
            d0,
            d1,
            d1 - d0,
        });
    }
    std::sort(changes.begin(), changes.end(), [](const BondChange& a, const BondChange& b) {
        return std::abs(a.delta) > std::abs(b.delta);
    });
    if (changes.size() > 5) {
        changes.resize(5);
    }
    return changes;
}

std::vector<AngleChange> largest_angle_changes(const sbox::chem::MolecularSystem& initial,
                                               const sbox::chem::MolecularSystem& final) {
    std::vector<AngleChange> changes;
    for (int center = 0; center < initial.num_atoms(); ++center) {
        const std::vector<int> neigh = initial.neighbors(center);
        for (std::size_t i = 0; i < neigh.size(); ++i) {
            for (std::size_t j = i + 1; j < neigh.size(); ++j) {
                const int a = neigh[i];
                const int c = neigh[j];
                if (a >= final.num_atoms() || center >= final.num_atoms() || c >= final.num_atoms()) {
                    continue;
                }
                const double ang0 = initial.angle(a, center, c) * kRadToDeg;
                const double ang1 = final.angle(a, center, c) * kRadToDeg;
                changes.push_back({
                    atom_label(initial, a) + "-" + atom_label(initial, center) + "-" + atom_label(initial, c),
                    ang0,
                    ang1,
                    ang1 - ang0,
                });
            }
        }
    }
    std::sort(changes.begin(), changes.end(), [](const AngleChange& a, const AngleChange& b) {
        return std::abs(a.delta) > std::abs(b.delta);
    });
    if (changes.size() > 5) {
        changes.resize(5);
    }
    return changes;
}

void sync_frame_to_molecule(AppState::TrajectoryPlayerState& player,
                            const sbox::backend::JobResult& result,
                            sbox::chem::MolecularSystem& mol,
                            sbox::render::MolRenderer& renderer) {
    if (!result.has_trajectory || result.trajectory_frames.empty()) {
        return;
    }
    player.total_frames = static_cast<int>(result.trajectory_frames.size());
    player.current_frame = std::clamp(player.current_frame, 0, player.total_frames - 1);
    if (player.last_applied_frame == player.current_frame) {
        return;
    }
    mol = result.trajectory_frames[static_cast<std::size_t>(player.current_frame)];
    renderer.upload(mol);
    player.last_applied_frame = player.current_frame;
}

}  // namespace

void draw_optimization_panel(AppState& state,
                             const sbox::backend::JobResult& result,
                             sbox::chem::MolecularSystem& mol,
                             sbox::render::MolRenderer& renderer) {
    if (result.opt_history.empty() && !state.computation.job_running) {
        return;
    }

    if (!ImGui::Begin("Geometry Optimisation")) {
        ImGui::End();
        return;
    }

    AppState::TrajectoryPlayerState& player = state.optimization_player;

    if (!result.opt_history.empty()) {
        std::vector<float> xs(result.opt_history.size());
        std::vector<float> energies(result.opt_history.size());
        std::vector<float> gradients(result.opt_history.size());
        for (std::size_t i = 0; i < result.opt_history.size(); ++i) {
            xs[i] = static_cast<float>(result.opt_history[i].step);
            energies[i] = static_cast<float>(result.opt_history[i].energy);
            gradients[i] = static_cast<float>(std::max(result.opt_history[i].gradient_norm, 1.0e-12));
        }

        ImGui::TextUnformatted("Convergence");
        if (ImPlot::BeginPlot("Energy Convergence", ImVec2(-1.0f, 180.0f))) {
            ImPlot::SetupAxes("Step", "Energy (Hartree)");
            plot_line_styled("Energy", xs.data(), energies.data(), static_cast<int>(xs.size()), ImVec4(0.20f, 0.70f, 0.88f, 1.0f), 2.0f);
            const float final_x = xs.back();
            const float final_y = energies.back();
            ImPlot::PlotScatter("Final", &final_x, &final_y, 1,
                                {ImPlotProp_LineColor, ImVec4(0.95f, 0.78f, 0.22f, 1.0f),
                                 ImPlotProp_Marker, ImPlotMarker_Circle,
                                 ImPlotProp_MarkerSize, 6.0f});
            ImPlot::EndPlot();
        }

        if (ImPlot::BeginPlot("Gradient Convergence", ImVec2(-1.0f, 180.0f))) {
            ImPlot::SetupAxes("Step", "RMS Gradient (Hartree/Bohr)");
            ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
            plot_line_styled("Gradient", xs.data(), gradients.data(), static_cast<int>(xs.size()), ImVec4(0.90f, 0.42f, 0.32f, 1.0f), 2.0f);
            const float threshold_x[2] = {xs.front(), xs.back()};
            const float threshold_y[2] = {1.0e-5f, 1.0e-5f};
            plot_line_styled("Threshold", threshold_x, threshold_y, 2,
                             gradients.back() <= 1.0e-5 ? ImVec4(0.20f, 0.80f, 0.40f, 1.0f) : ImVec4(0.85f, 0.25f, 0.25f, 1.0f),
                             1.5f);
            ImPlot::EndPlot();
        }
    } else if (state.computation.job_running) {
        ImGui::TextUnformatted("Optimisation running...");
    }

    if (result.has_trajectory && !result.trajectory_frames.empty()) {
        ImGui::Separator();
        ImGui::TextUnformatted("Trajectory Playback");

        player.total_frames = static_cast<int>(result.trajectory_frames.size());
        player.current_frame = std::clamp(player.current_frame, 0, player.total_frames - 1);

        if (player.playing) {
            player.accumulator += ImGui::GetIO().DeltaTime * player.playback_speed;
            constexpr float frame_duration = 0.1f;
            while (player.accumulator >= frame_duration) {
                player.accumulator -= frame_duration;
                ++player.current_frame;
                if (player.current_frame >= player.total_frames) {
                    if (player.looping) {
                        player.current_frame = 0;
                    } else {
                        player.current_frame = player.total_frames - 1;
                        player.playing = false;
                        break;
                    }
                }
            }
        }

        bool frame_changed = false;
        if (ImGui::Button("|<")) {
            player.current_frame = 0;
            frame_changed = true;
        }
        ImGui::SameLine();
        if (ImGui::Button("<")) {
            player.current_frame = std::max(0, player.current_frame - 1);
            frame_changed = true;
        }
        ImGui::SameLine();
        if (ImGui::Button(player.playing ? "Pause" : "Play")) {
            player.playing = !player.playing;
        }
        ImGui::SameLine();
        if (ImGui::Button(">")) {
            player.current_frame = std::min(player.total_frames - 1, player.current_frame + 1);
            frame_changed = true;
        }
        ImGui::SameLine();
        if (ImGui::Button(">|")) {
            player.current_frame = player.total_frames - 1;
            frame_changed = true;
        }
        ImGui::SameLine();
        ImGui::Checkbox("Loop", &player.looping);

        if (ImGui::SliderInt("Step", &player.current_frame, 0, player.total_frames - 1)) {
            frame_changed = true;
        }
        ImGui::SliderFloat("Speed", &player.playback_speed, 0.1f, 5.0f, "%.1fx");

        if (frame_changed || player.last_applied_frame != player.current_frame) {
            sync_frame_to_molecule(player, result, mol, renderer);
        }

        ImGui::Separator();
        const int step = std::min(player.current_frame, static_cast<int>(result.opt_history.size()) - 1);
        const double current_energy = step >= 0 && !result.opt_history.empty() ? result.opt_history[static_cast<std::size_t>(step)].energy : result.total_energy;
        const double current_grad = step >= 0 && !result.opt_history.empty() ? result.opt_history[static_cast<std::size_t>(step)].gradient_norm : 0.0;
        const double start_energy = !result.opt_history.empty() ? result.opt_history.front().energy : current_energy;

        ImGui::Text("Step: %d / %d", player.current_frame + 1, player.total_frames);
        ImGui::Text("Energy: %.6f Hartree (%.2f eV)", current_energy, current_energy * kHartreeToEv);
        ImGui::Text("dE from start: %.6f Hartree (%.2f kcal/mol)",
                    current_energy - start_energy,
                    (current_energy - start_energy) * kHartreeToKcal);
        ImGui::Text("RMS Gradient: %.6e Hartree/Bohr", current_grad);
        if (player.current_frame == player.total_frames - 1 && result.optimization_converged) {
            ImGui::TextColored(ImVec4(0.24f, 0.82f, 0.42f, 1.0f), "Optimisation converged");
            ImGui::Text("Energy change: %.2f kcal/mol from starting geometry",
                        (current_energy - start_energy) * kHartreeToKcal);
        }

        if (result.trajectory_frames.size() >= 2) {
            ImGui::Separator();
            const auto bond_changes = largest_bond_changes(result.trajectory_frames.front(), result.trajectory_frames.back());
            const auto angle_changes = largest_angle_changes(result.trajectory_frames.front(), result.trajectory_frames.back());

            ImGui::TextUnformatted("Largest bond length changes");
            if (ImGui::BeginTable("OptBondChanges", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                ImGui::TableSetupColumn("Bond");
                ImGui::TableSetupColumn("Initial (A)");
                ImGui::TableSetupColumn("Final (A)");
                ImGui::TableSetupColumn("d (A)");
                ImGui::TableHeadersRow();
                for (const auto& change : bond_changes) {
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("%s", change.label.c_str());
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%.3f", change.initial);
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.3f", change.final);
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%+.3f", change.delta);
                }
                ImGui::EndTable();
            }

            ImGui::TextUnformatted("Largest angle changes");
            if (ImGui::BeginTable("OptAngleChanges", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                ImGui::TableSetupColumn("Angle");
                ImGui::TableSetupColumn("Initial (deg)");
                ImGui::TableSetupColumn("Final (deg)");
                ImGui::TableSetupColumn("d (deg)");
                ImGui::TableHeadersRow();
                for (const auto& change : angle_changes) {
                    ImGui::TableNextRow();
                    ImGui::TableSetColumnIndex(0); ImGui::Text("%s", change.label.c_str());
                    ImGui::TableSetColumnIndex(1); ImGui::Text("%.2f", change.initial);
                    ImGui::TableSetColumnIndex(2); ImGui::Text("%.2f", change.final);
                    ImGui::TableSetColumnIndex(3); ImGui::Text("%+.2f", change.delta);
                }
                ImGui::EndTable();
            }
        }

        ImGui::Separator();
        if (ImGui::Button("Export Trajectory (XYZ)")) {
            const std::string path = save_file_dialog("Export Trajectory", "xyz", "trajectory.xyz");
            if (!path.empty()) {
                write_xyz_trajectory(path, result.trajectory_frames);
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Export Optimised Geometry (XYZ)")) {
            const std::string path = save_file_dialog("Export Optimised Geometry", "xyz", "optimized.xyz");
            if (!path.empty()) {
                const auto& final_mol = result.trajectory_frames.back();
                sbox::io::write_xyz(path, final_mol);
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Apply Optimised Geometry")) {
            mol = result.trajectory_frames.back();
            renderer.upload(mol);
            player.current_frame = player.total_frames - 1;
            player.last_applied_frame = player.current_frame;
        }
    }

    ImGui::End();
}

}  // namespace sbox::ui
