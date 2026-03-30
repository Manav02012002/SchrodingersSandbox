#include "ui/results_panel.h"

#include "core/elements.h"
#include "io/cube_io.h"
#include "io/xyz_io.h"
#include "ui/file_dialog.h"

#include <imgui.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <filesystem>
#include <numeric>
#include <string>
#include <vector>

namespace sbox::ui {

namespace {

constexpr double kHartreeToEv = 27.2114;

ImVec4 charge_text_color(double charge) {
    if (charge > 0.01) {
        return ImVec4(0.92f, 0.35f, 0.35f, 1.0f);
    }
    if (charge < -0.01) {
        return ImVec4(0.35f, 0.55f, 0.95f, 1.0f);
    }
    return ImVec4(0.80f, 0.80f, 0.80f, 1.0f);
}

bool export_existing_file(const std::filesystem::path& source,
                          const char* filters,
                          const char* default_name) {
    if (!std::filesystem::exists(source)) {
        return false;
    }
    const std::string path = save_file_dialog("Export", filters, default_name);
    if (path.empty()) {
        return false;
    }
    std::filesystem::copy_file(source, path, std::filesystem::copy_options::overwrite_existing);
    return true;
}

}  // namespace

void draw_results_panel(const AppState& state,
                        const sbox::backend::JobResult& result,
                        const sbox::chem::MolecularSystem& mol) {
    (void)state;
    if (!result.converged()) {
        return;
    }

    if (!ImGui::Begin("Results")) {
        ImGui::End();
        return;
    }

    ImGui::Text("Total Energy: %.6f Hartree  (%.2f eV)", result.total_energy, result.total_energy * kHartreeToEv);

    const int homo = result.homo_index();
    const int lumo = result.lumo_index();
    if (result.has_mo_data && homo >= 0 && homo < result.mo_data.energies.size()) {
        const double homo_e = result.mo_data.energies(homo);
        ImGui::Text("HOMO Energy: %.4f Hartree  (%.2f eV)", homo_e, homo_e * kHartreeToEv);
    }
    if (result.has_mo_data && lumo >= 0 && lumo < result.mo_data.energies.size()) {
        const double lumo_e = result.mo_data.energies(lumo);
        ImGui::Text("LUMO Energy: %.4f Hartree  (%.2f eV)", lumo_e, lumo_e * kHartreeToEv);
    }
    if (homo >= 0 && lumo >= 0 && lumo < result.mo_data.energies.size()) {
        const double gap_h = result.mo_data.energies(lumo) - result.mo_data.energies(homo);
        ImGui::Text("HOMO-LUMO Gap: %.4f Hartree  (%.2f eV)", gap_h, gap_h * kHartreeToEv);
    }

    if (!result.mulliken_charges.empty()) {
        ImGui::Separator();
        ImGui::TextUnformatted("Mulliken Charges:");
        if (ImGui::BeginTable("MullikenCharges", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
            ImGui::TableSetupColumn("Atom");
            ImGui::TableSetupColumn("Element");
            ImGui::TableSetupColumn("Charge");
            ImGui::TableHeadersRow();

            for (int i = 0; i < mol.num_atoms() && i < static_cast<int>(result.mulliken_charges.size()); ++i) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::Text("%d", i + 1);
                ImGui::TableSetColumnIndex(1);
                ImGui::Text("%s", sbox::elements::get_element(mol.atom(i).Z).symbol);
                ImGui::TableSetColumnIndex(2);
                ImGui::TextColored(charge_text_color(result.mulliken_charges[static_cast<std::size_t>(i)]),
                                   "%+.4f",
                                   result.mulliken_charges[static_cast<std::size_t>(i)]);
            }
            ImGui::EndTable();
        }
        const double total_charge = std::accumulate(result.mulliken_charges.begin(), result.mulliken_charges.end(), 0.0);
        ImGui::Text("Total charge: %.4f", total_charge);
    }

    ImGui::Separator();
    ImGui::Text("Dipole Moment: %.2f Debye", result.dipole_moment.norm());
    ImGui::Text("  x: %.2f  y: %.2f  z: %.2f",
                result.dipole_moment.x(),
                result.dipole_moment.y(),
                result.dipole_moment.z());

    ImGui::Separator();
    ImGui::Text("SCF converged in %d iterations", static_cast<int>(result.scf_history.size()));
    ImGui::Text("Wall time: %.2f seconds", result.wall_time_seconds);
    if (!result.scf_history.empty()) {
        std::vector<float> energies;
        energies.reserve(result.scf_history.size());
        for (const auto& iter : result.scf_history) {
            energies.push_back(static_cast<float>(iter.energy));
        }
        ImGui::PlotLines("##result_scf_energy",
                         energies.data(),
                         static_cast<int>(energies.size()),
                         0,
                         "SCF Energy (Hartree)",
                         FLT_MAX,
                         FLT_MAX,
                         ImVec2(0.0f, 100.0f));

        std::vector<float> delta_e;
        delta_e.reserve(result.scf_history.size());
        for (const auto& iter : result.scf_history) {
            delta_e.push_back(static_cast<float>(std::log10(std::max(std::abs(iter.delta_energy), 1e-16))));
        }
        if (!delta_e.empty()) {
            ImGui::PlotLines("##result_scf_delta",
                             delta_e.data(),
                             static_cast<int>(delta_e.size()),
                             0,
                             "log10(|dE|)",
                             FLT_MAX,
                             FLT_MAX,
                             ImVec2(0.0f, 80.0f));
        }
    }

    if (result.mayer_bond_orders.size() > 0) {
        ImGui::Separator();
        ImGui::TextUnformatted("Bond Orders:");
        if (ImGui::BeginTable("BondOrders", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
            ImGui::TableSetupColumn("Bond");
            ImGui::TableSetupColumn("Atoms");
            ImGui::TableSetupColumn("Mayer Order");
            ImGui::TableHeadersRow();
            for (const sbox::chem::Bond& bond : mol.bonds()) {
                const double order = result.mayer_bond_orders(bond.atom_i, bond.atom_j);
                if (order <= 0.3) {
                    continue;
                }
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::Text("%d-%d", bond.atom_i + 1, bond.atom_j + 1);
                ImGui::TableSetColumnIndex(1);
                ImGui::Text("%s%d - %s%d",
                            sbox::elements::get_element(mol.atom(bond.atom_i).Z).symbol,
                            bond.atom_i + 1,
                            sbox::elements::get_element(mol.atom(bond.atom_j).Z).symbol,
                            bond.atom_j + 1);
                ImGui::TableSetColumnIndex(2);
                ImGui::Text("%.3f", order);
            }
            ImGui::EndTable();
        }
    }

    ImGui::Separator();
    const std::filesystem::path work_dir = result.work_dir;
    const std::filesystem::path molden_path = work_dir / "result.molden";
    const std::filesystem::path homo_cube_path = work_dir / "homo.cube";
    const std::filesystem::path density_cube_path = work_dir / "density.cube";

    const bool have_molden = std::filesystem::exists(molden_path);
    const bool have_cube = std::filesystem::exists(homo_cube_path) || std::filesystem::exists(density_cube_path);

    if (!have_molden) {
        ImGui::BeginDisabled();
    }
    if (ImGui::Button("Export Molden")) {
        export_existing_file(molden_path, "molden", "result.molden");
    }
    if (!have_molden) {
        ImGui::EndDisabled();
    }

    ImGui::SameLine();
    if (ImGui::Button("Export XYZ")) {
        const std::string path = save_file_dialog("Export XYZ", "xyz", "result.xyz");
        if (!path.empty()) {
            const sbox::chem::MolecularSystem& out_mol = result.has_optimized_geometry ? result.optimized_geometry : mol;
            sbox::io::write_xyz(path, out_mol);
        }
    }

    ImGui::SameLine();
    if (!have_cube) {
        ImGui::BeginDisabled();
    }
    if (ImGui::Button("Export Cube")) {
        const std::string path = save_file_dialog("Export Cube", "cube", "result.cube");
        if (!path.empty()) {
            if (result.has_homo_cube) {
                sbox::io::write_cube(path, result.homo_cube);
            } else if (result.has_density_cube) {
                sbox::io::write_cube(path, result.density_cube);
            } else if (std::filesystem::exists(homo_cube_path)) {
                std::filesystem::copy_file(homo_cube_path, path, std::filesystem::copy_options::overwrite_existing);
            } else if (std::filesystem::exists(density_cube_path)) {
                std::filesystem::copy_file(density_cube_path, path, std::filesystem::copy_options::overwrite_existing);
            }
        }
    }
    if (!have_cube) {
        ImGui::EndDisabled();
    }

    ImGui::SameLine();
    if (ImGui::Button("Copy Energy")) {
        const std::string energy_text = std::to_string(result.total_energy);
        ImGui::SetClipboardText(energy_text.c_str());
    }

    ImGui::End();
}

}  // namespace sbox::ui
