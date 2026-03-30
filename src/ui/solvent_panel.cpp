#include "ui/solvent_panel.h"

#include <imgui.h>
#include <implot.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace sbox::ui {

namespace {

struct SolventOption {
    const char* label;
    const char* value;
    double dielectric;
};

constexpr std::array<SolventOption, 10> kSolvents = {{
    {"Gas Phase", "", 1.0},
    {"Water (e=78.4)", "water", 78.4},
    {"DMSO (e=46.7)", "dmso", 46.7},
    {"Acetonitrile (e=35.7)", "acetonitrile", 35.7},
    {"Methanol (e=32.7)", "methanol", 32.7},
    {"Ethanol (e=24.3)", "ethanol", 24.3},
    {"DCM (e=8.9)", "dcm", 8.9},
    {"THF (e=7.6)", "thf", 7.6},
    {"Toluene (e=2.4)", "toluene", 2.4},
    {"Hexane (e=1.9)", "hexane", 1.9},
}};

double energy_ev(double hartree) {
    return hartree * 27.2114;
}

double dipole_mag(const sbox::backend::JobResult& result) {
    return result.dipole_moment.norm();
}

double homo_ev(const sbox::backend::JobResult& result) {
    const int idx = result.homo_index();
    return idx >= 0 ? result.mo_data.energies(idx) * 27.2114 : 0.0;
}

double lumo_ev(const sbox::backend::JobResult& result) {
    const int idx = result.lumo_index();
    return idx >= 0 ? result.mo_data.energies(idx) * 27.2114 : 0.0;
}

const char* solvation_classification(double kcalmol) {
    if (kcalmol < -10.0) return "strongly solvated";
    if (kcalmol < -3.0) return "moderately solvated";
    if (kcalmol < 0.0) return "weakly solvated";
    return "unfavorable";
}

ImVec4 delta_color(double delta, bool stabilization_green = false) {
    if (stabilization_green) {
        if (delta < 0.0) return ImVec4(0.30f, 0.85f, 0.42f, 1.0f);
        if (delta > 0.0) return ImVec4(0.92f, 0.36f, 0.30f, 1.0f);
    }
    if (delta > 0.0) return ImVec4(0.92f, 0.36f, 0.30f, 1.0f);
    if (delta < 0.0) return ImVec4(0.30f, 0.85f, 0.42f, 1.0f);
    return ImVec4(0.85f, 0.85f, 0.85f, 1.0f);
}

void draw_progress_block(const char* label, int job_id, sbox::backend::BackendManager& backend) {
    if (job_id < 0) {
        return;
    }
    const auto progress = backend.get_progress(job_id);
    ImGui::Text("%s: %s", label, progress.stage.empty() ? "pending" : progress.stage.c_str());
    if (progress.iteration > 0) {
        ImGui::SameLine();
        ImGui::TextDisabled("(iter %d)", progress.iteration);
    }
    if (!progress.message.empty()) {
        ImGui::TextWrapped("%s", progress.message.c_str());
    }
}

}  // namespace

void draw_solvent_panel(AppState& state, sbox::backend::BackendManager& backend) {
    if (!ImGui::Begin("Solvent Effects")) {
        ImGui::End();
        return;
    }

    AppState::SolventPanelState& panel = state.solvent;
    panel.selected_solvent_index = std::clamp(panel.selected_solvent_index, 0, static_cast<int>(kSolvents.size()) - 1);
    const SolventOption& solvent = kSolvents[static_cast<std::size_t>(panel.selected_solvent_index)];

    ImGui::TextUnformatted("Solvent Selection");
    if (ImGui::BeginCombo("Solvent", solvent.label)) {
        for (int i = 0; i < static_cast<int>(kSolvents.size()); ++i) {
            const bool selected = i == panel.selected_solvent_index;
            if (ImGui::Selectable(kSolvents[static_cast<std::size_t>(i)].label, selected)) {
                panel.selected_solvent_index = i;
            }
            if (selected) {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }
    ImGui::Text("Dielectric: %.1f", solvent.dielectric);
    const bool xtb = sbox::backend::method_is_xtb(state.computation.method);
    ImGui::Text("Solvation model: %s", xtb ? "GBSA (xTB)" : "ddCOSMO (PySCF)");

    ImGui::Separator();
    ImGui::TextUnformatted("Quick Comparison");
    const bool running = (panel.gas_job_id >= 0 && backend.is_running(panel.gas_job_id)) ||
                         (panel.solvent_job_id >= 0 && backend.is_running(panel.solvent_job_id));
    ImGui::BeginDisabled(!state.molecule_loaded || running || panel.selected_solvent_index == 0);
    if (ImGui::Button("Run Gas Phase + Solvent")) {
        panel.run_requested = true;
        panel.gas_done = false;
        panel.solvent_done = false;
        panel.gas_job_id = -1;
        panel.solvent_job_id = -1;
    }
    ImGui::EndDisabled();
    if (panel.selected_solvent_index == 0) {
        ImGui::TextDisabled("Select a solvent to compare against gas phase.");
    }

    if (running) {
        draw_progress_block("Gas", panel.gas_job_id, backend);
        draw_progress_block("Solvent", panel.solvent_job_id, backend);
    }

    if (panel.gas_done && panel.solvent_done && panel.gas_result.converged() && panel.solvent_result.converged()) {
        const double gas_e_h = panel.gas_result.total_energy;
        const double solv_e_h = panel.solvent_result.total_energy;
        const double delta_h = solv_e_h - gas_e_h;
        const double delta_kcal = delta_h * 627.509;
        const double gas_homo = homo_ev(panel.gas_result);
        const double solv_homo = homo_ev(panel.solvent_result);
        const double gas_lumo = lumo_ev(panel.gas_result);
        const double solv_lumo = lumo_ev(panel.solvent_result);

        if (ImGui::BeginTable("##solvent_compare", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
            ImGui::TableSetupColumn("Property");
            ImGui::TableSetupColumn("Gas Phase");
            ImGui::TableSetupColumn("Solvent");
            ImGui::TableSetupColumn("Delta");
            ImGui::TableHeadersRow();

            auto row = [&](const char* label, double gas, double solv, double delta, const char* fmt, bool highlight = false, bool stabilization = false) {
                ImGui::TableNextRow();
                if (highlight) {
                    ImGui::TableSetBgColor(ImGuiTableBgTarget_RowBg0, IM_COL32(42, 60, 48, 120));
                }
                ImGui::TableSetColumnIndex(0);
                ImGui::TextUnformatted(label);
                ImGui::TableSetColumnIndex(1);
                ImGui::Text(fmt, gas);
                ImGui::TableSetColumnIndex(2);
                ImGui::Text(fmt, solv);
                ImGui::TableSetColumnIndex(3);
                ImGui::TextColored(delta_color(delta, stabilization), fmt, delta);
            };

            row("Total Energy (Hartree)", gas_e_h, solv_e_h, delta_h, "%.8f");
            row("Total Energy (eV)", energy_ev(gas_e_h), energy_ev(solv_e_h), energy_ev(delta_h), "%.4f");
            row("Solvation Energy (kcal/mol)", 0.0, 0.0, delta_kcal, "%.2f", true, true);
            row("Dipole Moment (D)", dipole_mag(panel.gas_result), dipole_mag(panel.solvent_result),
                dipole_mag(panel.solvent_result) - dipole_mag(panel.gas_result), "%.4f");
            row("HOMO (eV)", gas_homo, solv_homo, solv_homo - gas_homo, "%.3f");
            row("LUMO (eV)", gas_lumo, solv_lumo, solv_lumo - gas_lumo, "%.3f");
            row("Gap (eV)", panel.gas_result.homo_lumo_gap_eV(), panel.solvent_result.homo_lumo_gap_eV(),
                panel.solvent_result.homo_lumo_gap_eV() - panel.gas_result.homo_lumo_gap_eV(), "%.3f");
            ImGui::EndTable();
        }

        ImGui::Separator();
        ImGui::Text("Solvation Free Energy: %.2f kcal/mol", delta_kcal);
        ImGui::SameLine();
        ImGui::TextDisabled("(%s)", solvation_classification(delta_kcal));

        if (!panel.gas_result.mulliken_charges.empty() &&
            panel.gas_result.mulliken_charges.size() == panel.solvent_result.mulliken_charges.size()) {
            ImGui::Separator();
            ImGui::TextUnformatted("Charge Comparison");
            if (ImGui::BeginTable("##charge_compare", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY, ImVec2(0.0f, 180.0f))) {
                ImGui::TableSetupColumn("Atom");
                ImGui::TableSetupColumn("Element");
                ImGui::TableSetupColumn("Gas Phase");
                ImGui::TableSetupColumn("Solvent");
                ImGui::TableSetupColumn("Delta");
                ImGui::TableHeadersRow();
                for (std::size_t i = 0; i < panel.gas_result.mulliken_charges.size(); ++i) {
                    const double gas_q = panel.gas_result.mulliken_charges[i];
                    const double solv_q = panel.solvent_result.mulliken_charges[i];
                    const double dq = solv_q - gas_q;
                    ImGui::TableNextRow();
                    if (std::abs(dq) > 0.05) {
                        ImGui::TableSetBgColor(ImGuiTableBgTarget_RowBg0, IM_COL32(70, 54, 36, 110));
                    }
                    ImGui::TableSetColumnIndex(0);
                    ImGui::Text("%zu", i + 1);
                    ImGui::TableSetColumnIndex(1);
                    ImGui::Text("A%zu", i + 1);
                    ImGui::TableSetColumnIndex(2);
                    ImGui::Text("%.4f", gas_q);
                    ImGui::TableSetColumnIndex(3);
                    ImGui::Text("%.4f", solv_q);
                    ImGui::TableSetColumnIndex(4);
                    ImGui::TextColored(delta_color(dq), "%.4f", dq);
                }
                ImGui::EndTable();
            }
        }

        if (panel.history.size() >= 3) {
            ImGui::Separator();
            ImGui::TextUnformatted("History");
            std::vector<double> xs(panel.history.size());
            std::vector<double> ys(panel.history.size());
            std::vector<const char*> labels(panel.history.size());
            for (std::size_t i = 0; i < panel.history.size(); ++i) {
                xs[i] = static_cast<double>(i);
                ys[i] = panel.history[i].solvation_energy_kcalmol;
                labels[i] = panel.history[i].solvent_name.c_str();
            }
            if (ImPlot::BeginPlot("##solvent_history", ImVec2(-1.0f, 220.0f))) {
                ImPlot::SetupAxes(nullptr, "Solvation Energy (kcal/mol)", ImPlotAxisFlags_NoTickLabels, 0);
                ImPlot::PlotBars("Delta G_solv", xs.data(), ys.data(), static_cast<int>(ys.size()), 0.55);
                for (int i = 0; i < static_cast<int>(labels.size()); ++i) {
                    ImPlot::Annotation(xs[static_cast<std::size_t>(i)], ys[static_cast<std::size_t>(i)], ImVec4(1,1,1,1), ImVec2(0,0), true, "%s", labels[static_cast<std::size_t>(i)]);
                }
                ImPlot::EndPlot();
            }
        }
    }

    ImGui::End();
}

}  // namespace sbox::ui
