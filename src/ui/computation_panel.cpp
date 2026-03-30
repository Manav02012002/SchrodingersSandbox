#include "ui/computation_panel.h"

#include <imgui.h>

#include <algorithm>
#include <array>
#include <cfloat>
#include <string>

namespace sbox::ui {

namespace {

using sbox::backend::BackendManager;
using sbox::backend::BasisSetType;
using sbox::backend::JobResult;
using sbox::backend::Method;

constexpr std::array<Method, 12> kMethods = {
    Method::HF,
    Method::UHF,
    Method::DFT_B3LYP,
    Method::DFT_PBE,
    Method::DFT_PBE0,
    Method::DFT_TPSS,
    Method::DFT_M06_2X,
    Method::MP2,
    Method::CCSD,
    Method::GFN2_XTB,
    Method::GFN1_XTB,
    Method::GFN_FF,
};

constexpr std::array<BasisSetType, 14> kBasisSets = {
    BasisSetType::STO_3G,
    BasisSetType::B3_21G,
    BasisSetType::B6_31G,
    BasisSetType::B6_31Gd,
    BasisSetType::B6_31Gdp,
    BasisSetType::B6_311Gdp,
    BasisSetType::cc_pVDZ,
    BasisSetType::cc_pVTZ,
    BasisSetType::cc_pVQZ,
    BasisSetType::def2_SVP,
    BasisSetType::def2_TZVP,
    BasisSetType::def2_TZVPP,
    BasisSetType::aug_cc_pVDZ,
    BasisSetType::aug_cc_pVTZ,
};

struct SolventOption {
    const char* label;
    const char* value;
};

constexpr std::array<SolventOption, 9> kSolvents = {{
    {"Gas Phase", ""},
    {"Water", "water"},
    {"DMSO", "dmso"},
    {"THF", "thf"},
    {"DCM", "dcm"},
    {"Acetonitrile", "acetonitrile"},
    {"Methanol", "methanol"},
    {"Ethanol", "ethanol"},
    {"Toluene", "toluene"},
}};

int solvent_index(const std::string& solvent) {
    for (int i = 0; i < static_cast<int>(kSolvents.size()); ++i) {
        if (solvent == kSolvents[static_cast<std::size_t>(i)].value) {
            return i;
        }
    }
    return 0;
}

void draw_status_line(const char* label, bool available) {
    ImGui::TextUnformatted(label);
    ImGui::SameLine();
    const ImVec4 color = available ? ImVec4(0.35f, 0.85f, 0.45f, 1.0f) : ImVec4(0.90f, 0.30f, 0.30f, 1.0f);
    ImGui::TextColored(color, "%s", available ? "available" : "not found");
}

double dipole_magnitude(const JobResult& result) {
    return result.dipole_moment.norm();
}

}  // namespace

void draw_computation_panel(AppState& state, BackendManager& backend) {
    if (!ImGui::Begin("Computation")) {
        ImGui::End();
        return;
    }

    AppState::ComputationState& comp = state.computation;

    if (ImGui::BeginCombo("Method", sbox::backend::method_display_name(comp.method))) {
        for (int i = 0; i < static_cast<int>(kMethods.size()); ++i) {
            const Method method = kMethods[static_cast<std::size_t>(i)];
            const bool selected = method == comp.method;
            if (ImGui::Selectable(sbox::backend::method_display_name(method), selected)) {
                comp.method = method;
            }
            if (selected) {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }

    if (sbox::backend::method_needs_basis(comp.method)) {
        if (ImGui::BeginCombo("Basis Set", sbox::backend::basis_display_name(comp.basis))) {
            for (int i = 0; i < static_cast<int>(kBasisSets.size()); ++i) {
                const BasisSetType basis = kBasisSets[static_cast<std::size_t>(i)];
                const bool selected = basis == comp.basis;
                if (ImGui::Selectable(sbox::backend::basis_display_name(basis), selected)) {
                    comp.basis = basis;
                }
                if (selected) {
                    ImGui::SetItemDefaultFocus();
                }
            }
            ImGui::EndCombo();
        }
    }

    ImGui::InputInt("Charge", &comp.charge);
    comp.charge = std::clamp(comp.charge, -10, 10);

    ImGui::InputInt("Multiplicity", &comp.multiplicity);
    comp.multiplicity = std::clamp(comp.multiplicity, 1, 10);

    ImGui::Checkbox("Optimize geometry", &comp.optimize);

    int solvent_idx = solvent_index(comp.solvent);
    if (ImGui::BeginCombo("Solvent", kSolvents[static_cast<std::size_t>(solvent_idx)].label)) {
        for (int i = 0; i < static_cast<int>(kSolvents.size()); ++i) {
            const bool selected = i == solvent_idx;
            if (ImGui::Selectable(kSolvents[static_cast<std::size_t>(i)].label, selected)) {
                comp.solvent = kSolvents[static_cast<std::size_t>(i)].value;
                solvent_idx = i;
            }
            if (selected) {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }

    ImGui::Separator();

    const bool run_disabled = !state.molecule_loaded || comp.job_running;
    if (run_disabled) {
        ImGui::BeginDisabled();
    }
    if (ImGui::Button("Run Calculation")) {
        comp.run_requested = true;
        comp.job_completed = false;
        comp.last_error.clear();
    }
    if (run_disabled) {
        ImGui::EndDisabled();
    }

    if (!state.molecule_loaded) {
        ImGui::TextDisabled("Load a molecule to enable computations.");
    }

    const BackendManager::Progress progress = comp.job_running ? backend.get_progress(comp.active_job_id) : BackendManager::Progress{};
    if (comp.job_running) {
        ImGui::Separator();
        ImGui::Text("Stage: %s", progress.stage.empty() ? "starting" : progress.stage.c_str());
        if (progress.stage == "scf") {
            ImGui::Text("Iteration: %d", progress.iteration);
            ImGui::Text("Energy: %.10f Hartree", progress.energy);
            if (progress.iteration > comp.last_progress_iteration) {
                comp.scf_plot_energies.push_back(static_cast<float>(progress.energy));
                comp.last_progress_iteration = progress.iteration;
            }
        } else if (!progress.message.empty()) {
            ImGui::TextWrapped("%s", progress.message.c_str());
        }

        float progress_fraction = 0.1f;
        if (progress.stage == "scf") {
            progress_fraction = std::clamp(progress.iteration / 200.0f, 0.0f, 1.0f);
        } else if (progress.stage == "done") {
            progress_fraction = 1.0f;
        }
        ImGui::ProgressBar(progress_fraction, ImVec2(-1.0f, 0.0f));
        if (ImGui::Button("Cancel")) {
            backend.cancel(comp.active_job_id);
        }
    }

    const JobResult* active_result = comp.active_job_id >= 0 ? backend.result(comp.active_job_id) : nullptr;
    if (comp.job_completed && active_result != nullptr) {
        ImGui::Separator();
        ImGui::Text("Total Energy: %.10f Hartree", active_result->total_energy);
        ImGui::Text("Total Energy: %.6f eV", active_result->total_energy * 27.2114);
        ImGui::Text("Dipole Magnitude: %.6f Debye", dipole_magnitude(*active_result));
        ImGui::Text("HOMO-LUMO Gap: %.6f eV", active_result->homo_lumo_gap_eV());
        ImGui::Text("SCF Iterations: %d", static_cast<int>(active_result->scf_history.size()));
        ImGui::Text("Wall Time: %.2f s", active_result->wall_time_seconds);
        if (ImGui::Button("Apply Results")) {
            comp.apply_results_requested = true;
        }

        comp.scf_plot_energies.clear();
        comp.scf_plot_energies.reserve(active_result->scf_history.size());
        for (const auto& iter : active_result->scf_history) {
            comp.scf_plot_energies.push_back(static_cast<float>(iter.energy));
        }
    }

    if (!comp.last_error.empty()) {
        ImGui::Separator();
        ImGui::TextColored(ImVec4(0.95f, 0.35f, 0.35f, 1.0f), "%s", comp.last_error.c_str());
    }

    if (!comp.scf_plot_energies.empty()) {
        ImGui::Separator();
        ImGui::PlotLines("##scf_convergence",
                         comp.scf_plot_energies.data(),
                         static_cast<int>(comp.scf_plot_energies.size()),
                         0,
                         "SCF Energy (Hartree)",
                         FLT_MAX,
                         FLT_MAX,
                         ImVec2(0.0f, 80.0f));
    }

    ImGui::Separator();
    draw_status_line("PySCF:", backend.can_run_pyscf());
    draw_status_line("xTB:", backend.can_run_xtb());

    ImGui::End();
}

}  // namespace sbox::ui
