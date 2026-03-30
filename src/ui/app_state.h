#pragma once

#include "backend/job_types.h"
#include "core/elements.h"
#include "core/slater.h"

#include <algorithm>
#include <string>
#include <vector>

namespace sbox::ui {

enum class ViewMode : int {
    AtomicOrbital = 0,
    MolecularOrbital = 1,
};

enum class PropertyView : int {
    Dashboard = 0,
    Population,
    BondOrders,
    MoDiagram,
    DOS,
    OrbitalComposition,
    IRSpectrum,
    ESPControls,
};

struct AppState {
    struct ComputationState {
        sbox::backend::Method method = sbox::backend::Method::HF;
        sbox::backend::BasisSetType basis = sbox::backend::BasisSetType::STO_3G;
        int charge = 0;
        int multiplicity = 1;
        bool optimize = false;
        std::string solvent;

        int active_job_id = -1;
        bool job_running = false;
        bool job_completed = false;
        std::string last_error;

        bool run_requested = false;
        bool apply_results_requested = false;
        int last_progress_iteration = 0;
        std::vector<float> scf_plot_energies;
    };

    int selected_Z = 1;
    int selected_orbital_index = -1;  // -1 means "last orbital" (valence)
    int selected_m = 0;
    bool needs_update = true;
    bool molecule_loaded = false;
    bool show_charges = false;
    bool show_bond_orders = false;
    bool show_dipole = false;
    bool show_esp_surface = false;
    float esp_density_iso = 0.005f;
    float esp_color_min = -0.05f;
    float esp_color_max = 0.05f;
    bool esp_auto_range = true;

    int render_mode = 0;  // 0=volume, 1=isosurface, 2=phase isosurface
    float iso_value = 0.01f;
    float gamma = 0.4f;
    ViewMode view_mode = ViewMode::AtomicOrbital;
    PropertyView property_view = PropertyView::Dashboard;
    int selected_mo = -1;
    int selected_vibrational_mode = -1;
    bool animate_vibrational_mode = false;
    int num_mo = 0;
    int homo_index = -1;
    float mol_bound_radius = 10.0f;
    int mol_render_mode = 0;
    int mol_num_basis = 0;
    double mol_total_energy_h = 0.0;
    double mol_homo_lumo_gap_ev = 0.0;
    bool mol_has_mo_summary = false;
    ComputationState computation;

    int current_n = 1;
    int current_l = 0;
    float current_Zeff = 1.0f;

    void update() {
        selected_Z = std::clamp(selected_Z, 1, 118);

        const auto& element = sbox::elements::get_element(selected_Z);
        const auto& config = element.config;

        if (config.empty()) {
            selected_orbital_index = -1;
            selected_m = 0;
            current_n = 1;
            current_l = 0;
            current_Zeff = 1.0f;
            needs_update = false;
            return;
        }

        if (selected_orbital_index < 0 || selected_orbital_index >= static_cast<int>(config.size())) {
            selected_orbital_index = static_cast<int>(config.size()) - 1;
        }

        const auto& subshell = config[static_cast<std::size_t>(selected_orbital_index)];
        const int n = subshell.n;
        const int l = subshell.l;
        const double zeff = sbox::slater::compute_zeff(selected_Z, config, n, l);

        selected_m = std::clamp(selected_m, -l, l);

        current_n = n;
        current_l = l;
        current_Zeff = static_cast<float>(zeff);
        needs_update = false;
    }
};

}  // namespace sbox::ui
