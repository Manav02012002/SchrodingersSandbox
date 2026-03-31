#pragma once

#include "analysis/crystal_field.h"
#include "backend/job_types.h"
#include "chem/coordination.h"
#include "core/elements.h"
#include "core/slater.h"
#include "ui/symmetry_overlay.h"

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
    NCI,
};

struct AppState {
    struct ComplexBuilderState {
        int metal_Z = 26;
        int oxidation_state = 2;
        sbox::chem::CoordinationGeometry geometry = sbox::chem::CoordinationGeometry::Octahedral;
        std::vector<std::string> site_ligands;
        bool all_same = true;
        std::string uniform_ligand = "water";
        bool auto_optimize = true;
        bool built = false;
    };

    struct SolventResult {
        std::string solvent_name;
        double dielectric = 0.0;
        double total_energy = 0.0;
        double solvation_energy_kcalmol = 0.0;
        double dipole = 0.0;
        double homo_eV = 0.0;
        double lumo_eV = 0.0;
    };

    struct SolventPanelState {
        int selected_solvent_index = 0;
        int gas_job_id = -1;
        int solvent_job_id = -1;
        bool gas_done = false;
        bool solvent_done = false;
        bool run_requested = false;
        sbox::backend::JobResult gas_result;
        sbox::backend::JobResult solvent_result;
        std::vector<SolventResult> history;
    };

    struct SpectrochemicalState {
        struct LigandResult {
            std::string ligand_name;
            int job_id = -1;
            bool complete = false;
            double energy = 0.0;
            double delta_eV = 0.0;
            sbox::analysis::DOrbitalEnergies d_orbitals;
            sbox::chem::MolecularSystem geometry;
        };

        int metal_Z = 26;
        int oxidation_state = 2;
        sbox::chem::CoordinationGeometry geometry = sbox::chem::CoordinationGeometry::Octahedral;
        sbox::backend::Method method = sbox::backend::Method::GFN2_XTB;
        std::vector<bool> ligand_selected;
        std::vector<LigandResult> results;
        bool series_running = false;
        bool series_complete = false;
    };

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

    struct TrajectoryPlayerState {
        int current_frame = 0;
        bool playing = false;
        bool looping = true;
        float playback_speed = 1.0f;
        float accumulator = 0.0f;
        int total_frames = 0;
        int last_applied_frame = -1;
    };

    struct GeometricConstraint {
        enum class Type { FreezeAtom, FixDistance, FixAngle, FixDihedral };
        Type type = Type::FreezeAtom;
        std::vector<int> atom_indices;
        double value = 0.0;
        bool active = true;
    };

    struct PESScanState {
        sbox::backend::JobSpec::ScanSpec::ScanCoordinate coord1;
        sbox::backend::JobSpec::ScanSpec::ScanCoordinate coord2;
        bool is_2d = false;
        bool coord1_set = false;
        bool coord2_set = false;
        int scan_job_id = -1;
        bool scan_running = false;
        bool scan_complete = false;
        sbox::backend::JobResult::ScanResult result;
        int viewed_point = -1;
    };

    struct ReactionPathState {
        struct TrackedCoordinate {
            std::string label;
            sbox::backend::JobSpec::ScanSpec::CoordType type = sbox::backend::JobSpec::ScanSpec::CoordType::Distance;
            std::vector<int> atoms;
            std::vector<double> values;
        };

        sbox::chem::MolecularSystem reactant;
        sbox::chem::MolecularSystem product;
        bool reactant_set = false;
        bool product_set = false;
        int num_images = 9;
        int max_neb_steps = 50;
        int neb_job_id = -1;
        bool neb_running = false;
        bool neb_complete = false;
        sbox::backend::JobResult::NEBResult result;
        TrajectoryPlayerState player;
        std::vector<TrackedCoordinate> tracked;
        bool smooth_interpolation = true;
        int ts_frequency_job_id = -1;
        bool ts_frequency_running = false;
        sbox::backend::JobResult ts_frequency_result;
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
    bool show_nci = false;
    bool show_symmetry_elements = false;
    bool show_complex_builder = false;
    bool show_spectrochemical = false;
    bool show_settings = false;
    int color_mode = 0;
    float esp_density_iso = 0.005f;
    float esp_color_min = -0.05f;
    float esp_color_max = 0.05f;
    bool esp_auto_range = true;
    float nci_rdg_iso = 0.3f;
    float nci_rho_cutoff = 0.05f;
    float nci_color_range = 0.04f;
    bool nci_compute_requested = false;
    std::vector<float> nci_plot_rdg;
    std::vector<float> nci_plot_sign_rho;
    std::vector<SymmetryElement> cached_symmetry_elements;
    bool symmetry_elements_dirty = true;

    int render_mode = 0;  // 0=volume, 1=isosurface, 2=phase isosurface
    float iso_value = 0.01f;
    float gamma = 0.4f;
    ViewMode view_mode = ViewMode::AtomicOrbital;
    PropertyView property_view = PropertyView::Dashboard;
    int selected_mo = -1;
    int selected_d_orbital = 0;
    int selected_crystal_field_metal = 0;
    bool d_orbital_high_spin = true;
    bool d_orbital_compare_free_ion = false;
    int selected_vibrational_mode = -1;
    bool animate_vibrational_mode = false;
    int num_mo = 0;
    int homo_index = -1;
    float mol_bound_radius = 10.0f;
    int mol_render_mode = 0;
    int lod_atoms_rendered = 0;
    int lod_atoms_culled = 0;
    int lod_bonds_rendered = 0;
    int mol_num_basis = 0;
    double mol_total_energy_h = 0.0;
    double mol_homo_lumo_gap_ev = 0.0;
    bool mol_has_mo_summary = false;
    ComputationState computation;
    TrajectoryPlayerState optimization_player;
    std::vector<GeometricConstraint> constraints;
    PESScanState pes;
    ReactionPathState reaction_path;
    ComplexBuilderState complex_builder;
    SolventPanelState solvent;
    SpectrochemicalState spectrochemical;

    int current_n = 1;
    int current_l = 0;
    float current_Zeff = 1.0f;
    std::string current_rendering_mode = "Rendering: Forward";
    std::string render_stats_summary;

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
