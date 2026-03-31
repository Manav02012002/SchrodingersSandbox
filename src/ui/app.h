#pragma once

#include "analysis/nci.h"
#include "backend/backend_manager.h"
#include "backend/python_env.h"
#include "chem/ligand_library.h"
#include "core/basis_set.h"
#include "core/update_checker.h"
#include "core/molecular_system.h"
#include "core/settings.h"
#include "editor/draw_mode.h"
#include "editor/erase_mode.h"
#include "editor/fragment_mode.h"
#include "editor/measure_mode.h"
#include "editor/picking.h"
#include "editor/select_mode.h"
#include "io/cube_io.h"
#include "io/fchk_io.h"
#include "io/pdb_io.h"
#include "io/sdf_io.h"
#include "io/trajectory_io.h"
#include "io/xyz_io.h"
#include "renderer/basis_texture.h"
#include "renderer/camera.h"
#include "renderer/esp_surface.h"
#include "renderer/gbuffer.h"
#include "renderer/lod_renderer.h"
#include "renderer/mol_renderer.h"
#include "renderer/post_process.h"
#include "renderer/ssao.h"
#include "renderer/screenshot.h"
#include "renderer/shader.h"
#include "renderer/shadow_map.h"
#include "renderer/volume_texture.h"
#include "renderer/window.h"
#include "ui/annotations.h"
#include "ui/annotation_editor.h"
#include "ui/app_state.h"
#include "ui/about_dialog.h"
#include "ui/context_menu.h"
#include "ui/editor_toolbar.h"
#include "ui/nci_panel.h"
#include "ui/results_panel.h"
#include "ui/setup_wizard.h"
#include "ui/symmetry_overlay.h"

#include <memory>
#include <optional>
#include <array>
#include <cstddef>
#include <string>

struct GLFWwindow;

namespace sbox {

class App {
public:
    struct GpuTimingState {
        enum class Pass : int { GBuffer = 0, SSAO, Shadows, Lighting, Orbitals, Post, Count };
        std::array<std::array<unsigned int, 2>, static_cast<std::size_t>(Pass::Count)> queries{};
        bool initialized = false;
    };

    App();
    ~App();

    App(const App&) = delete;
    App& operator=(const App&) = delete;

    void run();
    void load_file_by_extension(const std::string& path);
    void render_single_frame();
    void render_to_fbo(unsigned int fbo, int w, int h);
    void render_to_fbo(unsigned int fbo, int w, int h, bool transparent_background);
    ui::AppState& state();
    const ui::AppState& state() const;
    Camera& camera();
    const Camera& camera() const;
    const sbox::chem::MolecularSystem& current_molecule() const;
    const sbox::io::Trajectory& current_trajectory() const;
    bool has_trajectory() const;
    bool has_mo_data() const;
    void set_current_molecule_for_export(const sbox::chem::MolecularSystem& mol);

private:
    static void ScrollCallback(GLFWwindow* window, double xoffset, double yoffset);

    void initImGui();
    void shutdownImGui();
    void rebuild_imgui_scale();
    void on_content_scale_change(float x_scale, float y_scale);
    void ensureViewportTarget(int width, int height);
    void renderViewportToTexture();
    void renderForwardToTarget(unsigned int fbo, int width, int height, bool transparent_background = false);
    void renderDeferredToTarget(unsigned int fbo, int width, int height, bool transparent_background = false);
    void renderOrbitalPass(int width, int height);
    void renderESPPass(int width, int height);
    void renderNCIPass(int width, int height);
    void renderViewportToTarget(unsigned int fbo, int width, int height, bool transparent_background = false);
    void updateMaxDensityEstimate();
    [[nodiscard]] float computeMaxDensityEstimate() const;
    void loadMoldenFile(const std::string& path);
    void loadCubeFile(const std::string& path);
    void loadXYZFile(const std::string& path);
    void loadTrajectoryFile(const std::string& path);
    void loadSDFFile(const std::string& path);
    void loadPDBFile(const std::string& path);
    void loadFchkFile(const std::string& path);
    void loadProjectFile(const std::string& path);
    void uploadCurrentMoleculeToRenderers();
    void apply_settings();
    void applyMOData(const sbox::basis::MOData& mo_data, const std::string& name_hint);
    void applyBackendResult(const sbox::backend::JobResult& result);
    void loadESPSurface(const sbox::backend::JobResult& result);
    void detect_metal_center();
    void run_crystal_field_analysis();
    [[nodiscard]] sbox::backend::JobSpec makeJobSpecFromState() const;
    [[nodiscard]] int find_homo_index() const;
    [[nodiscard]] float compute_mol_bound_radius(const sbox::chem::MolecularSystem& mol) const;

    std::unique_ptr<Window> window_;
    std::unique_ptr<Shader> gradient_shader_;
    std::unique_ptr<Shader> orbital_shader_;
    std::unique_ptr<Shader> mo_shader_;
    std::unique_ptr<Shader> cube_shader_;
    std::unique_ptr<Shader> esp_shader_;
    std::unique_ptr<Shader> nci_shader_;
    std::unique_ptr<Shader> deferred_lighting_shader_;
    Camera camera_;
    ui::AppState state_;
    ui::AnnotationManager annotation_manager_;
    ui::EditorState editor_state_;
    sbox::SettingsManager settings_manager_;
    std::unique_ptr<sbox::UpdateChecker> update_checker_;
    sbox::backend::BackendManager backend_;
    sbox::backend::PythonEnvironment python_env_;
    ui::SetupWizardState wizard_state_;
    sbox::chem::LigandLibrary ligand_library_;
    sbox::analysis::DOrbitalEnergies current_d_orbitals_;
    bool has_d_orbital_analysis_ = false;
    int current_metal_index_ = -1;
    sbox::render::BasisTextures basis_textures_;
    sbox::render::ESPSurface esp_surface_;
    sbox::render::LODRenderer lod_renderer_;
    sbox::render::MolRenderer mol_renderer_;
    sbox::render::VolumeTexture nci_rdg_texture_;
    sbox::render::VolumeTexture nci_sign_texture_;
    sbox::render::PostProcess post_process_;
    sbox::render::GBuffer gbuffer_;
    sbox::render::SSAO ssao_;
    sbox::render::ShadowMap shadow_map_;
    sbox::render::VolumeTexture volume_texture_;

    sbox::basis::MOData current_mo_data_;
    sbox::chem::MolecularSystem current_molecule_;
    sbox::io::Trajectory current_trajectory_;
    sbox::io::PDBData current_pdb_data_;
    std::optional<sbox::backend::JobResult> latest_result_;
    std::optional<sbox::analysis::NCIGrid> nci_grid_;
    bool has_trajectory_ = false;
    bool has_mo_data_ = false;
    bool has_cube_data_ = false;
    bool use_cube_fallback_ = false;

    unsigned int viewport_fbo_ = 0;
    unsigned int viewport_color_tex_ = 0;
    unsigned int viewport_depth_rbo_ = 0;
    unsigned int fullscreen_vao_ = 0;
    int viewport_width_ = 1;
    int viewport_height_ = 1;

    float scroll_delta_ = 0.0f;
    float max_density_estimate_ = 1.0f;
    bool editor_consuming_left_drag_ = false;
    bool open_message_popup_ = false;
    std::string message_popup_title_;
    std::string message_popup_text_;
    bool open_metal_picker_popup_ = false;
    std::vector<int> detected_metal_choices_;
    bool show_about_ = false;
    bool show_shortcuts_ = false;
    bool show_update_dialog_ = false;
    bool show_export_dialog_ = false;
    bool update_notification_shown_ = false;
    int update_poll_frame_counter_ = 0;
    std::optional<sbox::UpdateInfo> pending_update_;
    ImGuiStyle base_imgui_style_{};
    bool base_imgui_style_initialized_ = false;
    float imgui_effective_scale_ = 1.0f;
    GpuTimingState gpu_timing_;

    int density_n_ = -1;
    int density_l_ = -1;
    int density_m_ = -999;
    float density_zeff_ = -1.0f;
};

}  // namespace sbox
#include <imgui.h>
