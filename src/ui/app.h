#pragma once

#include "backend/backend_manager.h"
#include "backend/python_env.h"
#include "core/basis_set.h"
#include "core/molecular_system.h"
#include "editor/draw_mode.h"
#include "editor/erase_mode.h"
#include "editor/fragment_mode.h"
#include "editor/measure_mode.h"
#include "editor/picking.h"
#include "editor/select_mode.h"
#include "io/cube_io.h"
#include "io/fchk_io.h"
#include "io/sdf_io.h"
#include "io/xyz_io.h"
#include "renderer/basis_texture.h"
#include "renderer/camera.h"
#include "renderer/esp_surface.h"
#include "renderer/mol_renderer.h"
#include "renderer/shader.h"
#include "renderer/volume_texture.h"
#include "renderer/window.h"
#include "ui/app_state.h"
#include "ui/context_menu.h"
#include "ui/editor_toolbar.h"
#include "ui/results_panel.h"
#include "ui/setup_wizard.h"

#include <memory>
#include <optional>
#include <string>

struct GLFWwindow;

namespace sbox {

class App {
public:
    App();
    ~App();

    App(const App&) = delete;
    App& operator=(const App&) = delete;

    void run();

private:
    static void ScrollCallback(GLFWwindow* window, double xoffset, double yoffset);

    void initImGui();
    void shutdownImGui();
    void ensureViewportTarget(int width, int height);
    void renderViewportToTexture();
    void updateMaxDensityEstimate();
    [[nodiscard]] float computeMaxDensityEstimate() const;
    void loadMoldenFile(const std::string& path);
    void loadCubeFile(const std::string& path);
    void loadXYZFile(const std::string& path);
    void loadSDFFile(const std::string& path);
    void loadFchkFile(const std::string& path);
    void applyMOData(const sbox::basis::MOData& mo_data, const std::string& name_hint);
    void applyBackendResult(const sbox::backend::JobResult& result);
    void loadESPSurface(const sbox::backend::JobResult& result);
    [[nodiscard]] sbox::backend::JobSpec makeJobSpecFromState() const;
    [[nodiscard]] int find_homo_index() const;
    [[nodiscard]] float compute_mol_bound_radius(const sbox::chem::MolecularSystem& mol) const;

    std::unique_ptr<Window> window_;
    std::unique_ptr<Shader> gradient_shader_;
    std::unique_ptr<Shader> orbital_shader_;
    std::unique_ptr<Shader> mo_shader_;
    std::unique_ptr<Shader> cube_shader_;
    std::unique_ptr<Shader> esp_shader_;
    Camera camera_;
    ui::AppState state_;
    ui::EditorState editor_state_;
    sbox::backend::BackendManager backend_;
    sbox::backend::PythonEnvironment python_env_;
    ui::SetupWizardState wizard_state_;
    sbox::render::BasisTextures basis_textures_;
    sbox::render::ESPSurface esp_surface_;
    sbox::render::MolRenderer mol_renderer_;
    sbox::render::VolumeTexture volume_texture_;

    sbox::basis::MOData current_mo_data_;
    sbox::chem::MolecularSystem current_molecule_;
    std::optional<sbox::backend::JobResult> latest_result_;
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

    int density_n_ = -1;
    int density_l_ = -1;
    int density_m_ = -999;
    float density_zeff_ = -1.0f;
};

}  // namespace sbox
