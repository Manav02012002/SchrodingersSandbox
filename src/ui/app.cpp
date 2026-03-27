#include "ui/app.h"

#include "core/gaussian_eval.h"
#include "core/hydrogen.h"
#include "core/molden_parser.h"
#include "ui/file_dialog.h"
#include "ui/mo_diagram.h"
#include "ui/molecule_info.h"
#include "ui/panels.h"

#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_internal.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

namespace sbox {

namespace {

std::unique_ptr<Shader> try_load_shader(const std::string& vertex,
                                        const std::string& fragment,
                                        const char* label) {
    try {
        auto shader = std::make_unique<Shader>(vertex, fragment);
        int info_log_len = 0;
        glGetProgramiv(shader->id(), GL_INFO_LOG_LENGTH, &info_log_len);
        if (info_log_len > 1) {
            std::vector<char> info_log(static_cast<std::size_t>(info_log_len));
            glGetProgramInfoLog(shader->id(), info_log_len, nullptr, info_log.data());
            std::cerr << label << " shader compile/link log: " << info_log.data() << '\n';
        } else {
            std::cerr << label << " shader compile/link log: <empty>\n";
        }
        return shader;
    } catch (const std::exception& ex) {
        std::cerr << label << " shader unavailable: " << ex.what() << '\n';
        return nullptr;
    }
}

void update_bound_radius_from_molecule(ui::AppState& state, const sbox::chem::MolecularSystem& molecule) {
    if (molecule.num_atoms() == 0) {
        state.mol_bound_radius = 10.0f;
        return;
    }

    const Eigen::Vector3d center = molecule.center_of_mass();
    double max_dist = 0.0;
    for (const sbox::chem::Atom& atom : molecule.atoms()) {
        max_dist = std::max(max_dist, (atom.position - center).norm());
    }
    state.mol_bound_radius = static_cast<float>(max_dist + 10.0);
}

void clear_mo_summary(ui::AppState& state) {
    state.num_mo = 0;
    state.homo_index = -1;
    state.selected_mo = -1;
    state.mol_num_basis = 0;
    state.mol_total_energy_h = 0.0;
    state.mol_homo_lumo_gap_ev = 0.0;
    state.mol_has_mo_summary = false;
}

std::optional<sbox::basis::MOData> mo_data_from_fchk(const sbox::io::FchkData& fchk) {
    if (fchk.shell_types.empty() || fchk.shell_to_atom_map.empty() || fchk.primitives_per_shell.empty() ||
        fchk.primitive_exponents.empty() || fchk.contraction_coefficients.empty() || fchk.mo_coefficients.size() == 0) {
        return std::nullopt;
    }

    if (fchk.shell_types.size() != fchk.shell_to_atom_map.size() ||
        fchk.shell_types.size() != fchk.primitives_per_shell.size()) {
        return std::nullopt;
    }

    sbox::basis::MOData mo_data;
    mo_data.atomic_numbers = fchk.atomic_numbers;
    mo_data.energies = fchk.mo_energies;
    mo_data.occupations = fchk.occupations;
    mo_data.coefficients = fchk.mo_coefficients;
    mo_data.total_energy = fchk.total_energy;
    mo_data.atom_positions.reserve(static_cast<std::size_t>(fchk.num_atoms));
    for (int i = 0; i < fchk.num_atoms; ++i) {
        mo_data.atom_positions.emplace_back(fchk.coordinates[static_cast<std::size_t>(3 * i + 0)],
                                            fchk.coordinates[static_cast<std::size_t>(3 * i + 1)],
                                            fchk.coordinates[static_cast<std::size_t>(3 * i + 2)]);
    }

    bool spherical = true;
    for (int shell_type : fchk.shell_types) {
        if (shell_type > 1) {
            spherical = false;
        }
        if (shell_type == -1) {
            continue;
        }
    }
    mo_data.basis.spherical = spherical;

    std::size_t primitive_offset = 0;
    for (std::size_t shell_index = 0; shell_index < fchk.shell_types.size(); ++shell_index) {
        const int shell_type = fchk.shell_types[shell_index];
        const int atom_index = fchk.shell_to_atom_map[shell_index] - 1;
        const int nprim = fchk.primitives_per_shell[shell_index];
        if (atom_index < 0 || atom_index >= fchk.num_atoms || nprim < 0) {
            return std::nullopt;
        }
        if (primitive_offset + static_cast<std::size_t>(nprim) > fchk.primitive_exponents.size() ||
            primitive_offset + static_cast<std::size_t>(nprim) > fchk.contraction_coefficients.size()) {
            return std::nullopt;
        }

        auto fill_shell = [&](sbox::basis::BasisShell& shell, const std::vector<double>& coeffs) {
            shell.atom_index = atom_index;
            shell.primitives.reserve(static_cast<std::size_t>(nprim));
            for (int p = 0; p < nprim; ++p) {
                shell.primitives.push_back({fchk.primitive_exponents[primitive_offset + static_cast<std::size_t>(p)],
                                            coeffs[primitive_offset + static_cast<std::size_t>(p)]});
            }
        };

        if (shell_type == -1) {
            if (primitive_offset + static_cast<std::size_t>(nprim) > fchk.sp_contraction_coefficients.size()) {
                return std::nullopt;
            }
            sbox::basis::BasisShell shell_s;
            shell_s.angular_momentum = 0;
            fill_shell(shell_s, fchk.contraction_coefficients);
            mo_data.basis.shells.push_back(shell_s);

            sbox::basis::BasisShell shell_p;
            shell_p.angular_momentum = 1;
            fill_shell(shell_p, fchk.sp_contraction_coefficients);
            mo_data.basis.shells.push_back(shell_p);
        } else {
            sbox::basis::BasisShell shell;
            shell.angular_momentum = std::abs(shell_type);
            if (shell.angular_momentum < 0 || shell.angular_momentum > 3) {
                return std::nullopt;
            }
            fill_shell(shell, fchk.contraction_coefficients);
            mo_data.basis.shells.push_back(shell);
        }

        primitive_offset += static_cast<std::size_t>(nprim);
    }

    if (mo_data.basis.num_basis_functions() != mo_data.coefficients.rows()) {
        return std::nullopt;
    }

    return mo_data;
}

}  // namespace

App::App() {
    window_ = std::make_unique<Window>(1600, 1000, "Schrödinger's Sandbox");
    gradient_shader_ = std::make_unique<Shader>("data/shaders/fullscreen_quad.vert", "data/shaders/test_gradient.frag");
    orbital_shader_ = try_load_shader("data/shaders/orbital_raymarch.vert", "data/shaders/orbital_raymarch.frag", "Orbital");
    mo_shader_ = try_load_shader("data/shaders/mo_raymarch.vert", "data/shaders/mo_raymarch.frag", "MO");
    cube_shader_ = try_load_shader("data/shaders/cube_raymarch.vert", "data/shaders/cube_raymarch.frag", "Cube");

    glGenVertexArrays(1, &fullscreen_vao_);

    glfwSetWindowUserPointer(window_->handle(), this);
    glfwSetScrollCallback(window_->handle(), &App::ScrollCallback);

    initImGui();
    ensureViewportTarget(viewport_width_, viewport_height_);
}

App::~App() {
    shutdownImGui();

    if (fullscreen_vao_ != 0U) {
        glDeleteVertexArrays(1, &fullscreen_vao_);
        fullscreen_vao_ = 0;
    }

    if (viewport_depth_rbo_ != 0U) {
        glDeleteRenderbuffers(1, &viewport_depth_rbo_);
        viewport_depth_rbo_ = 0;
    }

    if (viewport_color_tex_ != 0U) {
        glDeleteTextures(1, &viewport_color_tex_);
        viewport_color_tex_ = 0;
    }

    if (viewport_fbo_ != 0U) {
        glDeleteFramebuffers(1, &viewport_fbo_);
        viewport_fbo_ = 0;
    }
}

void App::run() {
    static bool dock_layout_built = false;

    while (!window_->shouldClose()) {
        window_->pollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        if (ImGui::BeginMainMenuBar()) {
            if (ImGui::BeginMenu("File")) {
                if (ImGui::MenuItem("Open Molden...")) {
                    try {
                        const std::string path = ui::open_file_dialog("Open Molden File", "molden");
                        if (!path.empty()) {
                            loadMoldenFile(path);
                        }
                    } catch (const std::exception& ex) {
                        std::cerr << "Failed to load Molden file: " << ex.what() << '\n';
                    }
                }
                if (ImGui::MenuItem("Open XYZ...")) {
                    try {
                        const std::string path = ui::open_file_dialog("Open XYZ File", "xyz");
                        if (!path.empty()) {
                            loadXYZFile(path);
                        }
                    } catch (const std::exception& ex) {
                        std::cerr << "Failed to load XYZ file: " << ex.what() << '\n';
                    }
                }
                if (ImGui::MenuItem("Open SDF...")) {
                    try {
                        const std::string path = ui::open_file_dialog("Open SDF File", "sdf,mol");
                        if (!path.empty()) {
                            loadSDFFile(path);
                        }
                    } catch (const std::exception& ex) {
                        std::cerr << "Failed to load SDF file: " << ex.what() << '\n';
                    }
                }
                if (ImGui::MenuItem("Open Cube...")) {
                    try {
                        const std::string path = ui::open_file_dialog("Open Cube File", "cube");
                        if (!path.empty()) {
                            loadCubeFile(path);
                        }
                    } catch (const std::exception& ex) {
                        std::cerr << "Failed to load Cube file: " << ex.what() << '\n';
                    }
                }
                if (ImGui::MenuItem("Open FCHK...")) {
                    try {
                        const std::string path = ui::open_file_dialog("Open FCHK File", "fchk,fch");
                        if (!path.empty()) {
                            loadFchkFile(path);
                        }
                    } catch (const std::exception& ex) {
                        std::cerr << "Failed to load FCHK file: " << ex.what() << '\n';
                    }
                }
                ImGui::Separator();
                if (ImGui::MenuItem("Save Project...")) {
                    (void)ui::open_file_dialog("Save Project", "sbox");
                }
                ImGui::Separator();
                if (ImGui::MenuItem("Quit")) {
                    glfwSetWindowShouldClose(window_->handle(), GLFW_TRUE);
                }
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("View")) {
                if (ImGui::MenuItem("Reset Layout")) {
                    dock_layout_built = false;
                    const ImGuiID reset_dockspace_id = ImGui::GetID("MainDockSpace");
                    ImGui::DockBuilderRemoveNode(reset_dockspace_id);
                }
                ImGui::EndMenu();
            }
            ImGui::EndMainMenuBar();
        }

        const ImGuiID dockspace_id = ImGui::GetID("MainDockSpace");
        ImGui::DockSpaceOverViewport(dockspace_id, ImGui::GetMainViewport());

        if (!dock_layout_built) {
            ImGuiDockNode* dockspace_node = ImGui::DockBuilderGetNode(dockspace_id);
            const bool has_existing_layout =
                dockspace_node != nullptr &&
                (dockspace_node->ChildNodes[0] != nullptr ||
                 dockspace_node->ChildNodes[1] != nullptr ||
                 dockspace_node->Windows.Size > 0);

            if (!has_existing_layout) {
                ImGui::DockBuilderRemoveNode(dockspace_id);
                ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
                ImGui::DockBuilderSetNodeSize(dockspace_id, ImGui::GetMainViewport()->WorkSize);

                ImGuiID left_node = dockspace_id;
                ImGuiID center_right_node = 0;
                ImGui::DockBuilderSplitNode(left_node, ImGuiDir_Left, 0.20f, &left_node, &center_right_node);

                ImGuiID top_center_right_node = center_right_node;
                ImGuiID bottom_center_right_node = 0;
                ImGui::DockBuilderSplitNode(top_center_right_node,
                                            ImGuiDir_Down,
                                            0.30f,
                                            &bottom_center_right_node,
                                            &top_center_right_node);

                ImGuiID center_node = top_center_right_node;
                ImGuiID right_node = 0;
                ImGui::DockBuilderSplitNode(center_node, ImGuiDir_Left, 0.6875f, &center_node, &right_node);

                ImGuiID orbital_node = left_node;
                ImGuiID properties_node = 0;
                ImGui::DockBuilderSplitNode(orbital_node, ImGuiDir_Down, 0.45f, &properties_node, &orbital_node);

                ImGui::DockBuilderDockWindow("Orbital Browser", orbital_node);
                ImGui::DockBuilderDockWindow("Properties", properties_node);
                ImGui::DockBuilderDockWindow("Energy Levels", properties_node);
                ImGui::DockBuilderDockWindow("3D Viewport", center_node);
                ImGui::DockBuilderDockWindow("Periodic Table", bottom_center_right_node);
                ImGui::DockBuilderFinish(dockspace_id);
                (void)right_node;
            }

            dock_layout_built = true;
        }

        if (state_.view_mode == ui::ViewMode::AtomicOrbital && state_.needs_update) {
            state_.update();
        }

        ui::draw_periodic_table(state_);
        ui::draw_orbital_browser(state_);
        if (state_.view_mode == ui::ViewMode::MolecularOrbital) {
            ui::draw_molecule_info(state_, current_molecule_);
        } else {
            ui::draw_properties(state_);
        }
        if (state_.view_mode == ui::ViewMode::MolecularOrbital && has_mo_data_) {
            ui::draw_mo_diagram(state_, current_mo_data_);
        } else {
            ui::draw_energy_diagram(state_);
        }

        const ui::ViewportPanelState viewport_state = ui::draw_viewport(state_, viewport_color_tex_);
        ui::draw_status_bar(state_);
        const int new_width = std::max(1, static_cast<int>(std::floor(viewport_state.size.x)));
        const int new_height = std::max(1, static_cast<int>(std::floor(viewport_state.size.y)));
        ensureViewportTarget(new_width, new_height);

        camera_.setViewportSize(static_cast<float>(viewport_width_), static_cast<float>(viewport_height_));

        const bool left_drag = glfwGetMouseButton(window_->handle(), GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        const bool middle_drag = glfwGetMouseButton(window_->handle(), GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS;

        const ImVec2 mouse_delta = ImGui::GetIO().MouseDelta;
        camera_.handleRotationDrag(mouse_delta.x, mouse_delta.y, viewport_state.hovered && left_drag);
        camera_.handlePanDrag(mouse_delta.x, mouse_delta.y, viewport_state.hovered && middle_drag);
        camera_.handleScroll(scroll_delta_, viewport_state.hovered);
        scroll_delta_ = 0.0f;

        if (state_.view_mode == ui::ViewMode::AtomicOrbital && state_.needs_update) {
            state_.update();
        }
        updateMaxDensityEstimate();

        renderViewportToTexture();

        ImGui::Render();

        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0, 0, window_->framebufferWidth(), window_->framebufferHeight());
        glClearColor(0.07f, 0.08f, 0.10f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        window_->swapBuffers();
    }
}

void App::ScrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    (void)xoffset;
    auto* self = static_cast<App*>(glfwGetWindowUserPointer(window));
    if (self != nullptr) {
        self->scroll_delta_ += static_cast<float>(yoffset);
    }
}

void App::initImGui() {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;

    ImGuiStyle& style = ImGui::GetStyle();
    ImGui::StyleColorsDark(&style);
    style.WindowRounding = 4.0f;
    style.FrameRounding = 3.0f;
    style.GrabRounding = 3.0f;
    style.FramePadding = ImVec2(6.0f, 4.0f);
    style.ChildBorderSize = 1.0f;
    style.ScrollbarSize = 12.0f;
    style.ScrollbarRounding = 6.0f;

    ImVec4* colors = style.Colors;
    const ImVec4 accent_active(0.15f, 0.55f, 0.65f, 1.0f);
    const ImVec4 accent(0.13f, 0.44f, 0.52f, 1.0f);
    const ImVec4 accent_hovered(0.18f, 0.62f, 0.74f, 1.0f);

    colors[ImGuiCol_WindowBg] = ImVec4(0.10f, 0.10f, 0.14f, 1.0f);
    colors[ImGuiCol_TitleBg] = ImVec4(0.08f, 0.08f, 0.12f, 1.0f);
    colors[ImGuiCol_TitleBgActive] = ImVec4(0.12f, 0.12f, 0.18f, 1.0f);
    colors[ImGuiCol_FrameBgActive] = accent_active;
    colors[ImGuiCol_ButtonActive] = accent_active;
    colors[ImGuiCol_HeaderActive] = accent_active;
    colors[ImGuiCol_Button] = accent;
    colors[ImGuiCol_ButtonHovered] = accent_hovered;
    colors[ImGuiCol_Header] = accent;
    colors[ImGuiCol_HeaderHovered] = accent_hovered;
    colors[ImGuiCol_CheckMark] = accent_hovered;
    colors[ImGuiCol_SliderGrab] = accent;
    colors[ImGuiCol_SliderGrabActive] = accent_active;
    colors[ImGuiCol_ResizeGrip] = accent;
    colors[ImGuiCol_ResizeGripHovered] = accent_hovered;
    colors[ImGuiCol_ResizeGripActive] = accent_active;
    colors[ImGuiCol_Tab] = ImVec4(0.10f, 0.26f, 0.30f, 1.0f);
    colors[ImGuiCol_TabHovered] = accent_hovered;
    colors[ImGuiCol_TabActive] = accent_active;
    colors[ImGuiCol_TabUnfocused] = ImVec4(0.08f, 0.18f, 0.22f, 1.0f);
    colors[ImGuiCol_TabUnfocusedActive] = ImVec4(0.12f, 0.38f, 0.45f, 1.0f);

    io.Fonts->Clear();
    constexpr float kFontSize = 15.0f;
    const std::array<const char*, 8> font_candidates = {
        "/Library/Fonts/Inter-Regular.ttf",
        "/System/Library/Fonts/Supplemental/Inter.ttc",
        "/System/Library/Fonts/Inter.ttc",
        "/usr/share/fonts/truetype/inter/Inter-Regular.ttf",
        "/usr/share/fonts/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/local/share/fonts/Inter-Regular.ttf",
        "/usr/local/share/fonts/dejavu/DejaVuSans.ttf",
    };

    ImFont* ui_font = nullptr;
    for (const char* path : font_candidates) {
        if (std::filesystem::exists(path)) {
            ui_font = io.Fonts->AddFontFromFileTTF(path, kFontSize);
            if (ui_font != nullptr) {
                break;
            }
        }
    }

    if (ui_font == nullptr) {
        ImFontConfig font_cfg;
        font_cfg.SizePixels = kFontSize;
        ui_font = io.Fonts->AddFontDefault(&font_cfg);
    }
    io.FontDefault = ui_font;

    if (!ImGui_ImplGlfw_InitForOpenGL(window_->handle(), true)) {
        throw std::runtime_error("Failed to initialize ImGui GLFW backend");
    }
    if (!ImGui_ImplOpenGL3_Init("#version 410")) {
        throw std::runtime_error("Failed to initialize ImGui OpenGL3 backend");
    }
}

void App::shutdownImGui() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void App::ensureViewportTarget(int width, int height) {
    width = std::max(width, 1);
    height = std::max(height, 1);

    if (viewport_fbo_ != 0U && width == viewport_width_ && height == viewport_height_) {
        return;
    }

    viewport_width_ = width;
    viewport_height_ = height;

    if (viewport_fbo_ == 0U) {
        glGenFramebuffers(1, &viewport_fbo_);
    }
    if (viewport_color_tex_ == 0U) {
        glGenTextures(1, &viewport_color_tex_);
    }
    if (viewport_depth_rbo_ == 0U) {
        glGenRenderbuffers(1, &viewport_depth_rbo_);
    }

    glBindTexture(GL_TEXTURE_2D, viewport_color_tex_);
    glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_RGBA8,
                 viewport_width_,
                 viewport_height_,
                 0,
                 GL_RGBA,
                 GL_UNSIGNED_BYTE,
                 nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindRenderbuffer(GL_RENDERBUFFER, viewport_depth_rbo_);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, viewport_width_, viewport_height_);

    glBindFramebuffer(GL_FRAMEBUFFER, viewport_fbo_);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, viewport_color_tex_, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, viewport_depth_rbo_);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        throw std::runtime_error("Viewport framebuffer is incomplete");
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

int App::find_homo_index() const {
    int homo = -1;
    for (int i = 0; i < current_mo_data_.occupations.size(); ++i) {
        if (current_mo_data_.occupations(i) > 0.5) {
            homo = i;
        }
    }
    if (homo < 0 && current_mo_data_.coefficients.cols() > 0) {
        homo = 0;
    }
    return homo;
}

void App::loadMoldenFile(const std::string& path) {
    current_mo_data_ = sbox::molden::parse_molden_file(path);
    current_molecule_.clear();
    current_molecule_.set_name(std::filesystem::path(path).filename().string());

    const std::size_t atom_count = std::min(current_mo_data_.atom_positions.size(), current_mo_data_.atomic_numbers.size());
    for (std::size_t i = 0; i < atom_count; ++i) {
        current_molecule_.add_atom({current_mo_data_.atomic_numbers[i], current_mo_data_.atom_positions[i], "", 0});
    }
    current_molecule_.perceive_bonds();
    mol_renderer_.upload(current_molecule_);

    use_cube_fallback_ = !basis_textures_.upload(current_mo_data_);
    has_mo_data_ = !use_cube_fallback_;
    has_cube_data_ = false;
    state_.mol_has_mo_summary = true;
    state_.mol_num_basis = current_mo_data_.basis.num_basis_functions();
    state_.mol_total_energy_h = current_mo_data_.total_energy;

    update_bound_radius_from_molecule(state_, current_molecule_);

    state_.view_mode = ui::ViewMode::MolecularOrbital;
    state_.num_mo = static_cast<int>(current_mo_data_.coefficients.cols());
    state_.homo_index = find_homo_index();
    state_.selected_mo = state_.homo_index;
    if (state_.homo_index >= 0 && state_.homo_index + 1 < current_mo_data_.energies.size()) {
        state_.mol_homo_lumo_gap_ev =
            (current_mo_data_.energies(state_.homo_index + 1) - current_mo_data_.energies(state_.homo_index)) * 27.2114;
    } else {
        state_.mol_homo_lumo_gap_ev = 0.0;
    }
}

void App::loadCubeFile(const std::string& path) {
    const sbox::io::CubeData cube = sbox::io::read_cube(path);
    current_molecule_.clear();
    current_molecule_.set_name(std::filesystem::path(path).filename().string());
    for (std::size_t i = 0; i < cube.atom_Z.size() && i < cube.atom_pos.size(); ++i) {
        current_molecule_.add_atom({cube.atom_Z[i], cube.atom_pos[i], "", 0});
    }
    current_molecule_.perceive_bonds();
    mol_renderer_.upload(current_molecule_);

    if (!volume_texture_.upload(cube)) {
        throw std::runtime_error("Failed to upload cube volume texture");
    }

    has_cube_data_ = true;
    use_cube_fallback_ = true;
    has_mo_data_ = false;
    state_.view_mode = ui::ViewMode::MolecularOrbital;
    update_bound_radius_from_molecule(state_, current_molecule_);
    clear_mo_summary(state_);
}

void App::loadXYZFile(const std::string& path) {
    current_molecule_ = sbox::io::read_xyz(path);
    mol_renderer_.upload(current_molecule_);
    current_mo_data_ = sbox::basis::MOData{};
    has_mo_data_ = false;
    has_cube_data_ = false;
    use_cube_fallback_ = false;
    state_.view_mode = ui::ViewMode::MolecularOrbital;
    update_bound_radius_from_molecule(state_, current_molecule_);
    clear_mo_summary(state_);
}

void App::loadSDFFile(const std::string& path) {
    current_molecule_ = sbox::io::read_sdf(path);
    mol_renderer_.upload(current_molecule_);
    current_mo_data_ = sbox::basis::MOData{};
    has_mo_data_ = false;
    has_cube_data_ = false;
    use_cube_fallback_ = false;
    state_.view_mode = ui::ViewMode::MolecularOrbital;
    update_bound_radius_from_molecule(state_, current_molecule_);
    clear_mo_summary(state_);
}

void App::loadFchkFile(const std::string& path) {
    const sbox::io::FchkData fchk = sbox::io::read_fchk(path);
    current_molecule_.clear();
    current_molecule_.set_name(fchk.title);
    current_molecule_.set_charge(fchk.charge);
    current_molecule_.set_multiplicity(fchk.multiplicity);
    for (int i = 0; i < fchk.num_atoms; ++i) {
        current_molecule_.add_atom({fchk.atomic_numbers[static_cast<std::size_t>(i)],
                                    Eigen::Vector3d(fchk.coordinates[static_cast<std::size_t>(3 * i + 0)],
                                                    fchk.coordinates[static_cast<std::size_t>(3 * i + 1)],
                                                    fchk.coordinates[static_cast<std::size_t>(3 * i + 2)]),
                                    "",
                                    0});
    }
    current_molecule_.perceive_bonds();
    mol_renderer_.upload(current_molecule_);

    current_mo_data_ = sbox::basis::MOData{};
    const std::optional<sbox::basis::MOData> maybe_mo = mo_data_from_fchk(fchk);
    if (maybe_mo.has_value()) {
        current_mo_data_ = *maybe_mo;
        use_cube_fallback_ = !basis_textures_.upload(current_mo_data_);
        has_mo_data_ = !use_cube_fallback_;
    } else {
        has_mo_data_ = false;
        use_cube_fallback_ = false;
    }
    has_cube_data_ = false;
    state_.view_mode = ui::ViewMode::MolecularOrbital;
    update_bound_radius_from_molecule(state_, current_molecule_);
    state_.num_mo = static_cast<int>(fchk.mo_energies.size());
    state_.homo_index = -1;
    for (int i = 0; i < fchk.occupations.size(); ++i) {
        if (fchk.occupations(i) > 0.5) {
            state_.homo_index = i;
        }
    }
    state_.selected_mo = state_.homo_index;
    state_.mol_num_basis = fchk.num_basis;
    state_.mol_total_energy_h = fchk.total_energy;
    state_.mol_has_mo_summary = true;
    if (state_.homo_index >= 0 && state_.homo_index + 1 < fchk.mo_energies.size()) {
        state_.mol_homo_lumo_gap_ev =
            (fchk.mo_energies(state_.homo_index + 1) - fchk.mo_energies(state_.homo_index)) * 27.2114;
    } else {
        state_.mol_homo_lumo_gap_ev = 0.0;
    }
}

float App::computeMaxDensityEstimate() const {
    if (state_.view_mode == ui::ViewMode::MolecularOrbital) {
        if (use_cube_fallback_ && has_cube_data_ && volume_texture_.is_uploaded()) {
            const float max_psi = volume_texture_.max_abs_value();
            return std::max(max_psi * max_psi, 1e-6f);
        }

        if (has_mo_data_ && current_mo_data_.coefficients.cols() > 0) {
            const int mo_index = state_.selected_mo >= 0 ? state_.selected_mo : find_homo_index();
            double max_density = 0.0;
            Eigen::Vector3d center = current_molecule_.num_atoms() > 0 ? current_molecule_.center_of_mass() : Eigen::Vector3d::Zero();
            const std::array<Eigen::Vector3d, 7> offsets = {
                Eigen::Vector3d::Zero(),
                Eigen::Vector3d(1.0, 0.0, 0.0),
                Eigen::Vector3d(-1.0, 0.0, 0.0),
                Eigen::Vector3d(0.0, 1.0, 0.0),
                Eigen::Vector3d(0.0, -1.0, 0.0),
                Eigen::Vector3d(0.0, 0.0, 1.0),
                Eigen::Vector3d(0.0, 0.0, -1.0),
            };

            max_density = std::max(max_density,
                                   sbox::basis::evaluate_mo_density_at_point(current_mo_data_, mo_index, center));
            for (const sbox::chem::Atom& atom : current_molecule_.atoms()) {
                for (const Eigen::Vector3d& offset : offsets) {
                    max_density = std::max(max_density,
                                           sbox::basis::evaluate_mo_density_at_point(current_mo_data_,
                                                                                     mo_index,
                                                                                     atom.position + offset));
                }
            }

            if (max_density <= 0.0 || std::isnan(max_density)) {
                return 1.0f;
            }
            return static_cast<float>(std::max(max_density, 1e-6));
        }

        return 1.0f;
    }

    const int n = state_.current_n;
    const int l = state_.current_l;
    const int m = state_.selected_m;
    const double zeff = static_cast<double>(state_.current_Zeff);

    if (n <= 0 || zeff <= 0.0) {
        return 1.0f;
    }

    double max_density = 0.0;
    if (l == 0) {
        max_density = sbox::hydrogen::probability_density(n, l, m, zeff, 0.0, 0.0, 0.0);
    } else {
        const int samples = 100;
        const double max_r = 4.0 * static_cast<double>(n * n) / zeff;
        for (int i = 0; i < samples; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(samples - 1);
            const double r = max_r * t;
            const double density_z = sbox::hydrogen::probability_density(n, l, m, zeff, 0.0, 0.0, r);
            const double density_x = sbox::hydrogen::probability_density(n, l, m, zeff, r, 0.0, 0.0);
            max_density = std::max(max_density, std::max(density_z, density_x));
        }
    }

    if (max_density <= 0.0 || std::isnan(max_density)) {
        max_density = 1.0;
    } else if (max_density <= 1e-9) {
        max_density = 1e-6;
    }
    return static_cast<float>(max_density);
}

void App::updateMaxDensityEstimate() {
    if (state_.view_mode == ui::ViewMode::MolecularOrbital) {
        max_density_estimate_ = computeMaxDensityEstimate();
        return;
    }

    const bool changed = (density_n_ != state_.current_n) ||
                         (density_l_ != state_.current_l) ||
                         (density_m_ != state_.selected_m) ||
                         (std::abs(density_zeff_ - state_.current_Zeff) > 1e-5f);

    if (!changed) {
        return;
    }

    max_density_estimate_ = computeMaxDensityEstimate();
    density_n_ = state_.current_n;
    density_l_ = state_.current_l;
    density_m_ = state_.selected_m;
    density_zeff_ = state_.current_Zeff;
}

void App::renderViewportToTexture() {
    glBindFramebuffer(GL_FRAMEBUFFER, viewport_fbo_);
    glViewport(0, 0, viewport_width_, viewport_height_);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.04f, 0.055f, 0.09f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (state_.view_mode == ui::ViewMode::AtomicOrbital) {
        Shader* active_shader = orbital_shader_ ? orbital_shader_.get() : gradient_shader_.get();
        if (active_shader == nullptr) {
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            return;
        }

        active_shader->bind();
        if (orbital_shader_) {
            active_shader->setUniform("u_n", state_.current_n);
            active_shader->setUniform("u_l", state_.current_l);
            active_shader->setUniform("u_m", state_.selected_m);
            active_shader->setUniform("u_Zeff", state_.current_Zeff);
            active_shader->setUniform("u_inv_vp", camera_.inv_view_projection());
            active_shader->setUniform("u_camera_pos", camera_.camera_position());
            active_shader->setUniform("u_iso_value", state_.iso_value);
            active_shader->setUniform("u_render_mode", state_.render_mode);
            active_shader->setUniform("u_gamma", state_.gamma);

            const int res_loc = glGetUniformLocation(active_shader->id(), "u_resolution");
            if (res_loc >= 0) {
                glUniform2f(res_loc, static_cast<float>(viewport_width_), static_cast<float>(viewport_height_));
            }

            const int density_loc = glGetUniformLocation(active_shader->id(), "u_max_density");
            if (density_loc >= 0) {
                glUniform1f(density_loc, max_density_estimate_);
            }
        }

        glBindVertexArray(fullscreen_vao_);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glBindVertexArray(0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        return;
    }

    if (mol_renderer_.has_data()) {
        mol_renderer_.render(camera_.viewMatrix(),
                             camera_.projectionMatrix(),
                             camera_.camera_position(),
                             static_cast<sbox::render::MolRenderMode>(state_.mol_render_mode));
    }

    Shader* active = nullptr;
    if (use_cube_fallback_ && has_cube_data_ && cube_shader_) {
        active = cube_shader_.get();
    } else if (has_mo_data_ && mo_shader_) {
        active = mo_shader_.get();
    }

    if (active != nullptr) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        active->bind();
        active->setUniform("u_inv_vp", camera_.inv_view_projection());
        active->setUniform("u_camera_pos", camera_.camera_position());
        active->setUniform("u_render_mode", state_.render_mode);
        active->setUniform("u_iso_value", state_.iso_value);
        active->setUniform("u_gamma", state_.gamma);
        active->setUniform("u_max_density", max_density_estimate_);
        active->setUniform("u_bound_radius", state_.mol_bound_radius);

        const int res_loc = glGetUniformLocation(active->id(), "u_resolution");
        if (res_loc >= 0) {
            glUniform2f(res_loc, static_cast<float>(viewport_width_), static_cast<float>(viewport_height_));
        }

        if (!use_cube_fallback_) {
            const int mo_idx = state_.selected_mo >= 0 ? state_.selected_mo : find_homo_index();
            active->setUniform("u_mo_index", mo_idx);
            active->setUniform("u_num_shells", basis_textures_.num_shells());
            active->setUniform("u_num_basis", basis_textures_.num_basis());
            active->setUniform("u_num_mo", basis_textures_.num_mo());
            basis_textures_.bind(active->id(), 1);
        } else {
            volume_texture_.bind(active->id(), 5);
            const Eigen::Vector3f orig = volume_texture_.origin();
            const Eigen::Matrix3f w2g = volume_texture_.world_to_grid();
            const int origin_loc = glGetUniformLocation(active->id(), "u_grid_origin");
            if (origin_loc >= 0) {
                glUniform3f(origin_loc, orig.x(), orig.y(), orig.z());
            }
            const int w2g_loc = glGetUniformLocation(active->id(), "u_world_to_grid");
            if (w2g_loc >= 0) {
                glUniformMatrix3fv(w2g_loc, 1, GL_FALSE, w2g.data());
            }
            const int dims_loc = glGetUniformLocation(active->id(), "u_grid_dims");
            if (dims_loc >= 0) {
                glUniform3i(dims_loc, volume_texture_.nx(), volume_texture_.ny(), volume_texture_.nz());
            }
        }

        glBindVertexArray(fullscreen_vao_);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glBindVertexArray(0);

        if (!use_cube_fallback_) {
            basis_textures_.unbind();
        } else {
            volume_texture_.unbind();
        }

        glDisable(GL_BLEND);
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

}  // namespace sbox
