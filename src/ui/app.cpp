#include "ui/app.h"

#include "core/hydrogen.h"
#include "ui/panels.h"

#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_internal.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace sbox {

App::App() {
    window_ = std::make_unique<Window>(1600, 1000, "Schrödinger's Sandbox");
    gradient_shader_ = std::make_unique<Shader>("data/shaders/fullscreen_quad.vert", "data/shaders/test_gradient.frag");

    try {
        orbital_shader_ = std::make_unique<Shader>("data/shaders/orbital_raymarch.vert", "data/shaders/orbital_raymarch.frag");

        int info_log_len = 0;
        glGetProgramiv(orbital_shader_->id(), GL_INFO_LOG_LENGTH, &info_log_len);
        if (info_log_len > 1) {
            std::vector<char> info_log(static_cast<std::size_t>(info_log_len));
            glGetProgramInfoLog(orbital_shader_->id(), info_log_len, nullptr, info_log.data());
            std::cerr << "Orbital shader compile/link log: " << info_log.data() << '\n';
        } else {
            std::cerr << "Orbital shader compile/link log: <empty>\n";
        }
    } catch (const std::exception& ex) {
        std::cerr << "Orbital shader unavailable (full compile/link error log emitted above), "
                     "falling back to gradient shader: "
                  << ex.what() << '\n';
        orbital_shader_.reset();
    }

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

        if (state_.needs_update) {
            state_.update();
        }

        ui::draw_periodic_table(state_);
        ui::draw_orbital_browser(state_);
        ui::draw_properties(state_);
        ui::draw_energy_diagram(state_);

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

        if (state_.needs_update) {
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

    glBindFramebuffer(GL_FRAMEBUFFER, viewport_fbo_);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, viewport_color_tex_, 0);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        throw std::runtime_error("Viewport framebuffer is incomplete");
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

float App::computeMaxDensityEstimate() const {
    const int n = state_.current_n;
    const int l = state_.current_l;
    const int m = state_.selected_m;
    const double zeff = static_cast<double>(state_.current_Zeff);

    if (n <= 0 || zeff <= 0.0) {
        return 1.0f;
    }

    double max_density = 0.0;

    if (l == 0) {
        const double r = static_cast<double>(n * n) / zeff;
        max_density = sbox::hydrogen::probability_density(n, l, m, zeff, 0.0, 0.0, r);
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
    std::fprintf(stderr, "max_density=%.8f\n", max_density);

    return static_cast<float>(max_density);
}

void App::updateMaxDensityEstimate() {
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
    glClearColor(0.04f, 0.055f, 0.09f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

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

        if (!printed_uniforms_once_) {
            const Eigen::Vector3f cam_pos = camera_.camera_position();
            std::fprintf(stderr,
                         "UNIFORMS: n=%d l=%d m=%d Zeff=%.4f mode=%d iso=%.4f gamma=%.4f max_density=%.8f\n",
                         state_.current_n,
                         state_.current_l,
                         state_.selected_m,
                         state_.current_Zeff,
                         state_.render_mode,
                         state_.iso_value,
                         state_.gamma,
                         max_density_estimate_);
            std::fprintf(stderr, "CAMERA: pos=(%.2f,%.2f,%.2f)\n", cam_pos.x(), cam_pos.y(), cam_pos.z());
            printed_uniforms_once_ = true;
        }

        for (GLenum err = glGetError(); err != GL_NO_ERROR; err = glGetError()) {
            std::cerr << "Post-uniform error: 0x" << std::hex << err << std::dec << '\n';
        }
    }

    glBindVertexArray(fullscreen_vao_);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

}  // namespace sbox
