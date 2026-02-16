#include "ui/app.h"

#include "core/hydrogen.h"
#include "ui/panels.h"

#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>

#include <algorithm>
#include <cstdio>
#include <cmath>
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
    while (!window_->shouldClose()) {
        window_->pollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::DockSpaceOverViewport();

        if (state_.needs_update) {
            state_.update();
        }

        ui::draw_periodic_table(state_);
        ui::draw_orbital_browser(state_);
        ui::draw_properties(state_);

        const ui::ViewportPanelState viewport_state = ui::draw_viewport(state_, viewport_color_tex_);
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

    ImGui::StyleColorsDark();

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
