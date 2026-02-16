#pragma once

#include "renderer/camera.h"
#include "renderer/shader.h"
#include "renderer/window.h"
#include "ui/app_state.h"

#include <memory>

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

    std::unique_ptr<Window> window_;
    std::unique_ptr<Shader> gradient_shader_;
    std::unique_ptr<Shader> orbital_shader_;
    Camera camera_;
    ui::AppState state_;

    unsigned int viewport_fbo_ = 0;
    unsigned int viewport_color_tex_ = 0;
    unsigned int fullscreen_vao_ = 0;
    int viewport_width_ = 1;
    int viewport_height_ = 1;

    float scroll_delta_ = 0.0f;
    float max_density_estimate_ = 1.0f;

    int density_n_ = -1;
    int density_l_ = -1;
    int density_m_ = -999;
    float density_zeff_ = -1.0f;
    bool printed_uniforms_once_ = false;
};

}  // namespace sbox
