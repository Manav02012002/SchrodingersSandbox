// HiDPI manual test checklist:
// 1. On Retina Mac: verify text is sharp, not blurry
// 2. On Retina Mac: verify orbital rendering is sharp (not pixelated)
// 3. On Retina Mac: verify mouse picking works (click on atoms hits correctly)
// 4. On Retina Mac: verify screenshots at viewport resolution produce 2x pixel images
// 5. On non-Retina display: verify everything still works at scale=1.0
// 6. Moving window between Retina and non-Retina: verify resize and re-render

#include "renderer/window.h"
#include "core/logging.h"

#include <glad/gl.h>
#include <GLFW/glfw3.h>

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <cmath>

namespace sbox {

namespace {

std::unordered_map<GLFWwindow*, Window*>& WindowMap() {
    static std::unordered_map<GLFWwindow*, Window*> map;
    return map;
}

}  // namespace

Window::Window(int width, int height, const std::string& title) {
    glfwSetErrorCallback(&Window::ErrorCallback);
    if (glfwInit() == GLFW_FALSE) {
        throw std::runtime_error("Failed to initialize GLFW");
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
#endif

    window_ = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
    if (window_ == nullptr) {
        glfwTerminate();
        throw std::runtime_error("Failed to create GLFW window");
    }

    glfwMakeContextCurrent(window_);
    glfwSwapInterval(1);

    if (gladLoadGL(reinterpret_cast<GLADloadfunc>(glfwGetProcAddress)) == 0) {
        glfwDestroyWindow(window_);
        glfwTerminate();
        throw std::runtime_error("Failed to load OpenGL functions via glad2");
    }

    WindowMap()[window_] = this;
    glfwSetFramebufferSizeCallback(window_, &Window::FramebufferSizeCallback);
    glfwGetFramebufferSize(window_, &framebuffer_width_, &framebuffer_height_);
}

Window::~Window() {
    if (window_ != nullptr) {
        WindowMap().erase(window_);
        glfwDestroyWindow(window_);
        window_ = nullptr;
    }
    glfwTerminate();
}

bool Window::shouldClose() const {
    return glfwWindowShouldClose(window_) != 0;
}

void Window::pollEvents() const {
    glfwPollEvents();
}

void Window::swapBuffers() const {
    glfwSwapBuffers(window_);
}

float Window::content_scale() const {
    float x_scale = 1.0f;
    float y_scale = 1.0f;
    glfwGetWindowContentScale(window_, &x_scale, &y_scale);
    (void)y_scale;
    return x_scale;
}

int Window::framebuffer_to_window_ratio() const {
    int window_width = 1;
    int window_height = 1;
    glfwGetWindowSize(window_, &window_width, &window_height);
    (void)window_height;
    if (window_width <= 0) {
        return 1;
    }
    return std::max(1, static_cast<int>(std::lround(static_cast<double>(framebuffer_width_) /
                                                    static_cast<double>(window_width))));
}

void Window::ErrorCallback(int error, const char* description) {
    SBOX_LOG_ERROR("GLFW error [%d]: %s", error, description != nullptr ? description : "unknown");
}

void Window::FramebufferSizeCallback(GLFWwindow* window, int width, int height) {
    Window* self = FromGlfwWindow(window);
    if (self != nullptr) {
        self->framebuffer_width_ = width;
        self->framebuffer_height_ = height;
    }
}

Window* Window::FromGlfwWindow(GLFWwindow* window) {
    const auto it = WindowMap().find(window);
    if (it == WindowMap().end()) {
        return nullptr;
    }
    return it->second;
}

}  // namespace sbox
