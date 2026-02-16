#pragma once

#include <string>

struct GLFWwindow;

namespace sbox {

class Window {
public:
    Window(int width, int height, const std::string& title);
    ~Window();

    Window(const Window&) = delete;
    Window& operator=(const Window&) = delete;

    [[nodiscard]] GLFWwindow* handle() const { return window_; }
    [[nodiscard]] bool shouldClose() const;

    void pollEvents() const;
    void swapBuffers() const;

    [[nodiscard]] int framebufferWidth() const { return framebuffer_width_; }
    [[nodiscard]] int framebufferHeight() const { return framebuffer_height_; }

private:
    static void ErrorCallback(int error, const char* description);
    static void FramebufferSizeCallback(GLFWwindow* window, int width, int height);
    static Window* FromGlfwWindow(GLFWwindow* window);

    GLFWwindow* window_ = nullptr;
    int framebuffer_width_ = 0;
    int framebuffer_height_ = 0;
};

}  // namespace sbox
