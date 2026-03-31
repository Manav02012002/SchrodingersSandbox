// HiDPI manual test checklist:
// 1. On Retina Mac: verify text is sharp, not blurry
// 2. On Retina Mac: verify orbital rendering is sharp (not pixelated)
// 3. On Retina Mac: verify mouse picking works (click on atoms hits correctly)
// 4. On Retina Mac: verify screenshots at viewport resolution produce 2x pixel images
// 5. On non-Retina display: verify everything still works at scale=1.0
// 6. Moving window between Retina and non-Retina: verify resize and re-render

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
    [[nodiscard]] float content_scale() const;
    [[nodiscard]] int framebuffer_to_window_ratio() const;

private:
    static void ErrorCallback(int error, const char* description);
    static void FramebufferSizeCallback(GLFWwindow* window, int width, int height);
    static Window* FromGlfwWindow(GLFWwindow* window);

    GLFWwindow* window_ = nullptr;
    int framebuffer_width_ = 0;
    int framebuffer_height_ = 0;
};

}  // namespace sbox
