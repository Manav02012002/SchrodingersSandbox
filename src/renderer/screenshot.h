#pragma once

#include <functional>
#include <string>

namespace sbox::render {

bool save_screenshot(const std::string& filepath, unsigned int fbo, int width, int height);

bool save_screenshot_highres(
    const std::string& filepath,
    int render_width,
    int render_height,
    std::function<void(unsigned int fbo, int w, int h)> render_fn);

bool save_screenshot_transparent(
    const std::string& filepath,
    unsigned int fbo,
    int width,
    int height);

}  // namespace sbox::render
