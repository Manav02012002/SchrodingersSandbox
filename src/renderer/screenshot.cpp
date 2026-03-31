#include "renderer/screenshot.h"

#include <glad/gl.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <string>
#include <vector>

namespace sbox::render {

namespace {

void flip_rgba_rows(std::vector<unsigned char>& pixels, int width, int height) {
    const std::size_t row_bytes = static_cast<std::size_t>(width) * 4u;
    std::vector<unsigned char> temp(row_bytes);
    for (int y = 0; y < height / 2; ++y) {
        unsigned char* top = pixels.data() + static_cast<std::size_t>(y) * row_bytes;
        unsigned char* bottom = pixels.data() + static_cast<std::size_t>(height - 1 - y) * row_bytes;
        std::copy(top, top + row_bytes, temp.data());
        std::copy(bottom, bottom + row_bytes, top);
        std::copy(temp.data(), temp.data() + row_bytes, bottom);
    }
}

std::string lowercase_extension(const std::string& filepath) {
    std::string ext = std::filesystem::path(filepath).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return ext;
}

bool write_image(const std::string& filepath, int width, int height, const unsigned char* pixels) {
    const std::string ext = lowercase_extension(filepath);
    if (ext == ".jpg" || ext == ".jpeg") {
        return stbi_write_jpg(filepath.c_str(), width, height, 4, pixels, 95) != 0;
    }
    if (ext == ".bmp") {
        return stbi_write_bmp(filepath.c_str(), width, height, 4, pixels) != 0;
    }
    return stbi_write_png(filepath.c_str(), width, height, 4, pixels, width * 4) != 0;
}

bool read_and_write(const std::string& filepath, unsigned int fbo, int width, int height) {
    if (width <= 0 || height <= 0) {
        return false;
    }

    GLint prev_fbo = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    std::vector<unsigned char> pixels(static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4u);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
    glBindFramebuffer(GL_FRAMEBUFFER, static_cast<unsigned int>(prev_fbo));

    flip_rgba_rows(pixels, width, height);
    return write_image(filepath, width, height, pixels.data());
}

}  // namespace

bool save_screenshot(const std::string& filepath, unsigned int fbo, int width, int height) {
    return read_and_write(filepath, fbo, width, height);
}

bool save_screenshot_highres(
    const std::string& filepath,
    int render_width,
    int render_height,
    std::function<void(unsigned int fbo, int w, int h)> render_fn) {
    if (render_width <= 0 || render_height <= 0 || !render_fn) {
        return false;
    }

    unsigned int temp_fbo = 0;
    unsigned int color_tex = 0;
    unsigned int depth_rbo = 0;

    glGenFramebuffers(1, &temp_fbo);
    glGenTextures(1, &color_tex);
    glGenRenderbuffers(1, &depth_rbo);

    glBindTexture(GL_TEXTURE_2D, color_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, render_width, render_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glBindRenderbuffer(GL_RENDERBUFFER, depth_rbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, render_width, render_height);

    glBindFramebuffer(GL_FRAMEBUFFER, temp_fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, color_tex, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depth_rbo);
    const bool complete = glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;

    bool ok = false;
    if (complete) {
        render_fn(temp_fbo, render_width, render_height);
        ok = read_and_write(filepath, temp_fbo, render_width, render_height);
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glDeleteRenderbuffers(1, &depth_rbo);
    glDeleteTextures(1, &color_tex);
    glDeleteFramebuffers(1, &temp_fbo);
    return ok;
}

bool save_screenshot_transparent(
    const std::string& filepath,
    unsigned int fbo,
    int width,
    int height) {
    return read_and_write(filepath, fbo, width, height);
}

}  // namespace sbox::render
