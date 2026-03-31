#include "renderer/post_process.h"

#include "core/paths.h"

#include <glad/gl.h>

#include <algorithm>
#include <stdexcept>

namespace sbox::render {

namespace {

void allocate_color_target(unsigned int texture, int width, int height) {
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}

}

PostProcess::PostProcess() = default;

PostProcess::~PostProcess() {
    destroy();
}

void PostProcess::init(int width, int height) {
    destroy();

    fxaa_shader_ = std::make_unique<Shader>(sbox::get_shader_path("fullscreen_quad.vert"),
                                            sbox::get_shader_path("fxaa.frag"));
    tonemap_shader_ = std::make_unique<Shader>(sbox::get_shader_path("fullscreen_quad.vert"),
                                               sbox::get_shader_path("tonemap.frag"));
    glGenVertexArrays(1, &fullscreen_vao_);
    glGenFramebuffers(1, &ping_fbo_);
    glGenTextures(1, &ping_tex_);
    glGenFramebuffers(1, &pong_fbo_);
    glGenTextures(1, &pong_tex_);
    use_ping_ = true;
    initialized_ = true;
    ensure_fbos(width, height);
}

void PostProcess::resize(int width, int height) {
    if (!initialized_) {
        init(width, height);
        return;
    }
    ensure_fbos(width, height);
}

void PostProcess::apply_fxaa(unsigned int input_texture, int width, int height) {
    if (!initialized_) {
        init(width, height);
    } else {
        ensure_fbos(width, height);
    }

    glBindFramebuffer(GL_FRAMEBUFFER, current_output_fbo());
    glViewport(0, 0, width_, height_);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    fxaa_shader_->bind();
    fxaa_shader_->setUniform("u_input", 0);
    const int texel_loc = glGetUniformLocation(fxaa_shader_->id(), "u_texel_size");
    if (texel_loc >= 0) {
        glUniform2f(texel_loc, 1.0f / static_cast<float>(width_), 1.0f / static_cast<float>(height_));
    }

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, input_texture);
    glBindVertexArray(fullscreen_vao_);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    swap_buffers();
}

void PostProcess::apply_tonemap(unsigned int input_texture, int width, int height, float exposure, float gamma) {
    if (!initialized_) {
        init(width, height);
    } else {
        ensure_fbos(width, height);
    }

    glBindFramebuffer(GL_FRAMEBUFFER, current_output_fbo());
    glViewport(0, 0, width_, height_);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    tonemap_shader_->bind();
    tonemap_shader_->setUniform("u_input", 0);
    tonemap_shader_->setUniform("u_exposure", exposure);
    tonemap_shader_->setUniform("u_gamma", gamma);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, input_texture);
    glBindVertexArray(fullscreen_vao_);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    swap_buffers();
}

unsigned int PostProcess::output_texture() const {
    return current_input();
}

void PostProcess::blit_to_screen(unsigned int fullscreen_vao) {
    (void)fullscreen_vao;

    GLint draw_fbo = 0;
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &draw_fbo);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, use_ping_ ? pong_fbo_ : ping_fbo_);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, static_cast<unsigned int>(draw_fbo));
    glBlitFramebuffer(0, 0, width_, height_,
                      0, 0, width_, height_,
                      GL_COLOR_BUFFER_BIT,
                      GL_NEAREST);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
}

bool PostProcess::is_initialized() const {
    return initialized_;
}

void PostProcess::destroy() {
    initialized_ = false;
    width_ = 0;
    height_ = 0;
    use_ping_ = true;
    tonemap_shader_.reset();
    fxaa_shader_.reset();

    if (pong_tex_ != 0U) {
        glDeleteTextures(1, &pong_tex_);
        pong_tex_ = 0;
    }
    if (pong_fbo_ != 0U) {
        glDeleteFramebuffers(1, &pong_fbo_);
        pong_fbo_ = 0;
    }
    if (ping_tex_ != 0U) {
        glDeleteTextures(1, &ping_tex_);
        ping_tex_ = 0;
    }
    if (ping_fbo_ != 0U) {
        glDeleteFramebuffers(1, &ping_fbo_);
        ping_fbo_ = 0;
    }
    if (fullscreen_vao_ != 0U) {
        glDeleteVertexArrays(1, &fullscreen_vao_);
        fullscreen_vao_ = 0;
    }
}

void PostProcess::ensure_fbos(int width, int height) {
    width = std::max(width, 1);
    height = std::max(height, 1);
    if (width_ == width && height_ == height) {
        return;
    }

    width_ = width;
    height_ = height;

    glBindFramebuffer(GL_FRAMEBUFFER, ping_fbo_);
    allocate_color_target(ping_tex_, width_, height_);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ping_tex_, 0);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        throw std::runtime_error("Post-process ping framebuffer is incomplete");
    }

    glBindFramebuffer(GL_FRAMEBUFFER, pong_fbo_);
    allocate_color_target(pong_tex_, width_, height_);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, pong_tex_, 0);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        throw std::runtime_error("Post-process pong framebuffer is incomplete");
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

unsigned int PostProcess::current_input() const {
    return use_ping_ ? pong_tex_ : ping_tex_;
}

unsigned int PostProcess::current_output_fbo() const {
    return use_ping_ ? ping_fbo_ : pong_fbo_;
}

void PostProcess::swap_buffers() {
    use_ping_ = !use_ping_;
}

}  // namespace sbox::render
