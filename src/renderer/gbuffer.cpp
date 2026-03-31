#include "renderer/gbuffer.h"

#include <glad/gl.h>

#include <array>
#include <stdexcept>

namespace sbox::render {

namespace {

void configure_gbuffer_texture(unsigned int texture) {
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}

}

GBuffer::GBuffer() = default;

GBuffer::~GBuffer() {
    destroy();
}

void GBuffer::resize(int width, int height) {
    width = width < 1 ? 1 : width;
    height = height < 1 ? 1 : height;

    if (width == width_ && height == height_ && fbo_ != 0U) {
        return;
    }

    destroy();

    width_ = width;
    height_ = height;

    glGenFramebuffers(1, &fbo_);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo_);

    glGenTextures(1, &position_tex_);
    configure_gbuffer_texture(position_tex_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, width_, height_, 0, GL_RGB, GL_FLOAT, nullptr);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, position_tex_, 0);

    glGenTextures(1, &normal_tex_);
    configure_gbuffer_texture(normal_tex_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, width_, height_, 0, GL_RGB, GL_FLOAT, nullptr);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, normal_tex_, 0);

    glGenTextures(1, &albedo_tex_);
    configure_gbuffer_texture(albedo_tex_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width_, height_, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, albedo_tex_, 0);

    glGenTextures(1, &depth_tex_);
    configure_gbuffer_texture(depth_tex_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, width_, height_, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depth_tex_, 0);

    constexpr std::array<unsigned int, 3> attachments = {
        GL_COLOR_ATTACHMENT0,
        GL_COLOR_ATTACHMENT1,
        GL_COLOR_ATTACHMENT2,
    };
    glDrawBuffers(static_cast<int>(attachments.size()), attachments.data());

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        destroy();
        throw std::runtime_error("G-buffer framebuffer is incomplete");
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void GBuffer::bind_for_writing() const {
    glBindFramebuffer(GL_FRAMEBUFFER, fbo_);
}

void GBuffer::bind_position_texture(int unit) const {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, position_tex_);
}

void GBuffer::bind_normal_texture(int unit) const {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, normal_tex_);
}

void GBuffer::bind_albedo_texture(int unit) const {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, albedo_tex_);
}

void GBuffer::bind_depth_texture(int unit) const {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, depth_tex_);
}

void GBuffer::unbind() const {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

unsigned int GBuffer::fbo() const {
    return fbo_;
}

unsigned int GBuffer::position_texture() const {
    return position_tex_;
}

unsigned int GBuffer::normal_texture() const {
    return normal_tex_;
}

unsigned int GBuffer::albedo_texture() const {
    return albedo_tex_;
}

unsigned int GBuffer::depth_texture() const {
    return depth_tex_;
}

int GBuffer::width() const {
    return width_;
}

int GBuffer::height() const {
    return height_;
}

void GBuffer::destroy() {
    if (depth_rbo_ != 0U) {
        glDeleteRenderbuffers(1, &depth_rbo_);
        depth_rbo_ = 0;
    }
    if (depth_tex_ != 0U) {
        glDeleteTextures(1, &depth_tex_);
        depth_tex_ = 0;
    }
    if (albedo_tex_ != 0U) {
        glDeleteTextures(1, &albedo_tex_);
        albedo_tex_ = 0;
    }
    if (normal_tex_ != 0U) {
        glDeleteTextures(1, &normal_tex_);
        normal_tex_ = 0;
    }
    if (position_tex_ != 0U) {
        glDeleteTextures(1, &position_tex_);
        position_tex_ = 0;
    }
    if (fbo_ != 0U) {
        glDeleteFramebuffers(1, &fbo_);
        fbo_ = 0;
    }
    width_ = 0;
    height_ = 0;
}

}  // namespace sbox::render
