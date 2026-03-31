#include "renderer/ssao.h"

#include "core/paths.h"

#include <glad/gl.h>

#include <algorithm>
#include <random>
#include <stdexcept>
#include <vector>

namespace sbox::render {

namespace {

float random_float(std::mt19937& rng, float min_value, float max_value) {
    std::uniform_real_distribution<float> distribution(min_value, max_value);
    return distribution(rng);
}

void configure_single_channel_texture(unsigned int texture, int width, int height) {
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R16F, width, height, 0, GL_RED, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}

}  // namespace

SSAO::SSAO() = default;

SSAO::~SSAO() {
    destroy();
}

void SSAO::init(int width, int height) {
    destroy();

    ssao_shader_ = std::make_unique<Shader>(sbox::get_shader_path("fullscreen_quad.vert"),
                                            sbox::get_shader_path("ssao.frag"));
    blur_shader_ = std::make_unique<Shader>(sbox::get_shader_path("fullscreen_quad.vert"),
                                            sbox::get_shader_path("ssao_blur.frag"));

    glGenVertexArrays(1, &fullscreen_vao_);
    glGenFramebuffers(1, &ssao_fbo_);
    glGenTextures(1, &ssao_tex_);
    glGenFramebuffers(1, &blur_fbo_);
    glGenTextures(1, &blur_tex_);

    generate_kernel();
    generate_noise_texture();
    initialized_ = true;
    resize(width, height);
}

void SSAO::resize(int width, int height) {
    width = std::max(width, 1);
    height = std::max(height, 1);

    if (!initialized_) {
        width_ = width;
        height_ = height;
        return;
    }

    if (width_ == width && height_ == height) {
        return;
    }

    width_ = width;
    height_ = height;

    glBindFramebuffer(GL_FRAMEBUFFER, ssao_fbo_);
    configure_single_channel_texture(ssao_tex_, width_, height_);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ssao_tex_, 0);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        throw std::runtime_error("SSAO framebuffer is incomplete");
    }

    glBindFramebuffer(GL_FRAMEBUFFER, blur_fbo_);
    configure_single_channel_texture(blur_tex_, width_, height_);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, blur_tex_, 0);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        throw std::runtime_error("SSAO blur framebuffer is incomplete");
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void SSAO::compute(const GBuffer& gbuffer, const Eigen::Matrix4f& projection) {
    if (!initialized_) {
        init(gbuffer.width(), gbuffer.height());
    } else if (width_ != gbuffer.width() || height_ != gbuffer.height()) {
        resize(gbuffer.width(), gbuffer.height());
    }

    glBindFramebuffer(GL_FRAMEBUFFER, ssao_fbo_);
    glViewport(0, 0, width_, height_);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    ssao_shader_->bind();
    ssao_shader_->setUniform("u_position", 0);
    ssao_shader_->setUniform("u_normal", 1);
    ssao_shader_->setUniform("u_noise", 2);
    ssao_shader_->setUniform("u_projection", projection);
    ssao_shader_->setUniform("u_kernel_size", static_cast<int>(kernel_.size()));
    ssao_shader_->setUniform("u_radius", radius_);
    ssao_shader_->setUniform("u_bias", 0.025f);
    ssao_shader_->setUniform("u_power", power_);

    const int noise_scale_loc = glGetUniformLocation(ssao_shader_->id(), "u_noise_scale");
    if (noise_scale_loc >= 0) {
        glUniform2f(noise_scale_loc, static_cast<float>(width_) / 4.0f, static_cast<float>(height_) / 4.0f);
    }

    for (std::size_t i = 0; i < kernel_.size(); ++i) {
        const std::string uniform_name = "u_samples[" + std::to_string(i) + "]";
        const int sample_loc = glGetUniformLocation(ssao_shader_->id(), uniform_name.c_str());
        if (sample_loc >= 0) {
            glUniform3f(sample_loc, kernel_[i].x(), kernel_[i].y(), kernel_[i].z());
        }
    }

    gbuffer.bind_position_texture(0);
    gbuffer.bind_normal_texture(1);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, noise_tex_);

    glBindVertexArray(fullscreen_vao_);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    glBindFramebuffer(GL_FRAMEBUFFER, blur_fbo_);
    glViewport(0, 0, width_, height_);
    glClear(GL_COLOR_BUFFER_BIT);

    blur_shader_->bind();
    blur_shader_->setUniform("u_ssao_input", 0);

    const int texel_loc = glGetUniformLocation(blur_shader_->id(), "u_texel_size");
    if (texel_loc >= 0) {
        glUniform2f(texel_loc, 1.0f / static_cast<float>(width_), 1.0f / static_cast<float>(height_));
    }

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, ssao_tex_);
    glBindVertexArray(fullscreen_vao_);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    glBindTexture(GL_TEXTURE_2D, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void SSAO::bind_ao_texture(int unit) const {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, blur_tex_);
}

bool SSAO::is_initialized() const {
    return initialized_;
}

void SSAO::set_parameters(float radius, float power) {
    radius_ = std::max(radius, 0.01f);
    power_ = std::max(power, 0.01f);
}

void SSAO::destroy() {
    initialized_ = false;
    width_ = 0;
    height_ = 0;
    kernel_.clear();
    blur_shader_.reset();
    ssao_shader_.reset();

    if (noise_tex_ != 0U) {
        glDeleteTextures(1, &noise_tex_);
        noise_tex_ = 0;
    }
    if (blur_tex_ != 0U) {
        glDeleteTextures(1, &blur_tex_);
        blur_tex_ = 0;
    }
    if (blur_fbo_ != 0U) {
        glDeleteFramebuffers(1, &blur_fbo_);
        blur_fbo_ = 0;
    }
    if (ssao_tex_ != 0U) {
        glDeleteTextures(1, &ssao_tex_);
        ssao_tex_ = 0;
    }
    if (ssao_fbo_ != 0U) {
        glDeleteFramebuffers(1, &ssao_fbo_);
        ssao_fbo_ = 0;
    }
    if (fullscreen_vao_ != 0U) {
        glDeleteVertexArrays(1, &fullscreen_vao_);
        fullscreen_vao_ = 0;
    }
}

void SSAO::generate_kernel(int num_samples) {
    kernel_.clear();
    kernel_.reserve(static_cast<std::size_t>(num_samples));

    std::mt19937 rng(42);
    for (int i = 0; i < num_samples; ++i) {
        Eigen::Vector3f sample(random_float(rng, -1.0f, 1.0f),
                               random_float(rng, -1.0f, 1.0f),
                               random_float(rng, 0.0f, 1.0f));
        if (sample.squaredNorm() < 1.0e-6f) {
            sample = Eigen::Vector3f(0.0f, 0.0f, 1.0f);
        } else {
            sample.normalize();
        }
        sample *= random_float(rng, 0.0f, 1.0f);

        float scale = static_cast<float>(i) / static_cast<float>(num_samples);
        scale = 0.1f + scale * scale * 0.9f;
        sample *= scale;

        kernel_.push_back(sample);
    }
}

void SSAO::generate_noise_texture() {
    std::mt19937 rng(42);
    std::vector<Eigen::Vector3f> noise(16);
    for (Eigen::Vector3f& value : noise) {
        value = Eigen::Vector3f(random_float(rng, -1.0f, 1.0f), random_float(rng, -1.0f, 1.0f), 0.0f);
        if (value.head<2>().squaredNorm() < 1.0e-6f) {
            value = Eigen::Vector3f(1.0f, 0.0f, 0.0f);
        } else {
            value.normalize();
        }
    }

    glGenTextures(1, &noise_tex_);
    glBindTexture(GL_TEXTURE_2D, noise_tex_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, 4, 4, 0, GL_RGB, GL_FLOAT, noise.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glBindTexture(GL_TEXTURE_2D, 0);
}

}  // namespace sbox::render
