#pragma once

#include "renderer/gbuffer.h"
#include "renderer/shader.h"

#include <Eigen/Core>

#include <memory>
#include <vector>

namespace sbox::render {

class SSAO {
public:
    SSAO();
    ~SSAO();

    void init(int width, int height);
    void resize(int width, int height);
    void compute(const GBuffer& gbuffer, const Eigen::Matrix4f& projection);
    void bind_ao_texture(int unit = 4) const;
    void set_parameters(float radius, float power);

    bool is_initialized() const;

private:
    void destroy();
    void generate_kernel(int num_samples = 32);
    void generate_noise_texture();

    unsigned int ssao_fbo_ = 0;
    unsigned int ssao_tex_ = 0;
    unsigned int blur_fbo_ = 0;
    unsigned int blur_tex_ = 0;
    unsigned int noise_tex_ = 0;

    std::vector<Eigen::Vector3f> kernel_;

    std::unique_ptr<Shader> ssao_shader_;
    std::unique_ptr<Shader> blur_shader_;

    unsigned int fullscreen_vao_ = 0;
    int width_ = 0;
    int height_ = 0;
    bool initialized_ = false;
    float radius_ = 1.5f;
    float power_ = 2.0f;
};

}  // namespace sbox::render
