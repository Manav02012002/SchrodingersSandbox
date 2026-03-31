#pragma once

#include "renderer/shader.h"

#include <memory>

namespace sbox::render {

class PostProcess {
public:
    PostProcess();
    ~PostProcess();

    void init(int width, int height);
    void resize(int width, int height);

    void apply_fxaa(unsigned int input_texture, int width, int height);
    void apply_tonemap(unsigned int input_texture, int width, int height,
                       float exposure = 1.0f, float gamma = 2.2f);

    unsigned int output_texture() const;
    void blit_to_screen(unsigned int fullscreen_vao);
    bool is_initialized() const;

private:
    void destroy();
    void ensure_fbos(int width, int height);
    unsigned int current_input() const;
    unsigned int current_output_fbo() const;
    void swap_buffers();

    unsigned int ping_fbo_ = 0;
    unsigned int ping_tex_ = 0;
    unsigned int pong_fbo_ = 0;
    unsigned int pong_tex_ = 0;
    bool use_ping_ = true;

    std::unique_ptr<Shader> fxaa_shader_;
    std::unique_ptr<Shader> tonemap_shader_;
    unsigned int fullscreen_vao_ = 0;

    int width_ = 0;
    int height_ = 0;
    bool initialized_ = false;
};

}  // namespace sbox::render
