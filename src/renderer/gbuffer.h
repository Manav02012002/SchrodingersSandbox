#pragma once

namespace sbox::render {

class GBuffer {
public:
    GBuffer();
    ~GBuffer();

    void resize(int width, int height);

    void bind_for_writing() const;
    void bind_position_texture(int unit = 0) const;
    void bind_normal_texture(int unit = 1) const;
    void bind_albedo_texture(int unit = 2) const;
    void bind_depth_texture(int unit = 3) const;

    void unbind() const;

    unsigned int fbo() const;
    unsigned int position_texture() const;
    unsigned int normal_texture() const;
    unsigned int albedo_texture() const;
    unsigned int depth_texture() const;

    int width() const;
    int height() const;

private:
    void destroy();

    unsigned int fbo_ = 0;
    unsigned int position_tex_ = 0;
    unsigned int normal_tex_ = 0;
    unsigned int albedo_tex_ = 0;
    unsigned int depth_rbo_ = 0;
    unsigned int depth_tex_ = 0;
    int width_ = 0;
    int height_ = 0;
};

}  // namespace sbox::render
