#pragma once

#include "core/molecular_system.h"
#include "renderer/mol_renderer.h"
#include "renderer/shader.h"

#include <Eigen/Core>

#include <memory>
#include <vector>

namespace sbox::render {

class ShadowMap {
public:
    ShadowMap();
    ~ShadowMap();

    void init(int resolution = 2048);

    void compute(MolRenderer& mol_renderer,
                 const sbox::chem::MolecularSystem& mol,
                 const Eigen::Vector3f& scene_center,
                 float scene_radius,
                 const Eigen::Vector3f& light_dir);

    void bind_shadow_texture(int unit = 5) const;
    Eigen::Matrix4f light_vp_matrix() const;
    bool is_initialized() const;
    int resolution() const;

private:
    void destroy();
    void create_sphere_mesh();

    unsigned int shadow_fbo_ = 0;
    unsigned int shadow_tex_ = 0;
    unsigned int sphere_vao_ = 0;
    unsigned int sphere_vbo_ = 0;
    unsigned int sphere_ebo_ = 0;
    unsigned int instance_vbo_ = 0;
    int sphere_index_count_ = 0;
    int resolution_ = 2048;
    Eigen::Matrix4f light_vp_ = Eigen::Matrix4f::Identity();
    bool initialized_ = false;

    std::unique_ptr<Shader> depth_shader_;
};

}  // namespace sbox::render
