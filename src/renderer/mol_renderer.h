#pragma once

#include "core/molecular_system.h"
#include "renderer/shader.h"

#include <Eigen/Core>

#include <memory>
#include <vector>

namespace sbox::render {

enum class MolRenderMode : int {
    BallAndStick = 0,
    SpaceFilling = 1,
    Wireframe = 2,
    StickOnly = 3
};

class MolRenderer {
public:
    MolRenderer();
    ~MolRenderer();

    MolRenderer(const MolRenderer&) = delete;
    MolRenderer& operator=(const MolRenderer&) = delete;

    void upload(const sbox::chem::MolecularSystem& mol);

    void render(const Eigen::Matrix4f& view_matrix,
                const Eigen::Matrix4f& proj_matrix,
                const Eigen::Vector3f& camera_pos,
                MolRenderMode mode = MolRenderMode::BallAndStick);

    bool has_data() const;
    int num_atoms() const;
    int num_bonds() const;

private:
    unsigned int atom_vao_ = 0;
    unsigned int atom_vbo_ = 0;
    unsigned int atom_instance_vbo_ = 0;
    int atom_count_ = 0;

    unsigned int bond_vao_ = 0;
    unsigned int bond_vbo_ = 0;
    unsigned int bond_instance_vbo_ = 0;
    int bond_count_ = 0;

    std::unique_ptr<sbox::Shader> atom_shader_;
    std::unique_ptr<sbox::Shader> bond_shader_;

    std::vector<float> atom_instances_;
    std::vector<float> bond_instances_;
};

}  // namespace sbox::render
