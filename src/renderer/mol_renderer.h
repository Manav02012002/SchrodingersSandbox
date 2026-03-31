#pragma once

#include "core/molecular_system.h"
#include "editor/picking.h"
#include "io/pdb_io.h"
#include "renderer/shader.h"

#include <Eigen/Core>

#include <imgui.h>

#include <memory>
#include <vector>

namespace sbox::render {

enum class MolRenderMode : int {
    BallAndStick = 0,
    SpaceFilling = 1,
    Wireframe = 2,
    StickOnly = 3
};

enum class ColorMode : int {
    CPK = 0,
    ByChain = 1,
    ByResidue = 2,
    ByBFactor = 3,
    ByCharge = 4,
    BySecondary = 5,
    Custom = 6,
};

Eigen::Vector3f chain_color(int chain_index);
Eigen::Vector3f residue_color(const std::string& residue_name);
Eigen::Vector3f bfactor_color(float b_factor, float b_min, float b_max);
void set_atom_radius_scale(float scale);
void set_bond_radius_scale(float scale);

class MolRenderer {
public:
    MolRenderer();
    ~MolRenderer();

    MolRenderer(const MolRenderer&) = delete;
    MolRenderer& operator=(const MolRenderer&) = delete;

    void upload(const sbox::chem::MolecularSystem& mol,
                ColorMode color_mode = ColorMode::CPK,
                const sbox::io::PDBData* pdb_data = nullptr,
                const std::vector<double>* charges = nullptr);

    void render(const Eigen::Matrix4f& view_matrix,
                const Eigen::Matrix4f& proj_matrix,
                const Eigen::Vector3f& camera_pos,
                MolRenderMode mode = MolRenderMode::BallAndStick);
    void render_gbuffer(const Eigen::Matrix4f& view_matrix,
                        const Eigen::Matrix4f& proj_matrix,
                        const Eigen::Vector3f& camera_pos,
                        MolRenderMode mode = MolRenderMode::BallAndStick);
    void render_selection(const Eigen::Matrix4f& view_matrix,
                          const Eigen::Matrix4f& proj_matrix,
                          const Eigen::Vector3f& camera_pos,
                          const sbox::chem::MolecularSystem& mol,
                          const sbox::editor::Selection& selection,
                          MolRenderMode mode);
    void render_atom_labels(const sbox::chem::MolecularSystem& mol,
                            const Eigen::Matrix4f& vp_matrix,
                            const ImVec2& viewport_pos,
                            const ImVec2& viewport_size,
                            bool show_indices = false,
                            bool show_symbols = true,
                            bool show_hydrogens = false);

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
    std::unique_ptr<sbox::Shader> gbuffer_atom_shader_;
    std::unique_ptr<sbox::Shader> gbuffer_bond_shader_;

    std::vector<float> atom_instances_;
    std::vector<float> bond_instances_;

    void render_highlights(const Eigen::Matrix4f& view_matrix,
                           const Eigen::Matrix4f& proj_matrix,
                           const Eigen::Vector3f& camera_pos,
                           const std::vector<float>& atom_data,
                           const std::vector<float>& bond_data,
                           const Eigen::Vector3f& highlight_color,
                           MolRenderMode mode);
};

}  // namespace sbox::render
