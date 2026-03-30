#pragma once

#include "core/molecular_system.h"
#include "io/pdb_io.h"
#include "renderer/mol_renderer.h"
#include "renderer/shader.h"

#include <Eigen/Core>

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace sbox::render {

class SpatialHash {
public:
    explicit SpatialHash(float cell_size = 10.0f);

    void build(const sbox::chem::MolecularSystem& mol);
    void clear();

    std::vector<int> query_sphere(const Eigen::Vector3f& center, float radius) const;
    std::vector<int> query_frustum(const Eigen::Matrix4f& vp_matrix) const;

    int count() const;

private:
    struct Cell {
        std::vector<int> atom_indices;
        int ix = 0;
        int iy = 0;
        int iz = 0;
    };

    float cell_size_;
    std::unordered_map<int64_t, Cell> cells_;
    std::vector<Eigen::Vector3f> positions_;

    int64_t hash_position(float x, float y, float z) const;
};

class LODRenderer {
public:
    LODRenderer();
    ~LODRenderer();

    void rebuild(const sbox::chem::MolecularSystem& mol,
                 ColorMode color_mode = ColorMode::CPK,
                 const sbox::io::PDBData* pdb_data = nullptr,
                 const std::vector<double>* charges = nullptr);

    void render(const Eigen::Matrix4f& view_matrix,
                const Eigen::Matrix4f& proj_matrix,
                const Eigen::Vector3f& camera_pos,
                MolRenderMode mode,
                const sbox::chem::MolecularSystem& mol,
                ColorMode color_mode = ColorMode::CPK,
                const sbox::io::PDBData* pdb_data = nullptr,
                const std::vector<double>* charges = nullptr);

    int atoms_rendered() const;
    int atoms_culled() const;
    int bonds_rendered() const;

private:
    SpatialHash spatial_hash_;
    MolRenderer full_renderer_;
    unsigned int point_vao_ = 0;
    unsigned int point_vbo_ = 0;
    std::unique_ptr<sbox::Shader> point_shader_;

    int atoms_rendered_ = 0;
    int atoms_culled_ = 0;
    int bonds_rendered_ = 0;
    ColorMode color_mode_ = ColorMode::CPK;
    const sbox::io::PDBData* pdb_data_ = nullptr;
    const std::vector<double>* charges_ = nullptr;
};

}  // namespace sbox::render
