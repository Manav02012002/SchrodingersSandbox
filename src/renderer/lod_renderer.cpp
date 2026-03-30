#include "renderer/lod_renderer.h"

#include "core/elements.h"

#include <glad/gl.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace sbox::render {

namespace {

constexpr float kLodThreshold = 50.0f;
constexpr int kLodAtomCountThreshold = 200;

Eigen::Vector3f cpk_color(int Z) {
    switch (Z) {
    case 1: return Eigen::Vector3f(1.0f, 1.0f, 1.0f);
    case 6: return Eigen::Vector3f(0.56f, 0.56f, 0.56f);
    case 7: return Eigen::Vector3f(0.19f, 0.31f, 0.97f);
    case 8: return Eigen::Vector3f(1.0f, 0.05f, 0.05f);
    case 9: return Eigen::Vector3f(0.56f, 0.88f, 0.31f);
    case 15: return Eigen::Vector3f(1.0f, 0.50f, 0.0f);
    case 16: return Eigen::Vector3f(1.0f, 1.0f, 0.19f);
    case 17: return Eigen::Vector3f(0.12f, 0.94f, 0.12f);
    case 26: return Eigen::Vector3f(0.88f, 0.40f, 0.20f);
    default: return Eigen::Vector3f(0.55f, 0.60f, 0.68f);
    }
}

Eigen::Vector3f charge_color(double charge, double max_abs_charge) {
    if (max_abs_charge <= 1.0e-8) {
        return Eigen::Vector3f(1.0f, 1.0f, 1.0f);
    }
    const float t = static_cast<float>(std::clamp(charge / max_abs_charge, -1.0, 1.0));
    if (t < 0.0f) {
        return Eigen::Vector3f(1.0f, 1.0f, 1.0f) * (t + 1.0f) + Eigen::Vector3f(0.15f, 0.35f, 0.95f) * (-t);
    }
    return Eigen::Vector3f(1.0f, 1.0f, 1.0f) * (1.0f - t) + Eigen::Vector3f(0.95f, 0.20f, 0.15f) * t;
}

float point_size_for_element(int Z) {
    switch (Z) {
    case 1: return 6.0f;
    case 6: return 9.0f;
    case 7: return 9.0f;
    case 8: return 8.5f;
    case 16: return 10.0f;
    default: return 7.5f;
    }
}

struct Plane {
    Eigen::Vector3f normal = Eigen::Vector3f::Zero();
    float d = 0.0f;
};

Plane normalize_plane(const Eigen::Vector4f& coeffs) {
    Plane plane;
    plane.normal = coeffs.head<3>();
    const float len = plane.normal.norm();
    if (len > 1.0e-6f) {
        plane.normal /= len;
        plane.d = coeffs.w() / len;
    } else {
        plane.d = coeffs.w();
    }
    return plane;
}

std::array<Plane, 6> extract_frustum_planes(const Eigen::Matrix4f& vp) {
    const Eigen::Vector4f row0 = vp.row(0);
    const Eigen::Vector4f row1 = vp.row(1);
    const Eigen::Vector4f row2 = vp.row(2);
    const Eigen::Vector4f row3 = vp.row(3);
    return {{
        normalize_plane(row3 + row0),
        normalize_plane(row3 - row0),
        normalize_plane(row3 + row1),
        normalize_plane(row3 - row1),
        normalize_plane(row3 + row2),
        normalize_plane(row3 - row2),
    }};
}

bool aabb_intersects_frustum(const Eigen::Vector3f& aabb_min,
                             const Eigen::Vector3f& aabb_max,
                             const std::array<Plane, 6>& planes) {
    for (const Plane& plane : planes) {
        Eigen::Vector3f positive = aabb_min;
        if (plane.normal.x() >= 0.0f) positive.x() = aabb_max.x();
        if (plane.normal.y() >= 0.0f) positive.y() = aabb_max.y();
        if (plane.normal.z() >= 0.0f) positive.z() = aabb_max.z();
        if (plane.normal.dot(positive) + plane.d < 0.0f) {
            return false;
        }
    }
    return true;
}

}  // namespace

SpatialHash::SpatialHash(float cell_size)
    : cell_size_(cell_size) {}

void SpatialHash::build(const sbox::chem::MolecularSystem& mol) {
    clear();
    positions_.reserve(static_cast<std::size_t>(mol.num_atoms()));
    for (int i = 0; i < mol.num_atoms(); ++i) {
        const Eigen::Vector3f pos = mol.atom(i).position.cast<float>();
        positions_.push_back(pos);
        const int ix = static_cast<int>(std::floor(pos.x() / cell_size_));
        const int iy = static_cast<int>(std::floor(pos.y() / cell_size_));
        const int iz = static_cast<int>(std::floor(pos.z() / cell_size_));
        const int64_t hash = hash_position(pos.x(), pos.y(), pos.z());
        Cell& cell = cells_[hash];
        cell.ix = ix;
        cell.iy = iy;
        cell.iz = iz;
        cell.atom_indices.push_back(i);
    }
}

void SpatialHash::clear() {
    cells_.clear();
    positions_.clear();
}

std::vector<int> SpatialHash::query_sphere(const Eigen::Vector3f& center, float radius) const {
    std::vector<int> result;
    if (positions_.empty()) {
        return result;
    }

    const int min_x = static_cast<int>(std::floor((center.x() - radius) / cell_size_));
    const int min_y = static_cast<int>(std::floor((center.y() - radius) / cell_size_));
    const int min_z = static_cast<int>(std::floor((center.z() - radius) / cell_size_));
    const int max_x = static_cast<int>(std::floor((center.x() + radius) / cell_size_));
    const int max_y = static_cast<int>(std::floor((center.y() + radius) / cell_size_));
    const int max_z = static_cast<int>(std::floor((center.z() + radius) / cell_size_));
    const float radius_sq = radius * radius;

    for (int ix = min_x; ix <= max_x; ++ix) {
        for (int iy = min_y; iy <= max_y; ++iy) {
            for (int iz = min_z; iz <= max_z; ++iz) {
                const int64_t hash = static_cast<int64_t>(ix) +
                                     static_cast<int64_t>(iy) * 10007LL +
                                     static_cast<int64_t>(iz) * 100000007LL;
                const auto it = cells_.find(hash);
                if (it == cells_.end()) {
                    continue;
                }
                for (int atom_index : it->second.atom_indices) {
                    if ((positions_[static_cast<std::size_t>(atom_index)] - center).squaredNorm() <= radius_sq) {
                        result.push_back(atom_index);
                    }
                }
            }
        }
    }
    return result;
}

std::vector<int> SpatialHash::query_frustum(const Eigen::Matrix4f& vp_matrix) const {
    std::vector<int> result;
    const auto planes = extract_frustum_planes(vp_matrix);
    for (const auto& [hash, cell] : cells_) {
        (void)hash;
        const Eigen::Vector3f aabb_min(cell.ix * cell_size_, cell.iy * cell_size_, cell.iz * cell_size_);
        const Eigen::Vector3f aabb_max = aabb_min + Eigen::Vector3f::Constant(cell_size_);
        if (aabb_intersects_frustum(aabb_min, aabb_max, planes)) {
            result.insert(result.end(), cell.atom_indices.begin(), cell.atom_indices.end());
        }
    }
    return result;
}

int SpatialHash::count() const {
    return static_cast<int>(positions_.size());
}

int64_t SpatialHash::hash_position(float x, float y, float z) const {
    const int64_t ix = static_cast<int64_t>(std::floor(x / cell_size_));
    const int64_t iy = static_cast<int64_t>(std::floor(y / cell_size_));
    const int64_t iz = static_cast<int64_t>(std::floor(z / cell_size_));
    return ix + iy * 10007LL + iz * 100000007LL;
}

LODRenderer::LODRenderer() {
    point_shader_ = std::make_unique<sbox::Shader>("data/shaders/atom_point.vert", "data/shaders/atom_point.frag");

    glGenVertexArrays(1, &point_vao_);
    glGenBuffers(1, &point_vbo_);
    glBindVertexArray(point_vao_);
    glBindBuffer(GL_ARRAY_BUFFER, point_vbo_);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 7 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

LODRenderer::~LODRenderer() {
    if (point_vbo_ != 0) {
        glDeleteBuffers(1, &point_vbo_);
    }
    if (point_vao_ != 0) {
        glDeleteVertexArrays(1, &point_vao_);
    }
}

void LODRenderer::rebuild(const sbox::chem::MolecularSystem& mol,
                          ColorMode color_mode,
                          const sbox::io::PDBData* pdb_data,
                          const std::vector<double>* charges) {
    spatial_hash_.build(mol);
    color_mode_ = color_mode;
    pdb_data_ = pdb_data;
    charges_ = charges;
    full_renderer_.upload(mol, color_mode, pdb_data, charges);
}

void LODRenderer::render(const Eigen::Matrix4f& view_matrix,
                         const Eigen::Matrix4f& proj_matrix,
                         const Eigen::Vector3f& camera_pos,
                         MolRenderMode mode,
                         const sbox::chem::MolecularSystem& mol,
                         ColorMode color_mode,
                         const sbox::io::PDBData* pdb_data,
                         const std::vector<double>* charges) {
    atoms_rendered_ = 0;
    atoms_culled_ = 0;
    bonds_rendered_ = 0;

    if (mol.num_atoms() == 0) {
        return;
    }
    if (spatial_hash_.count() != mol.num_atoms() || color_mode != color_mode_ || pdb_data != pdb_data_ || charges != charges_) {
        rebuild(mol, color_mode, pdb_data, charges);
    }

    const Eigen::Matrix4f vp = proj_matrix * view_matrix;
    std::vector<int> visible = spatial_hash_.query_frustum(vp);
    std::sort(visible.begin(), visible.end());
    visible.erase(std::unique(visible.begin(), visible.end()), visible.end());

    std::unordered_set<int> visible_set(visible.begin(), visible.end());
    const float lod_threshold = mol.num_atoms() < kLodAtomCountThreshold
                                    ? std::numeric_limits<float>::max()
                                    : kLodThreshold;

    std::unordered_set<int> close_set;
    std::vector<float> far_points;
    far_points.reserve(static_cast<std::size_t>(visible.size()) * 7);

    sbox::chem::MolecularSystem close_mol;
    close_mol.set_name(mol.name());
    close_mol.set_charge(mol.charge());
    close_mol.set_multiplicity(mol.multiplicity());
    sbox::io::PDBData close_pdb;
    std::vector<double> close_charges;
    std::vector<int> index_map(static_cast<std::size_t>(mol.num_atoms()), -1);

    const bool has_pdb = pdb_data != nullptr && pdb_data->atoms.size() == mol.atoms().size();
    float b_min = std::numeric_limits<float>::max();
    float b_max = std::numeric_limits<float>::lowest();
    if (has_pdb) {
        for (const auto& atom : pdb_data->atoms) {
            b_min = std::min(b_min, static_cast<float>(atom.b_factor));
            b_max = std::max(b_max, static_cast<float>(atom.b_factor));
        }
        if (b_min == std::numeric_limits<float>::max()) {
            b_min = 0.0f;
            b_max = 1.0f;
        }
    }
    double max_abs_charge = 0.0;
    if (charges != nullptr) {
        for (double q : *charges) {
            max_abs_charge = std::max(max_abs_charge, std::abs(q));
        }
    }

    for (int atom_index : visible) {
        const sbox::chem::Atom& atom = mol.atom(atom_index);
        const Eigen::Vector3f pos = atom.position.cast<float>();
        const float distance = (pos - camera_pos).norm();
        Eigen::Vector3f color = cpk_color(atom.Z);
        if (color_mode == ColorMode::ByChain && has_pdb) {
            int chain_idx = 0;
            for (std::size_t i = 0; i < pdb_data->chains.size(); ++i) {
                bool found = false;
                for (int residue_idx : pdb_data->chains[i].residue_indices) {
                    const auto& residue = pdb_data->residues[static_cast<std::size_t>(residue_idx)];
                    if (std::find(residue.atom_indices.begin(), residue.atom_indices.end(), atom_index) != residue.atom_indices.end()) {
                        chain_idx = static_cast<int>(i);
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            color = chain_color(chain_idx);
        } else if (color_mode == ColorMode::ByResidue && has_pdb) {
            color = residue_color(pdb_data->atoms[static_cast<std::size_t>(atom_index)].residue_name);
        } else if (color_mode == ColorMode::ByBFactor && has_pdb) {
            color = bfactor_color(static_cast<float>(pdb_data->atoms[static_cast<std::size_t>(atom_index)].b_factor), b_min, b_max);
        } else if (color_mode == ColorMode::ByCharge && charges != nullptr && static_cast<std::size_t>(atom_index) < charges->size()) {
            color = charge_color((*charges)[static_cast<std::size_t>(atom_index)], max_abs_charge);
        } else if (color_mode == ColorMode::BySecondary) {
            color = Eigen::Vector3f(0.6f, 0.6f, 0.65f);
        }
        if (distance < lod_threshold) {
            close_set.insert(atom_index);
            index_map[static_cast<std::size_t>(atom_index)] = close_mol.add_atom(atom);
            if (has_pdb) {
                close_pdb.atoms.push_back(pdb_data->atoms[static_cast<std::size_t>(atom_index)]);
            }
            if (charges != nullptr && static_cast<std::size_t>(atom_index) < charges->size()) {
                close_charges.push_back((*charges)[static_cast<std::size_t>(atom_index)]);
            }
        } else {
            far_points.push_back(pos.x());
            far_points.push_back(pos.y());
            far_points.push_back(pos.z());
            far_points.push_back(color.x());
            far_points.push_back(color.y());
            far_points.push_back(color.z());
            far_points.push_back(point_size_for_element(atom.Z));
        }
    }

    for (const sbox::chem::Bond& bond : mol.bonds()) {
        if (close_set.count(bond.atom_i) > 0 && close_set.count(bond.atom_j) > 0) {
            close_mol.add_bond(index_map[static_cast<std::size_t>(bond.atom_i)],
                               index_map[static_cast<std::size_t>(bond.atom_j)],
                               bond.order);
            ++bonds_rendered_;
        }
    }

    atoms_rendered_ = static_cast<int>(visible.size());
    atoms_culled_ = mol.num_atoms() - atoms_rendered_;

    if (close_mol.num_atoms() > 0) {
        full_renderer_.upload(close_mol,
                              color_mode,
                              has_pdb ? &close_pdb : nullptr,
                              charges != nullptr ? &close_charges : nullptr);
        full_renderer_.render(view_matrix, proj_matrix, camera_pos, mode);
    }

    if (!far_points.empty() && point_shader_ != nullptr) {
        glEnable(GL_PROGRAM_POINT_SIZE);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glBindBuffer(GL_ARRAY_BUFFER, point_vbo_);
        glBufferData(GL_ARRAY_BUFFER,
                     static_cast<GLsizeiptr>(far_points.size() * sizeof(float)),
                     far_points.data(),
                     GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        point_shader_->bind();
        point_shader_->setUniform("u_view", view_matrix);
        point_shader_->setUniform("u_proj", proj_matrix);
        const int camera_loc = glGetUniformLocation(point_shader_->id(), "u_camera_pos");
        if (camera_loc >= 0) {
            glUniform3f(camera_loc, camera_pos.x(), camera_pos.y(), camera_pos.z());
        }

        glBindVertexArray(point_vao_);
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(far_points.size() / 7));
        glBindVertexArray(0);
    }
}

int LODRenderer::atoms_rendered() const {
    return atoms_rendered_;
}

int LODRenderer::atoms_culled() const {
    return atoms_culled_;
}

int LODRenderer::bonds_rendered() const {
    return bonds_rendered_;
}

}  // namespace sbox::render
