#pragma once

#include "core/molecular_system.h"

#include <Eigen/Core>

#include <imgui.h>

#include <vector>

namespace sbox::editor {

struct Ray {
    Eigen::Vector3f origin;
    Eigen::Vector3f direction;
};

struct PickResult {
    enum class Type { None, Atom, Bond };
    Type type = Type::None;
    int index = -1;
    float distance = 0.0f;
    Eigen::Vector3f hit_point = Eigen::Vector3f::Zero();
};

Ray screen_to_ray(float screen_x,
                  float screen_y,
                  const ImVec2& viewport_pos,
                  const ImVec2& viewport_size,
                  const Eigen::Matrix4f& inv_vp);

PickResult pick_atom(const Ray& ray,
                     const sbox::chem::MolecularSystem& mol,
                     float radius_scale = 1.0f);

PickResult pick_bond(const Ray& ray,
                     const sbox::chem::MolecularSystem& mol,
                     float radius = 0.3f);

PickResult pick(const Ray& ray, const sbox::chem::MolecularSystem& mol);

struct Selection {
    std::vector<int> atoms;
    std::vector<int> bonds;

    bool has_atom(int i) const;
    bool has_bond(int i) const;
    void toggle_atom(int i);
    void toggle_bond(int i);
    void clear();
    bool empty() const;
    int num_atoms() const;
    int num_bonds() const;
};

}  // namespace sbox::editor
