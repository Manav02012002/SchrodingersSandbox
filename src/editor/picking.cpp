#include "editor/picking.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sbox::editor {

namespace {

float atom_pick_radius(int Z) {
    switch (Z) {
    case 1: return 0.6f;
    case 6: return 1.0f;
    case 7: return 0.95f;
    case 8: return 0.9f;
    case 16: return 1.2f;
    default: return 0.8f;
    }
}

}  // namespace

Ray screen_to_ray(float screen_x,
                  float screen_y,
                  const ImVec2& viewport_pos,
                  const ImVec2& viewport_size,
                  const Eigen::Matrix4f& inv_vp) {
    const float rel_x = screen_x - viewport_pos.x;
    const float rel_y = screen_y - viewport_pos.y;
    const float ndc_x = (rel_x / viewport_size.x) * 2.0f - 1.0f;
    const float ndc_y = 1.0f - (rel_y / viewport_size.y) * 2.0f;

    Eigen::Vector4f near_clip = inv_vp * Eigen::Vector4f(ndc_x, ndc_y, -1.0f, 1.0f);
    Eigen::Vector4f far_clip = inv_vp * Eigen::Vector4f(ndc_x, ndc_y, 1.0f, 1.0f);
    near_clip /= near_clip.w();
    far_clip /= far_clip.w();

    Ray ray;
    ray.origin = near_clip.head<3>();
    ray.direction = (far_clip.head<3>() - near_clip.head<3>()).normalized();
    return ray;
}

PickResult pick_atom(const Ray& ray, const sbox::chem::MolecularSystem& mol, float radius_scale) {
    PickResult best;
    best.distance = std::numeric_limits<float>::infinity();

    for (int i = 0; i < mol.num_atoms(); ++i) {
        const Eigen::Vector3f center = mol.atom(i).position.cast<float>();
        const float radius = atom_pick_radius(mol.atom(i).Z) * radius_scale;
        const Eigen::Vector3f oc = ray.origin - center;
        const float b = oc.dot(ray.direction);
        const float c = oc.dot(oc) - radius * radius;
        const float discriminant = b * b - c;
        if (discriminant < 0.0f) {
            continue;
        }

        const float sqrt_disc = std::sqrt(discriminant);
        float t = -b - sqrt_disc;
        if (t < 0.0f) {
            t = -b + sqrt_disc;
        }
        if (t < 0.0f || t >= best.distance) {
            continue;
        }

        best.type = PickResult::Type::Atom;
        best.index = i;
        best.distance = t;
        best.hit_point = ray.origin + t * ray.direction;
    }

    if (best.type == PickResult::Type::None) {
        best.distance = 0.0f;
    }
    return best;
}

PickResult pick_bond(const Ray& ray, const sbox::chem::MolecularSystem& mol, float radius) {
    PickResult best;
    best.distance = std::numeric_limits<float>::infinity();

    for (int i = 0; i < mol.num_bonds(); ++i) {
        const sbox::chem::Bond& bond = mol.bond(i);
        const Eigen::Vector3f p0 = mol.atom(bond.atom_i).position.cast<float>();
        const Eigen::Vector3f p1 = mol.atom(bond.atom_j).position.cast<float>();
        const Eigen::Vector3f segment = p1 - p0;
        const float segment_len = segment.norm();
        if (segment_len <= 1.0e-6f) {
            continue;
        }

        const Eigen::Vector3f axis = segment / segment_len;
        const Eigen::Vector3f op_full = ray.origin - p0;
        const Eigen::Vector3f dp = ray.direction - ray.direction.dot(axis) * axis;
        const Eigen::Vector3f op = op_full - op_full.dot(axis) * axis;

        const float a = dp.dot(dp);
        if (a <= 1.0e-8f) {
            continue;
        }
        const float b = 2.0f * op.dot(dp);
        const float c = op.dot(op) - radius * radius;
        const float discriminant = b * b - 4.0f * a * c;
        if (discriminant < 0.0f) {
            continue;
        }

        const float sqrt_disc = std::sqrt(discriminant);
        float t = (-b - sqrt_disc) / (2.0f * a);
        if (t < 0.0f) {
            t = (-b + sqrt_disc) / (2.0f * a);
        }
        if (t < 0.0f || t >= best.distance) {
            continue;
        }

        const Eigen::Vector3f hit_point = ray.origin + t * ray.direction;
        const float axis_t = (hit_point - p0).dot(axis);
        if (axis_t < 0.0f || axis_t > segment_len) {
            continue;
        }

        best.type = PickResult::Type::Bond;
        best.index = i;
        best.distance = t;
        best.hit_point = hit_point;
    }

    if (best.type == PickResult::Type::None) {
        best.distance = 0.0f;
    }
    return best;
}

PickResult pick(const Ray& ray, const sbox::chem::MolecularSystem& mol) {
    const PickResult atom_hit = pick_atom(ray, mol);
    const PickResult bond_hit = pick_bond(ray, mol);

    if (atom_hit.type == PickResult::Type::None) {
        return bond_hit;
    }
    if (bond_hit.type == PickResult::Type::None) {
        return atom_hit;
    }
    if (std::abs(atom_hit.distance - bond_hit.distance) <= 0.5f || atom_hit.distance < bond_hit.distance) {
        return atom_hit;
    }
    return bond_hit;
}

bool Selection::has_atom(int i) const {
    return std::find(atoms.begin(), atoms.end(), i) != atoms.end();
}

bool Selection::has_bond(int i) const {
    return std::find(bonds.begin(), bonds.end(), i) != bonds.end();
}

void Selection::toggle_atom(int i) {
    const auto it = std::find(atoms.begin(), atoms.end(), i);
    if (it != atoms.end()) {
        atoms.erase(it);
    } else {
        atoms.push_back(i);
    }
}

void Selection::toggle_bond(int i) {
    const auto it = std::find(bonds.begin(), bonds.end(), i);
    if (it != bonds.end()) {
        bonds.erase(it);
    } else {
        bonds.push_back(i);
    }
}

void Selection::clear() {
    atoms.clear();
    bonds.clear();
}

bool Selection::empty() const {
    return atoms.empty() && bonds.empty();
}

int Selection::num_atoms() const {
    return static_cast<int>(atoms.size());
}

int Selection::num_bonds() const {
    return static_cast<int>(bonds.size());
}

}  // namespace sbox::editor
