#include "editor/measure_mode.h"

#include "core/elements.h"

#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace sbox::editor {

namespace {

constexpr double kBohrToAngstrom = 1.0 / 1.8897259886;
constexpr double kPi = 3.14159265358979323846;

ImVec2 world_to_screen(const Eigen::Vector3f& world_pos,
                       const Eigen::Matrix4f& vp_matrix,
                       const ImVec2& viewport_pos,
                       const ImVec2& viewport_size,
                       bool& visible) {
    const Eigen::Vector4f clip = vp_matrix * Eigen::Vector4f(world_pos.x(), world_pos.y(), world_pos.z(), 1.0f);
    visible = false;
    if (clip.w() <= 0.0f) {
        return {};
    }
    const Eigen::Vector3f ndc = clip.head<3>() / clip.w();
    if (ndc.x() < -1.0f || ndc.x() > 1.0f || ndc.y() < -1.0f || ndc.y() > 1.0f || ndc.z() < -1.0f || ndc.z() > 1.0f) {
        return {};
    }
    visible = true;
    return ImVec2(
        viewport_pos.x + (ndc.x() * 0.5f + 0.5f) * viewport_size.x,
        viewport_pos.y + (1.0f - (ndc.y() * 0.5f + 0.5f)) * viewport_size.y);
}

std::string atom_tag(const sbox::chem::MolecularSystem& mol, int index) {
    return std::string(sbox::elements::get_element(mol.atom(index).Z).symbol) + std::to_string(index + 1);
}

std::string distance_label(const sbox::chem::MolecularSystem& mol, int a, int b, double distance_bohr) {
    char buffer[128];
    std::snprintf(buffer, sizeof(buffer), "d(%s-%s) = %.3f A",
                  atom_tag(mol, a).c_str(),
                  atom_tag(mol, b).c_str(),
                  distance_bohr * kBohrToAngstrom);
    return buffer;
}

std::string angle_label(const sbox::chem::MolecularSystem& mol, int a, int b, int c, double angle_rad) {
    char buffer[160];
    std::snprintf(buffer, sizeof(buffer), "angle(%s-%s-%s) = %.1f deg",
                  atom_tag(mol, a).c_str(),
                  atom_tag(mol, b).c_str(),
                  atom_tag(mol, c).c_str(),
                  angle_rad * 180.0 / std::acos(-1.0));
    return buffer;
}

std::string dihedral_label(const sbox::chem::MolecularSystem& mol, int a, int b, int c, int d, double value_rad) {
    char buffer[192];
    std::snprintf(buffer, sizeof(buffer), "tau(%s-%s-%s-%s) = %.1f deg",
                  atom_tag(mol, a).c_str(),
                  atom_tag(mol, b).c_str(),
                  atom_tag(mol, c).c_str(),
                  atom_tag(mol, d).c_str(),
                  value_rad * 180.0 / std::acos(-1.0));
    return buffer;
}

void draw_centered_text(ImDrawList* draw_list, const ImVec2& pos, ImU32 color, const char* text) {
    const ImVec2 size = ImGui::CalcTextSize(text);
    draw_list->AddText(ImVec2(pos.x - 0.5f * size.x, pos.y - 0.5f * size.y), color, text);
}

}  // namespace

void MeasureMode::on_mouse_down(const Ray& ray, int button, bool shift,
                                sbox::chem::MolecularSystem& mol,
                                Selection& selection,
                                CommandStack& commands) {
    (void)shift;
    (void)selection;
    (void)commands;

    if (button == 1) {
        picked_atoms_.clear();
        return;
    }
    if (button != 0) {
        return;
    }

    const PickResult hit = pick_atom(ray, mol);
    if (hit.type != PickResult::Type::Atom) {
        return;
    }

    picked_atoms_.push_back(hit.index);

    if (picked_atoms_.size() == 2) {
        Measurement m;
        m.type = Measurement::Type::Distance;
        m.atom_indices = picked_atoms_;
        m.value = mol.distance(picked_atoms_[0], picked_atoms_[1]);
        m.label = distance_label(mol, picked_atoms_[0], picked_atoms_[1], m.value);
        measurements_.push_back(std::move(m));
        picked_atoms_.clear();
    } else if (picked_atoms_.size() == 3) {
        Measurement m;
        m.type = Measurement::Type::Angle;
        m.atom_indices = picked_atoms_;
        m.value = mol.angle(picked_atoms_[0], picked_atoms_[1], picked_atoms_[2]);
        m.label = angle_label(mol, picked_atoms_[0], picked_atoms_[1], picked_atoms_[2], m.value);
        measurements_.push_back(std::move(m));
        picked_atoms_.clear();
    } else if (picked_atoms_.size() == 4) {
        Measurement m;
        m.type = Measurement::Type::Dihedral;
        m.atom_indices = picked_atoms_;
        m.value = mol.dihedral(picked_atoms_[0], picked_atoms_[1], picked_atoms_[2], picked_atoms_[3]);
        m.label = dihedral_label(mol, picked_atoms_[0], picked_atoms_[1], picked_atoms_[2], picked_atoms_[3], m.value);
        measurements_.push_back(std::move(m));
        picked_atoms_.clear();
    }
}

void MeasureMode::on_mouse_up(const Ray& ray, int button,
                              sbox::chem::MolecularSystem& mol,
                              Selection& selection,
                              CommandStack& commands) {
    (void)ray;
    (void)button;
    (void)mol;
    (void)selection;
    (void)commands;
}

void MeasureMode::on_mouse_move(const Ray& ray, float dx, float dy, bool dragging,
                                sbox::chem::MolecularSystem& mol,
                                Selection& selection,
                                CommandStack& commands) {
    (void)ray;
    (void)dx;
    (void)dy;
    (void)dragging;
    (void)mol;
    (void)selection;
    (void)commands;
}

void MeasureMode::on_key(int key, bool ctrl, bool shift,
                         sbox::chem::MolecularSystem& mol,
                         Selection& selection,
                         CommandStack& commands) {
    (void)ctrl;
    (void)shift;
    (void)mol;
    (void)selection;
    (void)commands;

    if (key == GLFW_KEY_ESCAPE) {
        picked_atoms_.clear();
    } else if (key == GLFW_KEY_DELETE || key == GLFW_KEY_BACKSPACE) {
        measurements_.clear();
        picked_atoms_.clear();
    }
}

void MeasureMode::clear_measurements() {
    measurements_.clear();
    picked_atoms_.clear();
}

void MeasureMode::draw_overlay(ImDrawList* draw_list,
                               const sbox::chem::MolecularSystem& mol,
                               const Selection& selection,
                               const Eigen::Matrix4f& vp_matrix,
                               const ImVec2& viewport_pos,
                               const ImVec2& viewport_size) {
    (void)selection;
    const ImU32 active = IM_COL32(240, 220, 90, 255);
    const ImU32 persistent = IM_COL32(120, 220, 255, 255);

    for (std::size_t i = 0; i < picked_atoms_.size(); ++i) {
        const int atom_index = picked_atoms_[i];
        if (atom_index < 0 || atom_index >= mol.num_atoms()) {
            continue;
        }
        bool visible = false;
        const ImVec2 p = world_to_screen(mol.atom(atom_index).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible);
        if (!visible) {
            continue;
        }
        draw_list->AddCircle(p, 12.0f, active, 0, 2.5f);
        if (i > 0) {
            bool visible_prev = false;
            const ImVec2 prev = world_to_screen(mol.atom(picked_atoms_[i - 1]).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible_prev);
            if (visible_prev) {
                draw_list->AddLine(prev, p, active, 2.0f);
            }
        }
    }

    for (const Measurement& m : measurements_) {
        if (m.type == Measurement::Type::Distance && m.atom_indices.size() == 2) {
            bool va = false;
            bool vb = false;
            const ImVec2 a = world_to_screen(mol.atom(m.atom_indices[0]).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, va);
            const ImVec2 b = world_to_screen(mol.atom(m.atom_indices[1]).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, vb);
            if (va && vb) {
                draw_list->AddLine(a, b, persistent, 2.0f);
                draw_centered_text(draw_list, ImVec2(0.5f * (a.x + b.x), 0.5f * (a.y + b.y) - 10.0f), persistent, m.label.c_str());
            }
        } else if (m.type == Measurement::Type::Angle && m.atom_indices.size() == 3) {
            const int a_idx = m.atom_indices[0];
            const int b_idx = m.atom_indices[1];
            const int c_idx = m.atom_indices[2];
            bool va = false, vb = false, vc = false;
            const ImVec2 a = world_to_screen(mol.atom(a_idx).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, va);
            const ImVec2 b = world_to_screen(mol.atom(b_idx).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, vb);
            const ImVec2 c = world_to_screen(mol.atom(c_idx).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, vc);
            if (va && vb && vc) {
                draw_list->AddLine(a, b, persistent, 2.0f);
                draw_list->AddLine(b, c, persistent, 2.0f);

                const ImVec2 v1(a.x - b.x, a.y - b.y);
                const ImVec2 v2(c.x - b.x, c.y - b.y);
                const float l1 = std::sqrt(v1.x * v1.x + v1.y * v1.y);
                const float l2 = std::sqrt(v2.x * v2.x + v2.y * v2.y);
                if (l1 > 1.0f && l2 > 1.0f) {
                    const float start_angle = std::atan2(v1.y, v1.x);
                    float end_angle = std::atan2(v2.y, v2.x);
                    float delta = end_angle - start_angle;
                    if (delta > static_cast<float>(kPi)) delta -= static_cast<float>(2.0 * kPi);
                    if (delta < -static_cast<float>(kPi)) delta += static_cast<float>(2.0 * kPi);
                    const float radius = 20.0f;
                    ImVec2 prev(b.x + std::cos(start_angle) * radius, b.y + std::sin(start_angle) * radius);
                    for (int seg = 1; seg <= 20; ++seg) {
                        const float t = static_cast<float>(seg) / 20.0f;
                        const float ang = start_angle + delta * t;
                        const ImVec2 next(b.x + std::cos(ang) * radius, b.y + std::sin(ang) * radius);
                        draw_list->AddLine(prev, next, persistent, 1.5f);
                        prev = next;
                    }
                    draw_centered_text(draw_list, ImVec2(b.x, b.y - radius - 14.0f), persistent, m.label.c_str());
                }
            }
        } else if (m.type == Measurement::Type::Dihedral && m.atom_indices.size() == 4) {
            bool all_visible = true;
            ImVec2 pts[4];
            for (int i = 0; i < 4; ++i) {
                bool visible = false;
                pts[i] = world_to_screen(mol.atom(m.atom_indices[i]).position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible);
                all_visible = all_visible && visible;
            }
            if (all_visible) {
                draw_list->AddLine(pts[0], pts[1], persistent, 2.0f);
                draw_list->AddLine(pts[1], pts[2], persistent, 2.0f);
                draw_list->AddLine(pts[2], pts[3], persistent, 2.0f);
                draw_centered_text(draw_list, ImVec2(0.5f * (pts[1].x + pts[2].x), 0.5f * (pts[1].y + pts[2].y) - 12.0f), persistent, m.label.c_str());
            }
        }
    }
}

}  // namespace sbox::editor
