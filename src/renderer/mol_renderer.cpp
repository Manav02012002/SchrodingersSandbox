#include "renderer/mol_renderer.h"

#include "core/elements.h"

#include <glad/gl.h>

#include <array>
#include <cstddef>
#include <cstdio>
#include <stdexcept>

namespace sbox::render {

namespace {

Eigen::Vector3f cpk_color(int Z) {
    switch (Z) {
    case 1:
        return Eigen::Vector3f(1.0f, 1.0f, 1.0f);
    case 2:
        return Eigen::Vector3f(0.85f, 1.0f, 1.0f);
    case 3:
        return Eigen::Vector3f(0.8f, 0.5f, 1.0f);
    case 5:
        return Eigen::Vector3f(1.0f, 0.71f, 0.71f);
    case 6:
        return Eigen::Vector3f(0.56f, 0.56f, 0.56f);
    case 7:
        return Eigen::Vector3f(0.19f, 0.31f, 0.97f);
    case 8:
        return Eigen::Vector3f(1.0f, 0.05f, 0.05f);
    case 9:
        return Eigen::Vector3f(0.56f, 0.88f, 0.31f);
    case 10:
        return Eigen::Vector3f(0.70f, 0.89f, 0.96f);
    case 11:
        return Eigen::Vector3f(0.67f, 0.36f, 0.95f);
    case 12:
        return Eigen::Vector3f(0.54f, 1.0f, 0.0f);
    case 13:
        return Eigen::Vector3f(0.75f, 0.65f, 0.65f);
    case 14:
        return Eigen::Vector3f(0.94f, 0.78f, 0.63f);
    case 15:
        return Eigen::Vector3f(1.0f, 0.50f, 0.0f);
    case 16:
        return Eigen::Vector3f(1.0f, 1.0f, 0.19f);
    case 17:
        return Eigen::Vector3f(0.12f, 0.94f, 0.12f);
    case 18:
        return Eigen::Vector3f(0.50f, 0.82f, 0.89f);
    case 19:
        return Eigen::Vector3f(0.56f, 0.25f, 0.83f);
    case 20:
        return Eigen::Vector3f(0.24f, 1.0f, 0.0f);
    case 26:
        return Eigen::Vector3f(0.88f, 0.40f, 0.20f);
    case 29:
        return Eigen::Vector3f(0.78f, 0.50f, 0.20f);
    case 30:
        return Eigen::Vector3f(0.49f, 0.50f, 0.69f);
    case 35:
        return Eigen::Vector3f(0.65f, 0.16f, 0.16f);
    case 53:
        return Eigen::Vector3f(0.58f, 0.0f, 0.58f);
    case 79:
        return Eigen::Vector3f(1.0f, 0.82f, 0.14f);
    default:
        return Eigen::Vector3f(0.5f, 0.5f, 0.5f);
    }
}

float atom_render_radius(int Z) {
    switch (Z) {
    case 1:
        return 0.6f;
    case 6:
        return 1.0f;
    case 7:
        return 0.95f;
    case 8:
        return 0.9f;
    case 16:
        return 1.2f;
    default:
        return 0.8f;
    }
}

float vdw_radius(int Z) {
    switch (Z) {
    case 1:
        return 2.27f;
    case 2:
        return 2.65f;
    case 6:
        return 3.21f;
    case 7:
        return 2.91f;
    case 8:
        return 2.87f;
    case 9:
        return 2.78f;
    case 16:
        return 3.40f;
    case 17:
        return 3.31f;
    default:
        return 3.0f;
    }
}

constexpr float kBondRadiusBallAndStick = 0.15f;
constexpr float kBondRadiusStickOnly = 0.08f;

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

void create_quad_buffer(unsigned int* vao,
                        unsigned int* vertex_buffer,
                        unsigned int* instance_buffer,
                        const std::array<float, 12>& vertices) {
    glGenVertexArrays(1, vao);
    glGenBuffers(1, vertex_buffer);
    glGenBuffers(1, instance_buffer);

    glBindVertexArray(*vao);
    glBindBuffer(GL_ARRAY_BUFFER, *vertex_buffer);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(vertices.size() * sizeof(float)),
                 vertices.data(),
                 GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), nullptr);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void upload_atom_instances(unsigned int buffer, const std::vector<float>& data) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(data.size() * sizeof(float)),
                 data.empty() ? nullptr : data.data(),
                 GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void upload_bond_instances(unsigned int buffer, const std::vector<float>& data) {
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(data.size() * sizeof(float)),
                 data.empty() ? nullptr : data.data(),
                 GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

}  // namespace

MolRenderer::MolRenderer() {
    atom_shader_ = std::make_unique<sbox::Shader>("data/shaders/atom_impostor.vert", "data/shaders/atom_impostor.frag");
    bond_shader_ = std::make_unique<sbox::Shader>("data/shaders/bond_impostor.vert", "data/shaders/bond_impostor.frag");

    constexpr std::array<float, 12> quad_vertices = {
        -1.0f, -1.0f,
         1.0f, -1.0f,
        -1.0f,  1.0f,
         1.0f,  1.0f,
        -1.0f,  0.0f,
         1.0f,  0.0f,
    };

    create_quad_buffer(&atom_vao_, &atom_vbo_, &atom_instance_vbo_, quad_vertices);
    create_quad_buffer(&bond_vao_, &bond_vbo_, &bond_instance_vbo_, quad_vertices);

    glBindVertexArray(atom_vao_);
    glBindBuffer(GL_ARRAY_BUFFER, atom_instance_vbo_);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 8 * sizeof(float), nullptr);
    glVertexAttribDivisor(1, 1);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2,
                          4,
                          GL_FLOAT,
                          GL_FALSE,
                          8 * sizeof(float),
                          reinterpret_cast<void*>(4 * sizeof(float)));
    glVertexAttribDivisor(2, 1);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    glBindVertexArray(bond_vao_);
    glBindBuffer(GL_ARRAY_BUFFER, bond_instance_vbo_);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 13 * sizeof(float), nullptr);
    glVertexAttribDivisor(1, 1);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2,
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          13 * sizeof(float),
                          reinterpret_cast<void*>(3 * sizeof(float)));
    glVertexAttribDivisor(2, 1);
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3,
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          13 * sizeof(float),
                          reinterpret_cast<void*>(6 * sizeof(float)));
    glVertexAttribDivisor(3, 1);
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4,
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          13 * sizeof(float),
                          reinterpret_cast<void*>(9 * sizeof(float)));
    glVertexAttribDivisor(4, 1);
    glEnableVertexAttribArray(5);
    glVertexAttribPointer(5,
                          1,
                          GL_FLOAT,
                          GL_FALSE,
                          13 * sizeof(float),
                          reinterpret_cast<void*>(12 * sizeof(float)));
    glVertexAttribDivisor(5, 1);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

MolRenderer::~MolRenderer() {
    if (atom_instance_vbo_ != 0) {
        glDeleteBuffers(1, &atom_instance_vbo_);
    }
    if (atom_vbo_ != 0) {
        glDeleteBuffers(1, &atom_vbo_);
    }
    if (atom_vao_ != 0) {
        glDeleteVertexArrays(1, &atom_vao_);
    }

    if (bond_instance_vbo_ != 0) {
        glDeleteBuffers(1, &bond_instance_vbo_);
    }
    if (bond_vbo_ != 0) {
        glDeleteBuffers(1, &bond_vbo_);
    }
    if (bond_vao_ != 0) {
        glDeleteVertexArrays(1, &bond_vao_);
    }
}

void MolRenderer::upload(const sbox::chem::MolecularSystem& mol) {
    atom_instances_.clear();
    bond_instances_.clear();

    atom_count_ = mol.num_atoms();
    bond_count_ = mol.num_bonds();

    atom_instances_.reserve(static_cast<std::size_t>(atom_count_) * 8);
    for (const sbox::chem::Atom& atom : mol.atoms()) {
        const Eigen::Vector3f color = cpk_color(atom.Z);
        atom_instances_.push_back(static_cast<float>(atom.position.x()));
        atom_instances_.push_back(static_cast<float>(atom.position.y()));
        atom_instances_.push_back(static_cast<float>(atom.position.z()));
        atom_instances_.push_back(atom_render_radius(atom.Z));
        atom_instances_.push_back(color.x());
        atom_instances_.push_back(color.y());
        atom_instances_.push_back(color.z());
        atom_instances_.push_back(static_cast<float>(atom.Z));
    }

    bond_instances_.reserve(static_cast<std::size_t>(bond_count_) * 13);
    for (const sbox::chem::Bond& bond : mol.bonds()) {
        const sbox::chem::Atom& atom_a = mol.atom(bond.atom_i);
        const sbox::chem::Atom& atom_b = mol.atom(bond.atom_j);
        const Eigen::Vector3f color_a = cpk_color(atom_a.Z);
        const Eigen::Vector3f color_b = cpk_color(atom_b.Z);

        bond_instances_.push_back(static_cast<float>(atom_a.position.x()));
        bond_instances_.push_back(static_cast<float>(atom_a.position.y()));
        bond_instances_.push_back(static_cast<float>(atom_a.position.z()));
        bond_instances_.push_back(static_cast<float>(atom_b.position.x()));
        bond_instances_.push_back(static_cast<float>(atom_b.position.y()));
        bond_instances_.push_back(static_cast<float>(atom_b.position.z()));
        bond_instances_.push_back(color_a.x());
        bond_instances_.push_back(color_a.y());
        bond_instances_.push_back(color_a.z());
        bond_instances_.push_back(color_b.x());
        bond_instances_.push_back(color_b.y());
        bond_instances_.push_back(color_b.z());
        bond_instances_.push_back(kBondRadiusBallAndStick);
    }

    upload_atom_instances(atom_instance_vbo_, atom_instances_);
    upload_bond_instances(bond_instance_vbo_, bond_instances_);
}

void MolRenderer::render(const Eigen::Matrix4f& view_matrix,
                         const Eigen::Matrix4f& proj_matrix,
                         const Eigen::Vector3f& camera_pos,
                         MolRenderMode mode) {
    if (!has_data()) {
        return;
    }

    glEnable(GL_DEPTH_TEST);

    const bool render_atoms = (mode == MolRenderMode::BallAndStick || mode == MolRenderMode::SpaceFilling);
    const bool render_bonds = (mode == MolRenderMode::BallAndStick || mode == MolRenderMode::StickOnly || mode == MolRenderMode::Wireframe);

    if (render_atoms) {
        std::vector<float> atom_draw_data = atom_instances_;
        if (mode == MolRenderMode::SpaceFilling) {
            for (int i = 0; i < atom_count_; ++i) {
                const std::size_t base = static_cast<std::size_t>(i) * 8;
                atom_draw_data[base + 3] = vdw_radius(static_cast<int>(atom_draw_data[base + 7]));
            }
        }

        upload_atom_instances(atom_instance_vbo_, atom_draw_data);

        atom_shader_->bind();
        atom_shader_->setUniform("u_view", view_matrix);
        atom_shader_->setUniform("u_proj", proj_matrix);
        atom_shader_->setUniform("u_camera_pos", camera_pos);

        glBindVertexArray(atom_vao_);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, atom_count_);
        glBindVertexArray(0);
    }

    if (render_bonds && bond_count_ > 0) {
        std::vector<float> bond_draw_data = bond_instances_;
        const float bond_radius = (mode == MolRenderMode::StickOnly) ? kBondRadiusStickOnly : kBondRadiusBallAndStick;
        for (int i = 0; i < bond_count_; ++i) {
            bond_draw_data[static_cast<std::size_t>(i) * 13 + 12] = bond_radius;
        }
        upload_bond_instances(bond_instance_vbo_, bond_draw_data);

        bond_shader_->bind();
        bond_shader_->setUniform("u_view", view_matrix);
        bond_shader_->setUniform("u_proj", proj_matrix);
        bond_shader_->setUniform("u_camera_pos", camera_pos);
        bond_shader_->setUniform("u_wireframe", mode == MolRenderMode::Wireframe ? 1 : 0);

        glBindVertexArray(bond_vao_);
        if (mode == MolRenderMode::Wireframe) {
            glLineWidth(1.5f);
            glDrawArraysInstanced(GL_LINES, 4, 2, bond_count_);
        } else {
            glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, bond_count_);
        }
        glBindVertexArray(0);
    }
}

void MolRenderer::render_highlights(const Eigen::Matrix4f& view_matrix,
                                    const Eigen::Matrix4f& proj_matrix,
                                    const Eigen::Vector3f& camera_pos,
                                    const std::vector<float>& atom_data,
                                    const std::vector<float>& bond_data,
                                    const Eigen::Vector3f& highlight_color,
                                    MolRenderMode mode) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glDepthMask(GL_FALSE);
    glEnable(GL_DEPTH_TEST);

    if (!atom_data.empty() && mode != MolRenderMode::StickOnly && mode != MolRenderMode::Wireframe) {
        std::vector<float> draw_data = atom_data;
        for (std::size_t i = 0; i + 7 < draw_data.size(); i += 8) {
            draw_data[i + 3] *= 1.3f;
            draw_data[i + 4] = highlight_color.x();
            draw_data[i + 5] = highlight_color.y();
            draw_data[i + 6] = highlight_color.z();
        }
        upload_atom_instances(atom_instance_vbo_, draw_data);

        atom_shader_->bind();
        atom_shader_->setUniform("u_view", view_matrix);
        atom_shader_->setUniform("u_proj", proj_matrix);
        atom_shader_->setUniform("u_camera_pos", camera_pos);
        glBindVertexArray(atom_vao_);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, static_cast<int>(draw_data.size() / 8));
        glBindVertexArray(0);
    }

    if (!bond_data.empty() && mode != MolRenderMode::SpaceFilling) {
        std::vector<float> draw_data = bond_data;
        for (std::size_t i = 0; i + 12 < draw_data.size(); i += 13) {
            draw_data[i + 6] = highlight_color.x();
            draw_data[i + 7] = highlight_color.y();
            draw_data[i + 8] = highlight_color.z();
            draw_data[i + 9] = highlight_color.x();
            draw_data[i + 10] = highlight_color.y();
            draw_data[i + 11] = highlight_color.z();
            draw_data[i + 12] *= 1.5f;
        }
        upload_bond_instances(bond_instance_vbo_, draw_data);

        bond_shader_->bind();
        bond_shader_->setUniform("u_view", view_matrix);
        bond_shader_->setUniform("u_proj", proj_matrix);
        bond_shader_->setUniform("u_camera_pos", camera_pos);
        bond_shader_->setUniform("u_wireframe", 0);
        glBindVertexArray(bond_vao_);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, static_cast<int>(draw_data.size() / 13));
        glBindVertexArray(0);
    }

    upload_atom_instances(atom_instance_vbo_, atom_instances_);
    upload_bond_instances(bond_instance_vbo_, bond_instances_);
    glDepthMask(GL_TRUE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_BLEND);
}

void MolRenderer::render_selection(const Eigen::Matrix4f& view_matrix,
                                   const Eigen::Matrix4f& proj_matrix,
                                   const Eigen::Vector3f& camera_pos,
                                   const sbox::chem::MolecularSystem& mol,
                                   const sbox::editor::Selection& selection,
                                   MolRenderMode mode) {
    if (selection.empty() || !has_data()) {
        return;
    }

    std::vector<float> atom_data;
    atom_data.reserve(static_cast<std::size_t>(selection.atoms.size()) * 8);
    for (int atom_index : selection.atoms) {
        if (atom_index < 0 || atom_index >= mol.num_atoms()) {
            continue;
        }
        const sbox::chem::Atom& atom = mol.atom(atom_index);
        atom_data.push_back(static_cast<float>(atom.position.x()));
        atom_data.push_back(static_cast<float>(atom.position.y()));
        atom_data.push_back(static_cast<float>(atom.position.z()));
        atom_data.push_back(mode == MolRenderMode::SpaceFilling ? vdw_radius(atom.Z) : atom_render_radius(atom.Z));
        atom_data.push_back(0.15f);
        atom_data.push_back(0.55f);
        atom_data.push_back(0.65f);
        atom_data.push_back(static_cast<float>(atom.Z));
    }

    std::vector<float> bond_data;
    bond_data.reserve(static_cast<std::size_t>(selection.bonds.size()) * 13);
    for (int bond_index : selection.bonds) {
        if (bond_index < 0 || bond_index >= mol.num_bonds()) {
            continue;
        }
        const sbox::chem::Bond& bond = mol.bond(bond_index);
        const sbox::chem::Atom& atom_a = mol.atom(bond.atom_i);
        const sbox::chem::Atom& atom_b = mol.atom(bond.atom_j);
        bond_data.push_back(static_cast<float>(atom_a.position.x()));
        bond_data.push_back(static_cast<float>(atom_a.position.y()));
        bond_data.push_back(static_cast<float>(atom_a.position.z()));
        bond_data.push_back(static_cast<float>(atom_b.position.x()));
        bond_data.push_back(static_cast<float>(atom_b.position.y()));
        bond_data.push_back(static_cast<float>(atom_b.position.z()));
        bond_data.push_back(0.15f);
        bond_data.push_back(0.55f);
        bond_data.push_back(0.65f);
        bond_data.push_back(0.15f);
        bond_data.push_back(0.55f);
        bond_data.push_back(0.65f);
        bond_data.push_back(mode == MolRenderMode::StickOnly ? kBondRadiusStickOnly : kBondRadiusBallAndStick);
    }

    render_highlights(view_matrix, proj_matrix, camera_pos, atom_data, bond_data, Eigen::Vector3f(0.15f, 0.55f, 0.65f), mode);
}

void MolRenderer::render_atom_labels(const sbox::chem::MolecularSystem& mol,
                                     const Eigen::Matrix4f& vp_matrix,
                                     const ImVec2& viewport_pos,
                                     const ImVec2& viewport_size,
                                     bool show_indices,
                                     bool show_symbols) {
    if (!show_symbols && !show_indices) {
        return;
    }

    ImDrawList* draw_list = ImGui::GetForegroundDrawList();
    for (int i = 0; i < mol.num_atoms(); ++i) {
        const sbox::chem::Atom& atom = mol.atom(i);
        if (atom.Z == 1) {
            continue;
        }
        bool visible = false;
        const ImVec2 p = world_to_screen(atom.position.cast<float>(), vp_matrix, viewport_pos, viewport_size, visible);
        if (!visible) {
            continue;
        }

        char label[32];
        if (show_symbols && show_indices) {
            std::snprintf(label, sizeof(label), "%s%d", sbox::elements::get_element(atom.Z).symbol, i + 1);
        } else if (show_symbols) {
            std::snprintf(label, sizeof(label), "%s", sbox::elements::get_element(atom.Z).symbol);
        } else {
            std::snprintf(label, sizeof(label), "%d", i + 1);
        }
        const ImVec2 text_size = ImGui::CalcTextSize(label);
        draw_list->AddText(ImVec2(p.x - 0.5f * text_size.x, p.y - 0.5f * text_size.y), IM_COL32(245, 248, 255, 220), label);
    }
}

bool MolRenderer::has_data() const {
    return atom_count_ > 0 || bond_count_ > 0;
}

int MolRenderer::num_atoms() const {
    return atom_count_;
}

int MolRenderer::num_bonds() const {
    return bond_count_;
}

}  // namespace sbox::render
