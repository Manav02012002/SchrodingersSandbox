#include "renderer/shadow_map.h"

#include "core/paths.h"

#include <glad/gl.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>

namespace sbox::render {

namespace {

float atom_shadow_radius(int atomic_number) {
    switch (atomic_number) {
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

Eigen::Matrix4f look_at(const Eigen::Vector3f& eye,
                        const Eigen::Vector3f& center,
                        const Eigen::Vector3f& up) {
    const Eigen::Vector3f f = (center - eye).normalized();
    const Eigen::Vector3f s = f.cross(up).normalized();
    const Eigen::Vector3f u = s.cross(f);

    Eigen::Matrix4f mat = Eigen::Matrix4f::Identity();
    mat(0, 0) = s.x();
    mat(0, 1) = s.y();
    mat(0, 2) = s.z();
    mat(1, 0) = u.x();
    mat(1, 1) = u.y();
    mat(1, 2) = u.z();
    mat(2, 0) = -f.x();
    mat(2, 1) = -f.y();
    mat(2, 2) = -f.z();
    mat(0, 3) = -s.dot(eye);
    mat(1, 3) = -u.dot(eye);
    mat(2, 3) = f.dot(eye);
    return mat;
}

Eigen::Matrix4f orthographic(float left,
                             float right,
                             float bottom,
                             float top,
                             float z_near,
                             float z_far) {
    Eigen::Matrix4f mat = Eigen::Matrix4f::Identity();
    mat(0, 0) = 2.0f / (right - left);
    mat(1, 1) = 2.0f / (top - bottom);
    mat(2, 2) = -2.0f / (z_far - z_near);
    mat(0, 3) = -(right + left) / (right - left);
    mat(1, 3) = -(top + bottom) / (top - bottom);
    mat(2, 3) = -(z_far + z_near) / (z_far - z_near);
    return mat;
}

struct EdgeKey {
    uint32_t a = 0;
    uint32_t b = 0;

    bool operator<(const EdgeKey& other) const {
        if (a != other.a) {
            return a < other.a;
        }
        return b < other.b;
    }
};

uint32_t midpoint_index(uint32_t a,
                        uint32_t b,
                        std::vector<Eigen::Vector3f>& vertices,
                        std::map<EdgeKey, uint32_t>& cache) {
    const EdgeKey key{std::min(a, b), std::max(a, b)};
    const auto it = cache.find(key);
    if (it != cache.end()) {
        return it->second;
    }

    const Eigen::Vector3f midpoint = (vertices[a] + vertices[b]).normalized();
    vertices.push_back(midpoint);
    const uint32_t index = static_cast<uint32_t>(vertices.size() - 1);
    cache.emplace(key, index);
    return index;
}

}  // namespace

ShadowMap::ShadowMap() = default;

ShadowMap::~ShadowMap() {
    destroy();
}

void ShadowMap::init(int resolution) {
    resolution_ = std::max(resolution, 256);

    if (!depth_shader_) {
        depth_shader_ = std::make_unique<Shader>(sbox::get_shader_path("shadow_depth.vert"),
                                                 sbox::get_shader_path("shadow_depth.frag"));
    }
    if (sphere_vao_ == 0U) {
        create_sphere_mesh();
    }
    if (shadow_fbo_ == 0U) {
        glGenFramebuffers(1, &shadow_fbo_);
    }
    if (shadow_tex_ == 0U) {
        glGenTextures(1, &shadow_tex_);
    }

    glBindTexture(GL_TEXTURE_2D, shadow_tex_);
    glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_DEPTH_COMPONENT24,
                 resolution_,
                 resolution_,
                 0,
                 GL_DEPTH_COMPONENT,
                 GL_FLOAT,
                 nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    const std::array<float, 4> border_color = {1.0f, 1.0f, 1.0f, 1.0f};
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, border_color.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);

    glBindFramebuffer(GL_FRAMEBUFFER, shadow_fbo_);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, shadow_tex_, 0);
    glDrawBuffer(GL_NONE);
    glReadBuffer(GL_NONE);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        throw std::runtime_error("Shadow map framebuffer is incomplete");
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
    initialized_ = true;
}

void ShadowMap::compute(MolRenderer& mol_renderer,
                        const sbox::chem::MolecularSystem& mol,
                        const Eigen::Vector3f& scene_center,
                        float scene_radius,
                        const Eigen::Vector3f& light_dir) {
    if (!initialized_) {
        init(resolution_);
    }
    if (!mol_renderer.has_data() || mol.num_atoms() == 0) {
        return;
    }

    const float fit_radius = std::max(scene_radius, 1.0f);
    const Eigen::Vector3f light_direction = light_dir.normalized();
    const Eigen::Vector3f light_position = scene_center - light_direction * fit_radius * 2.0f;
    const Eigen::Vector3f up = std::abs(light_direction.dot(Eigen::Vector3f::UnitY())) > 0.95f
                                   ? Eigen::Vector3f::UnitZ()
                                   : Eigen::Vector3f::UnitY();
    const Eigen::Matrix4f view = look_at(light_position, scene_center, up);
    const Eigen::Matrix4f projection = orthographic(-fit_radius,
                                                    fit_radius,
                                                    -fit_radius,
                                                    fit_radius,
                                                    0.1f,
                                                    fit_radius * 4.0f);
    light_vp_ = projection * view;

    std::vector<float> instances;
    instances.reserve(static_cast<std::size_t>(mol.num_atoms()) * 4);
    for (const sbox::chem::Atom& atom : mol.atoms()) {
        instances.push_back(static_cast<float>(atom.position.x()));
        instances.push_back(static_cast<float>(atom.position.y()));
        instances.push_back(static_cast<float>(atom.position.z()));
        instances.push_back(atom_shadow_radius(atom.Z));
    }

    glBindBuffer(GL_ARRAY_BUFFER, instance_vbo_);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(instances.size() * sizeof(float)),
                 instances.data(),
                 GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, shadow_fbo_);
    glViewport(0, 0, resolution_, resolution_);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glClear(GL_DEPTH_BUFFER_BIT);

    depth_shader_->bind();
    depth_shader_->setUniform("u_light_vp", light_vp_);

    glBindVertexArray(sphere_vao_);
    glDrawElementsInstanced(GL_TRIANGLES,
                            sphere_index_count_,
                            GL_UNSIGNED_INT,
                            nullptr,
                            mol.num_atoms());
    glBindVertexArray(0);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void ShadowMap::bind_shadow_texture(int unit) const {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, shadow_tex_);
}

Eigen::Matrix4f ShadowMap::light_vp_matrix() const {
    return light_vp_;
}

bool ShadowMap::is_initialized() const {
    return initialized_;
}

int ShadowMap::resolution() const {
    return resolution_;
}

void ShadowMap::destroy() {
    initialized_ = false;
    depth_shader_.reset();

    if (instance_vbo_ != 0U) {
        glDeleteBuffers(1, &instance_vbo_);
        instance_vbo_ = 0;
    }
    if (sphere_ebo_ != 0U) {
        glDeleteBuffers(1, &sphere_ebo_);
        sphere_ebo_ = 0;
    }
    if (sphere_vbo_ != 0U) {
        glDeleteBuffers(1, &sphere_vbo_);
        sphere_vbo_ = 0;
    }
    if (sphere_vao_ != 0U) {
        glDeleteVertexArrays(1, &sphere_vao_);
        sphere_vao_ = 0;
    }
    if (shadow_tex_ != 0U) {
        glDeleteTextures(1, &shadow_tex_);
        shadow_tex_ = 0;
    }
    if (shadow_fbo_ != 0U) {
        glDeleteFramebuffers(1, &shadow_fbo_);
        shadow_fbo_ = 0;
    }
}

void ShadowMap::create_sphere_mesh() {
    const float t = (1.0f + std::sqrt(5.0f)) * 0.5f;
    std::vector<Eigen::Vector3f> vertices = {
        {-1.0f, t, 0.0f},  {1.0f, t, 0.0f},   {-1.0f, -t, 0.0f}, {1.0f, -t, 0.0f},
        {0.0f, -1.0f, t},  {0.0f, 1.0f, t},   {0.0f, -1.0f, -t}, {0.0f, 1.0f, -t},
        {t, 0.0f, -1.0f},  {t, 0.0f, 1.0f},   {-t, 0.0f, -1.0f}, {-t, 0.0f, 1.0f},
    };
    for (Eigen::Vector3f& vertex : vertices) {
        vertex.normalize();
    }

    std::vector<uint32_t> indices = {
        0, 11, 5,  0, 5,  1,  0, 1,  7,  0, 7, 10, 0, 10, 11,
        1, 5,  9,  5, 11, 4, 11, 10, 2, 10, 7,  6, 7, 1,  8,
        3, 9,  4,  3, 4,  2,  3, 2,  6,  3, 6, 8,  3, 8,  9,
        4, 9,  5,  2, 4, 11, 6, 2, 10, 8, 6,  7, 9, 8,  1,
    };

    std::map<EdgeKey, uint32_t> midpoint_cache;
    std::vector<uint32_t> subdivided;
    subdivided.reserve(indices.size() * 4);
    for (std::size_t i = 0; i + 2 < indices.size(); i += 3) {
        const uint32_t i0 = indices[i + 0];
        const uint32_t i1 = indices[i + 1];
        const uint32_t i2 = indices[i + 2];
        const uint32_t a = midpoint_index(i0, i1, vertices, midpoint_cache);
        const uint32_t b = midpoint_index(i1, i2, vertices, midpoint_cache);
        const uint32_t c = midpoint_index(i2, i0, vertices, midpoint_cache);

        subdivided.insert(subdivided.end(), {i0, a, c, i1, b, a, i2, c, b, a, b, c});
    }
    indices = std::move(subdivided);
    sphere_index_count_ = static_cast<int>(indices.size());

    glGenVertexArrays(1, &sphere_vao_);
    glGenBuffers(1, &sphere_vbo_);
    glGenBuffers(1, &sphere_ebo_);
    glGenBuffers(1, &instance_vbo_);

    glBindVertexArray(sphere_vao_);

    glBindBuffer(GL_ARRAY_BUFFER, sphere_vbo_);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(vertices.size() * sizeof(Eigen::Vector3f)),
                 vertices.data(),
                 GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Eigen::Vector3f), nullptr);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphere_ebo_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(indices.size() * sizeof(uint32_t)),
                 indices.data(),
                 GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, instance_vbo_);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), nullptr);
    glVertexAttribDivisor(1, 1);

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

}  // namespace sbox::render
