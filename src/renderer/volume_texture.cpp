#include "renderer/volume_texture.h"

#include <glad/gl.h>

#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace sbox::render {

namespace {

void set_sampler_uniform(unsigned int shader_id, const char* name, int value) {
    const int location = glGetUniformLocation(shader_id, name);
    if (location >= 0) {
        glUniform1i(location, value);
    }
}

}  // namespace

VolumeTexture::VolumeTexture() {
    glGenTextures(1, &texture_3d_);
    glBindTexture(GL_TEXTURE_3D, texture_3d_);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);
    const float border[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    glTexParameterfv(GL_TEXTURE_3D, GL_TEXTURE_BORDER_COLOR, border);
    glBindTexture(GL_TEXTURE_3D, 0);
}

VolumeTexture::~VolumeTexture() {
    if (texture_3d_ != 0) {
        glDeleteTextures(1, &texture_3d_);
    }
}

bool VolumeTexture::upload(const sbox::io::CubeData& cube) {
    return upload(cube.origin, cube.step_x, cube.step_y, cube.step_z, cube.nx, cube.ny, cube.nz, cube.data.data());
}

bool VolumeTexture::upload(const Eigen::Vector3d& origin,
                           const Eigen::Vector3d& step_x,
                           const Eigen::Vector3d& step_y,
                           const Eigen::Vector3d& step_z,
                           int nx,
                           int ny,
                           int nz,
                           const float* data) {
    if (nx <= 0 || ny <= 0 || nz <= 0 || data == nullptr) {
        uploaded_ = false;
        return false;
    }

    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
    origin_ = origin.cast<float>();
    grid_to_world_.col(0) = step_x.cast<float>();
    grid_to_world_.col(1) = step_y.cast<float>();
    grid_to_world_.col(2) = step_z.cast<float>();

    const float det = grid_to_world_.determinant();
    if (std::abs(det) < 1e-8f) {
        throw std::runtime_error("Volume grid transform is singular");
    }
    world_to_grid_ = grid_to_world_.inverse();

    max_abs_val_ = 0.0f;
    const std::size_t num_values =
        static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_) * static_cast<std::size_t>(nz_);
    for (std::size_t i = 0; i < num_values; ++i) {
        max_abs_val_ = std::max(max_abs_val_, std::abs(data[i]));
    }
    if (max_abs_val_ <= 0.0f) {
        max_abs_val_ = 1.0f;
    }

    glBindTexture(GL_TEXTURE_3D, texture_3d_);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);
    const float border[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    glTexParameterfv(GL_TEXTURE_3D, GL_TEXTURE_BORDER_COLOR, border);
    glTexImage3D(GL_TEXTURE_3D,
                 0,
                 GL_R32F,
                 nx_,
                 ny_,
                 nz_,
                 0,
                 GL_RED,
                 GL_FLOAT,
                 data);
    glBindTexture(GL_TEXTURE_3D, 0);

    uploaded_ = true;
    return true;
}

void VolumeTexture::bind(unsigned int shader_id, int texture_unit) const {
    bound_texture_unit_ = texture_unit;
    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(texture_unit));
    glBindTexture(GL_TEXTURE_3D, texture_3d_);
    set_sampler_uniform(shader_id, "u_volume", texture_unit);
}

void VolumeTexture::unbind() const {
    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(bound_texture_unit_));
    glBindTexture(GL_TEXTURE_3D, 0);
}

bool VolumeTexture::is_uploaded() const {
    return uploaded_;
}

int VolumeTexture::nx() const {
    return nx_;
}

int VolumeTexture::ny() const {
    return ny_;
}

int VolumeTexture::nz() const {
    return nz_;
}

Eigen::Vector3f VolumeTexture::origin() const {
    return origin_;
}

Eigen::Matrix3f VolumeTexture::grid_to_world() const {
    return grid_to_world_;
}

Eigen::Matrix3f VolumeTexture::world_to_grid() const {
    return world_to_grid_;
}

float VolumeTexture::max_abs_value() const {
    return max_abs_val_;
}

}  // namespace sbox::render
