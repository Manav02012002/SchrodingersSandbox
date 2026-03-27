#pragma once

#include "io/cube_io.h"

#include <Eigen/Core>

namespace sbox::render {

class VolumeTexture {
public:
    VolumeTexture();
    ~VolumeTexture();

    VolumeTexture(const VolumeTexture&) = delete;
    VolumeTexture& operator=(const VolumeTexture&) = delete;

    bool upload(const sbox::io::CubeData& cube);
    bool upload(const Eigen::Vector3d& origin,
                const Eigen::Vector3d& step_x,
                const Eigen::Vector3d& step_y,
                const Eigen::Vector3d& step_z,
                int nx,
                int ny,
                int nz,
                const float* data);

    void bind(unsigned int shader_id, int texture_unit = 5) const;
    void unbind() const;

    bool is_uploaded() const;
    int nx() const;
    int ny() const;
    int nz() const;

    Eigen::Vector3f origin() const;
    Eigen::Matrix3f grid_to_world() const;
    Eigen::Matrix3f world_to_grid() const;
    float max_abs_value() const;

private:
    unsigned int texture_3d_ = 0;
    int nx_ = 0;
    int ny_ = 0;
    int nz_ = 0;
    Eigen::Vector3f origin_ = Eigen::Vector3f::Zero();
    Eigen::Matrix3f grid_to_world_ = Eigen::Matrix3f::Identity();
    Eigen::Matrix3f world_to_grid_ = Eigen::Matrix3f::Identity();
    float max_abs_val_ = 1.0f;
    bool uploaded_ = false;

    mutable int bound_texture_unit_ = 5;
};

}  // namespace sbox::render
