#pragma once

#include "io/cube_io.h"

#include <Eigen/Core>

namespace sbox::render {

class ESPSurface {
public:
    ESPSurface();
    ~ESPSurface();

    ESPSurface(const ESPSurface&) = delete;
    ESPSurface& operator=(const ESPSurface&) = delete;

    bool upload(const sbox::io::CubeData& density_cube,
                const sbox::io::CubeData& esp_cube);

    bool is_uploaded() const;

    void bind(unsigned int shader_id, int density_unit = 6, int esp_unit = 7) const;
    void unbind() const;

    Eigen::Vector3f origin() const;
    Eigen::Matrix3f world_to_grid() const;
    float esp_min() const;
    float esp_max() const;

private:
    unsigned int density_tex_ = 0;
    unsigned int esp_tex_ = 0;
    Eigen::Vector3f origin_ = Eigen::Vector3f::Zero();
    Eigen::Matrix3f world_to_grid_ = Eigen::Matrix3f::Identity();
    float esp_min_ = 0.0f;
    float esp_max_ = 0.0f;
    bool uploaded_ = false;
    mutable int bound_density_unit_ = 6;
    mutable int bound_esp_unit_ = 7;
};

}  // namespace sbox::render
