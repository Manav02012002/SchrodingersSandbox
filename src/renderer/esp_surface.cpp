#include "renderer/esp_surface.h"

#include <glad/gl.h>

#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace sbox::render {

namespace {

void set_sampler_uniform(unsigned int shader_id, const char* name, int value) {
    const int location = glGetUniformLocation(shader_id, name);
    if (location >= 0) {
        glUniform1i(location, value);
    }
}

bool nearly_equal(float a, float b, float tol = 1.0e-4f) {
    return std::abs(a - b) <= tol;
}

bool nearly_equal_vec(const Eigen::Vector3f& a, const Eigen::Vector3f& b, float tol = 1.0e-4f) {
    return (a - b).cwiseAbs().maxCoeff() <= tol;
}

void configure_3d_texture(unsigned int tex, int nx, int ny, int nz, const float* data) {
    glBindTexture(GL_TEXTURE_3D, tex);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);
    const float border[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    glTexParameterfv(GL_TEXTURE_3D, GL_TEXTURE_BORDER_COLOR, border);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, nx, ny, nz, 0, GL_RED, GL_FLOAT, data);
    glBindTexture(GL_TEXTURE_3D, 0);
}

}  // namespace

ESPSurface::ESPSurface() {
    glGenTextures(1, &density_tex_);
    glGenTextures(1, &esp_tex_);
}

ESPSurface::~ESPSurface() {
    if (density_tex_ != 0) {
        glDeleteTextures(1, &density_tex_);
    }
    if (esp_tex_ != 0) {
        glDeleteTextures(1, &esp_tex_);
    }
}

bool ESPSurface::upload(const sbox::io::CubeData& density_cube,
                        const sbox::io::CubeData& esp_cube) {
    if (density_cube.nx != esp_cube.nx || density_cube.ny != esp_cube.ny || density_cube.nz != esp_cube.nz) {
        throw std::runtime_error("ESP cube grid dimensions do not match density cube");
    }

    const Eigen::Vector3f density_origin = density_cube.origin.cast<float>();
    const Eigen::Vector3f esp_origin = esp_cube.origin.cast<float>();
    if (!nearly_equal_vec(density_origin, esp_origin)) {
        throw std::runtime_error("ESP cube origin does not match density cube");
    }

    const Eigen::Matrix3f density_grid =
        (Eigen::Matrix3f() << density_cube.step_x.cast<float>(), density_cube.step_y.cast<float>(), density_cube.step_z.cast<float>()).finished();
    const Eigen::Matrix3f esp_grid =
        (Eigen::Matrix3f() << esp_cube.step_x.cast<float>(), esp_cube.step_y.cast<float>(), esp_cube.step_z.cast<float>()).finished();
    if (!nearly_equal_vec(density_grid.col(0), esp_grid.col(0)) ||
        !nearly_equal_vec(density_grid.col(1), esp_grid.col(1)) ||
        !nearly_equal_vec(density_grid.col(2), esp_grid.col(2))) {
        throw std::runtime_error("ESP cube grid transform does not match density cube");
    }

    origin_ = density_origin;
    const float det = density_grid.determinant();
    if (std::abs(det) < 1.0e-8f) {
        throw std::runtime_error("ESP surface grid transform is singular");
    }
    world_to_grid_ = density_grid.inverse();

    configure_3d_texture(density_tex_, density_cube.nx, density_cube.ny, density_cube.nz, density_cube.data.data());
    configure_3d_texture(esp_tex_, esp_cube.nx, esp_cube.ny, esp_cube.nz, esp_cube.data.data());

    if (esp_cube.data.empty()) {
        esp_min_ = -1.0f;
        esp_max_ = 1.0f;
    } else {
        const double mean = std::accumulate(esp_cube.data.begin(), esp_cube.data.end(), 0.0) /
                            static_cast<double>(esp_cube.data.size());
        double variance = 0.0;
        for (float value : esp_cube.data) {
            const double dv = static_cast<double>(value) - mean;
            variance += dv * dv;
        }
        variance /= static_cast<double>(esp_cube.data.size());
        const double sigma = std::sqrt(std::max(variance, 0.0));
        const double lo = mean - 3.0 * sigma;
        const double hi = mean + 3.0 * sigma;

        bool any = false;
        float min_v = 0.0f;
        float max_v = 0.0f;
        for (float value : esp_cube.data) {
            if (static_cast<double>(value) < lo || static_cast<double>(value) > hi) {
                continue;
            }
            if (!any) {
                min_v = value;
                max_v = value;
                any = true;
            } else {
                min_v = std::min(min_v, value);
                max_v = std::max(max_v, value);
            }
        }
        if (!any) {
            auto [min_it, max_it] = std::minmax_element(esp_cube.data.begin(), esp_cube.data.end());
            min_v = *min_it;
            max_v = *max_it;
        }
        if (nearly_equal(min_v, max_v, 1.0e-6f)) {
            min_v -= 1.0f;
            max_v += 1.0f;
        }
        esp_min_ = min_v;
        esp_max_ = max_v;
    }

    uploaded_ = true;
    return true;
}

bool ESPSurface::is_uploaded() const {
    return uploaded_;
}

void ESPSurface::bind(unsigned int shader_id, int density_unit, int esp_unit) const {
    bound_density_unit_ = density_unit;
    bound_esp_unit_ = esp_unit;
    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(density_unit));
    glBindTexture(GL_TEXTURE_3D, density_tex_);
    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(esp_unit));
    glBindTexture(GL_TEXTURE_3D, esp_tex_);
    set_sampler_uniform(shader_id, "u_density", density_unit);
    set_sampler_uniform(shader_id, "u_esp", esp_unit);
}

void ESPSurface::unbind() const {
    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(bound_density_unit_));
    glBindTexture(GL_TEXTURE_3D, 0);
    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(bound_esp_unit_));
    glBindTexture(GL_TEXTURE_3D, 0);
}

Eigen::Vector3f ESPSurface::origin() const {
    return origin_;
}

Eigen::Matrix3f ESPSurface::world_to_grid() const {
    return world_to_grid_;
}

float ESPSurface::esp_min() const {
    return esp_min_;
}

float ESPSurface::esp_max() const {
    return esp_max_;
}

}  // namespace sbox::render
