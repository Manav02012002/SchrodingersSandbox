#include "renderer/basis_texture.h"

#include <glad/gl.h>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace sbox::render {

namespace {

int num_basis_in_shell(int angular_momentum, bool spherical) {
    switch (angular_momentum) {
    case 0:
        return 1;
    case 1:
        return 3;
    case 2:
        return spherical ? 5 : 6;
    case 3:
        return spherical ? 7 : 10;
    default:
        throw std::runtime_error("Unsupported shell angular momentum: " + std::to_string(angular_momentum));
    }
}

void configure_texture(unsigned int texture) {
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}

void upload_rgba32f(unsigned int texture, int width, int height, const std::vector<float>& data) {
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_RGBA32F,
                 width,
                 height,
                 0,
                 GL_RGBA,
                 GL_FLOAT,
                 data.empty() ? nullptr : data.data());
}

void set_sampler_uniform(unsigned int shader_id, const char* name, int value) {
    const int location = glGetUniformLocation(shader_id, name);
    if (location >= 0) {
        glUniform1i(location, value);
    }
}

void ensure_texture(unsigned int* texture) {
    if (*texture != 0) {
        return;
    }
    glGenTextures(1, texture);
    configure_texture(*texture);
}

}  // namespace

BasisTextures::BasisTextures() = default;

BasisTextures::~BasisTextures() {
    if (tex_shells_ != 0) {
        glDeleteTextures(1, &tex_shells_);
    }
    if (tex_meta_ != 0) {
        glDeleteTextures(1, &tex_meta_);
    }
    if (tex_primitives_ != 0) {
        glDeleteTextures(1, &tex_primitives_);
    }
    if (tex_mo_coeffs_ != 0) {
        glDeleteTextures(1, &tex_mo_coeffs_);
    }
}

bool BasisTextures::pack(const sbox::basis::MOData& mo_data) {
    num_shells_ = static_cast<int>(mo_data.basis.shells.size());
    num_basis_ = mo_data.basis.num_basis_functions();
    num_mo_ = static_cast<int>(mo_data.coefficients.cols());
    spherical_ = mo_data.basis.spherical;
    uploaded_ = false;

    num_primitives_ = 0;
    for (const sbox::basis::BasisShell& shell : mo_data.basis.shells) {
        num_primitives_ += static_cast<int>(shell.primitives.size());
    }

    if (num_shells_ > MAX_SHELLS || num_primitives_ > MAX_PRIMITIVES || num_basis_ > MAX_BASIS || num_mo_ > MAX_MO) {
        shell_data_.clear();
        meta_data_.clear();
        prim_data_.clear();
        mo_data_.clear();
        return false;
    }

    if (mo_data.coefficients.rows() != num_basis_) {
        throw std::runtime_error("MO coefficient row count does not match basis function count");
    }

    shell_data_.assign(static_cast<std::size_t>(num_shells_) * 4, 0.0f);
    meta_data_.assign(static_cast<std::size_t>(num_shells_) * 4, 0.0f);
    prim_data_.assign(static_cast<std::size_t>(num_primitives_) * 4, 0.0f);

    const int mo_tex_height = std::max(1, (num_mo_ + 3) / 4);
    mo_data_.assign(static_cast<std::size_t>(num_basis_) * static_cast<std::size_t>(mo_tex_height) * 4, 0.0f);

    int primitive_offset = 0;
    int basis_offset = 0;
    for (int shell_index = 0; shell_index < num_shells_; ++shell_index) {
        const sbox::basis::BasisShell& shell = mo_data.basis.shells[static_cast<std::size_t>(shell_index)];
        if (shell.atom_index < 0 || shell.atom_index >= static_cast<int>(mo_data.atom_positions.size())) {
            throw std::runtime_error("Basis shell atom index out of range");
        }

        const Eigen::Vector3d& center = mo_data.atom_positions[static_cast<std::size_t>(shell.atom_index)];
        const int shell_data_offset = shell_index * 4;
        shell_data_[static_cast<std::size_t>(shell_data_offset + 0)] = static_cast<float>(center.x());
        shell_data_[static_cast<std::size_t>(shell_data_offset + 1)] = static_cast<float>(center.y());
        shell_data_[static_cast<std::size_t>(shell_data_offset + 2)] = static_cast<float>(center.z());
        shell_data_[static_cast<std::size_t>(shell_data_offset + 3)] = static_cast<float>(shell.angular_momentum);

        const int basis_in_shell = num_basis_in_shell(shell.angular_momentum, spherical_);
        meta_data_[static_cast<std::size_t>(shell_data_offset + 0)] = static_cast<float>(shell.primitives.size());
        meta_data_[static_cast<std::size_t>(shell_data_offset + 1)] = static_cast<float>(primitive_offset);
        meta_data_[static_cast<std::size_t>(shell_data_offset + 2)] = static_cast<float>(basis_in_shell);
        meta_data_[static_cast<std::size_t>(shell_data_offset + 3)] = static_cast<float>(basis_offset);

        for (const sbox::basis::GaussianPrimitive& primitive : shell.primitives) {
            const int primitive_data_offset = primitive_offset * 4;
            prim_data_[static_cast<std::size_t>(primitive_data_offset + 0)] = static_cast<float>(primitive.exponent);
            prim_data_[static_cast<std::size_t>(primitive_data_offset + 1)] = static_cast<float>(primitive.coefficient);
            prim_data_[static_cast<std::size_t>(primitive_data_offset + 2)] = 0.0f;
            prim_data_[static_cast<std::size_t>(primitive_data_offset + 3)] = 0.0f;
            ++primitive_offset;
        }

        basis_offset += basis_in_shell;
    }

    for (int mu = 0; mu < num_basis_; ++mu) {
        for (int mo = 0; mo < num_mo_; ++mo) {
            const int texel_y = mo / 4;
            const int channel = mo % 4;
            const std::size_t offset =
                (static_cast<std::size_t>(mu) * static_cast<std::size_t>(mo_tex_height)
                 + static_cast<std::size_t>(texel_y))
                * 4
                + static_cast<std::size_t>(channel);
            mo_data_[offset] = static_cast<float>(mo_data.coefficients(mu, mo));
        }
    }

    uploaded_ = true;
    return true;
}

bool BasisTextures::upload(const sbox::basis::MOData& mo_data) {
    if (!pack(mo_data)) {
        return false;
    }

    ensure_texture(&tex_shells_);
    ensure_texture(&tex_meta_);
    ensure_texture(&tex_primitives_);
    ensure_texture(&tex_mo_coeffs_);

    const int mo_tex_height = std::max(1, (num_mo_ + 3) / 4);
    upload_rgba32f(tex_shells_, std::max(1, num_shells_), 1, shell_data_);
    upload_rgba32f(tex_meta_, std::max(1, num_shells_), 1, meta_data_);
    upload_rgba32f(tex_primitives_, std::max(1, num_primitives_), 1, prim_data_);
    upload_rgba32f(tex_mo_coeffs_, std::max(1, num_basis_), mo_tex_height, mo_data_);
    glBindTexture(GL_TEXTURE_2D, 0);
    return true;
}

void BasisTextures::bind(unsigned int shader_id, int first_unit) const {
    bound_first_unit_ = first_unit;

    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(first_unit + 0));
    glBindTexture(GL_TEXTURE_2D, tex_shells_);
    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(first_unit + 1));
    glBindTexture(GL_TEXTURE_2D, tex_meta_);
    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(first_unit + 2));
    glBindTexture(GL_TEXTURE_2D, tex_primitives_);
    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(first_unit + 3));
    glBindTexture(GL_TEXTURE_2D, tex_mo_coeffs_);

    set_sampler_uniform(shader_id, "u_shell_desc", first_unit + 0);
    set_sampler_uniform(shader_id, "u_shell_meta", first_unit + 1);
    set_sampler_uniform(shader_id, "u_primitives", first_unit + 2);
    set_sampler_uniform(shader_id, "u_mo_coeffs", first_unit + 3);
}

void BasisTextures::unbind() const {
    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(bound_first_unit_ + 0));
    glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(bound_first_unit_ + 1));
    glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(bound_first_unit_ + 2));
    glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE0 + static_cast<unsigned int>(bound_first_unit_ + 3));
    glBindTexture(GL_TEXTURE_2D, 0);
}

int BasisTextures::num_shells() const {
    return num_shells_;
}

int BasisTextures::num_basis() const {
    return num_basis_;
}

int BasisTextures::num_mo() const {
    return num_mo_;
}

int BasisTextures::num_primitives() const {
    return num_primitives_;
}

bool BasisTextures::is_spherical() const {
    return spherical_;
}

bool BasisTextures::is_uploaded() const {
    return uploaded_;
}

const std::vector<float>& BasisTextures::shell_data() const {
    return shell_data_;
}

const std::vector<float>& BasisTextures::meta_data() const {
    return meta_data_;
}

const std::vector<float>& BasisTextures::primitive_data() const {
    return prim_data_;
}

const std::vector<float>& BasisTextures::mo_coeff_data() const {
    return mo_data_;
}

}  // namespace sbox::render
