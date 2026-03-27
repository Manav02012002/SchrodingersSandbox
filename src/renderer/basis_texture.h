#pragma once

#include "core/basis_set.h"

#include <vector>

namespace sbox::render {

inline constexpr int MAX_SHELLS = 512;
inline constexpr int MAX_PRIMITIVES = 2048;
inline constexpr int MAX_BASIS = 2048;
inline constexpr int MAX_MO = 512;

class BasisTextures {
public:
    BasisTextures();
    ~BasisTextures();

    BasisTextures(const BasisTextures&) = delete;
    BasisTextures& operator=(const BasisTextures&) = delete;

    bool pack(const sbox::basis::MOData& mo_data);
    bool upload(const sbox::basis::MOData& mo_data);
    void bind(unsigned int shader_id, int first_unit = 1) const;
    void unbind() const;

    int num_shells() const;
    int num_basis() const;
    int num_mo() const;
    int num_primitives() const;
    bool is_spherical() const;
    bool is_uploaded() const;

    const std::vector<float>& shell_data() const;
    const std::vector<float>& meta_data() const;
    const std::vector<float>& primitive_data() const;
    const std::vector<float>& mo_coeff_data() const;

private:
    unsigned int tex_shells_ = 0;
    unsigned int tex_meta_ = 0;
    unsigned int tex_primitives_ = 0;
    unsigned int tex_mo_coeffs_ = 0;

    int num_shells_ = 0;
    int num_basis_ = 0;
    int num_mo_ = 0;
    int num_primitives_ = 0;
    bool spherical_ = true;
    bool uploaded_ = false;

    std::vector<float> shell_data_;
    std::vector<float> meta_data_;
    std::vector<float> prim_data_;
    std::vector<float> mo_data_;

    mutable int bound_first_unit_ = 1;
};

}  // namespace sbox::render
