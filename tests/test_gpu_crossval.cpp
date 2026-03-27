#include "core/basis_set.h"
#include "core/gaussian_eval.h"
#include "renderer/basis_texture.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

#include <cmath>

namespace {

sbox::basis::MOData make_h2_sto3g() {
    sbox::basis::MOData mo_data;
    mo_data.basis.spherical = true;

    mo_data.atom_positions.push_back(Eigen::Vector3d(0.0, 0.0, -0.7));
    mo_data.atom_positions.push_back(Eigen::Vector3d(0.0, 0.0, 0.7));

    sbox::basis::BasisShell s1;
    s1.atom_index = 0;
    s1.angular_momentum = 0;
    s1.primitives = {{3.42525091, 0.15432897}, {0.62391373, 0.53532814}, {0.16885540, 0.44463454}};
    mo_data.basis.shells.push_back(s1);

    sbox::basis::BasisShell s2 = s1;
    s2.atom_index = 1;
    mo_data.basis.shells.push_back(s2);

    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    mo_data.coefficients = Eigen::MatrixXd(2, 2);
    mo_data.coefficients(0, 0) = inv_sqrt2;
    mo_data.coefficients(1, 0) = inv_sqrt2;
    mo_data.coefficients(0, 1) = inv_sqrt2;
    mo_data.coefficients(1, 1) = -inv_sqrt2;

    mo_data.energies = Eigen::VectorXd(2);
    mo_data.energies(0) = -0.6;
    mo_data.energies(1) = 0.7;

    mo_data.occupations = Eigen::VectorXd(2);
    mo_data.occupations(0) = 2.0;
    mo_data.occupations(1) = 0.0;

    mo_data.total_energy = -1.117;
    return mo_data;
}

}  // namespace

TEST(GpuCrossvalTest, BondingAndAntibondingMatchExpectedSymmetry) {
    const sbox::basis::MOData mo_data = make_h2_sto3g();

    const double bonding_mid = sbox::basis::evaluate_mo_at_point(mo_data, 0, Eigen::Vector3d(0.0, 0.0, 0.0));
    const double antibond_mid = sbox::basis::evaluate_mo_at_point(mo_data, 1, Eigen::Vector3d(0.0, 0.0, 0.0));

    EXPECT_GT(bonding_mid, 0.0);
    EXPECT_NEAR(antibond_mid, 0.0, 1e-12);
}

TEST(GpuCrossvalTest, BondingMoIsSymmetricAtBothAtoms) {
    const sbox::basis::MOData mo_data = make_h2_sto3g();

    const double at_a = sbox::basis::evaluate_mo_at_point(mo_data, 0, Eigen::Vector3d(0.0, 0.0, -0.7));
    const double at_b = sbox::basis::evaluate_mo_at_point(mo_data, 0, Eigen::Vector3d(0.0, 0.0, 0.7));

    EXPECT_GT(at_a, 0.0);
    EXPECT_NEAR(at_a, at_b, 1e-12);
}

TEST(GpuCrossvalTest, AntibondingMoIsAntisymmetricAtAtoms) {
    const sbox::basis::MOData mo_data = make_h2_sto3g();

    const double at_a = sbox::basis::evaluate_mo_at_point(mo_data, 1, Eigen::Vector3d(0.0, 0.0, -0.7));
    const double at_b = sbox::basis::evaluate_mo_at_point(mo_data, 1, Eigen::Vector3d(0.0, 0.0, 0.7));

    EXPECT_GT(at_a, 0.0);
    EXPECT_LT(at_b, 0.0);
    EXPECT_NEAR(std::abs(at_a), std::abs(at_b), 1e-12);
}

TEST(GpuCrossvalTest, FarFieldDecaysToZero) {
    const sbox::basis::MOData mo_data = make_h2_sto3g();
    const double far_val = sbox::basis::evaluate_mo_at_point(mo_data, 0, Eigen::Vector3d(0.0, 0.0, 50.0));
    EXPECT_NEAR(far_val, 0.0, 1e-8);
}

TEST(GpuCrossvalTest, BasisPackingMatchesExpectedCpuData) {
    const sbox::basis::MOData mo_data = make_h2_sto3g();
    sbox::render::BasisTextures textures;

    ASSERT_TRUE(textures.pack(mo_data));
    ASSERT_EQ(textures.num_shells(), 2);
    ASSERT_EQ(textures.num_primitives(), 6);
    ASSERT_EQ(textures.num_basis(), 2);
    ASSERT_EQ(textures.num_mo(), 2);

    const auto& shell_data = textures.shell_data();
    ASSERT_EQ(shell_data.size(), 8U);
    EXPECT_FLOAT_EQ(shell_data[0], 0.0f);
    EXPECT_FLOAT_EQ(shell_data[1], 0.0f);
    EXPECT_FLOAT_EQ(shell_data[2], -0.7f);
    EXPECT_FLOAT_EQ(shell_data[3], 0.0f);

    const auto& meta_data = textures.meta_data();
    ASSERT_EQ(meta_data.size(), 8U);
    EXPECT_FLOAT_EQ(meta_data[0], 3.0f);
    EXPECT_FLOAT_EQ(meta_data[1], 0.0f);
    EXPECT_FLOAT_EQ(meta_data[2], 1.0f);
    EXPECT_FLOAT_EQ(meta_data[3], 0.0f);

    const auto& prim_data = textures.primitive_data();
    ASSERT_EQ(prim_data.size(), 24U);
    EXPECT_FLOAT_EQ(prim_data[0], 3.42525091f);
    EXPECT_FLOAT_EQ(prim_data[1], 0.15432897f);
    EXPECT_FLOAT_EQ(prim_data[4], 0.62391373f);
    EXPECT_FLOAT_EQ(prim_data[5], 0.53532814f);

    const auto& mo_coeff = textures.mo_coeff_data();
    ASSERT_EQ(mo_coeff.size(), 8U);
    const float inv_sqrt2 = static_cast<float>(1.0 / std::sqrt(2.0));
    EXPECT_NEAR(mo_coeff[0], inv_sqrt2, 1e-6f);
    EXPECT_NEAR(mo_coeff[1], inv_sqrt2, 1e-6f);
    EXPECT_NEAR(mo_coeff[4], inv_sqrt2, 1e-6f);
    EXPECT_NEAR(mo_coeff[1 + 4], -inv_sqrt2, 1e-6f);
}
