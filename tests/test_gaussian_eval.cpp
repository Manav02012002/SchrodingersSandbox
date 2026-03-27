#include "core/gaussian_eval.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

#include <cmath>
#include <vector>

TEST(GaussianEvalTest, PrimitiveSTypeMatchesAnalyticValues) {
    EXPECT_DOUBLE_EQ(sbox::basis::evaluate_primitive(1.0, 1.0, 0, 0, 0, 0.0, 0.0, 0.0), 1.0);
    EXPECT_NEAR(sbox::basis::evaluate_primitive(1.0, 1.0, 0, 0, 0, 1.0, 0.0, 0.0), std::exp(-1.0), 1e-12);
}

TEST(GaussianEvalTest, PShellEvaluatesCartesianComponents) {
    sbox::basis::BasisShell shell;
    shell.atom_index = 0;
    shell.angular_momentum = 1;
    shell.primitives.push_back({1.0, 1.0});

    std::vector<double> values(3, 0.0);
    const double count = sbox::basis::evaluate_shell(
        shell, Eigen::Vector3d::Zero(), Eigen::Vector3d(1.0, 0.0, 0.0), 0, values);

    EXPECT_DOUBLE_EQ(count, 3.0);
    EXPECT_NEAR(values[0], std::exp(-1.0), 1e-12);
    EXPECT_DOUBLE_EQ(values[1], 0.0);
    EXPECT_DOUBLE_EQ(values[2], 0.0);
}

TEST(GaussianEvalTest, Sto3GHydrogenAtNucleusMatchesStoredContractions) {
    sbox::basis::BasisShell shell;
    shell.atom_index = 0;
    shell.angular_momentum = 0;
    shell.primitives.push_back({3.42525091, 0.15432897});
    shell.primitives.push_back({0.62391373, 0.53532814});
    shell.primitives.push_back({0.16885540, 0.44463454});

    std::vector<double> values(1, 0.0);
    (void)sbox::basis::evaluate_shell(shell, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 0, values);

    const double expected = 0.15432897 + 0.53532814 + 0.44463454;
    EXPECT_NEAR(values[0], expected, 1e-12);
}

TEST(GaussianEvalTest, TrivialSingleBasisMoMatchesBasisValue) {
    sbox::basis::MOData mo_data;
    mo_data.atom_positions.push_back(Eigen::Vector3d::Zero());
    mo_data.basis.spherical = true;

    sbox::basis::BasisShell shell;
    shell.atom_index = 0;
    shell.angular_momentum = 0;
    shell.primitives.push_back({1.0, 1.0});
    mo_data.basis.shells.push_back(shell);

    mo_data.coefficients = Eigen::MatrixXd::Ones(1, 1);

    const Eigen::Vector3d point(0.5, 0.0, 0.0);
    Eigen::VectorXd basis_values;
    sbox::basis::evaluate_basis_at_point(mo_data, point, basis_values);
    const double mo_value = sbox::basis::evaluate_mo_at_point(mo_data, 0, point);

    ASSERT_EQ(basis_values.size(), 1);
    EXPECT_NEAR(mo_value, basis_values(0), 1e-12);
}

TEST(GaussianEvalTest, SphericalDShellUsesMoldenOrdering) {
    sbox::basis::BasisShell shell;
    shell.atom_index = 0;
    shell.angular_momentum = 2;
    shell.primitives.push_back({0.0, 1.0});

    std::vector<double> values(5, 0.0);
    const double count = sbox::basis::evaluate_shell(
        shell, Eigen::Vector3d::Zero(), Eigen::Vector3d(0.0, 0.0, 1.0), 0, values);

    EXPECT_DOUBLE_EQ(count, 5.0);
    EXPECT_NEAR(values[0], 1.0, 1e-12);
    EXPECT_NEAR(values[1], 0.0, 1e-12);
    EXPECT_NEAR(values[2], 0.0, 1e-12);
    EXPECT_NEAR(values[3], 0.0, 1e-12);
    EXPECT_NEAR(values[4], 0.0, 1e-12);
}

TEST(GaussianEvalTest, SphericalFShellUsesMoldenOrdering) {
    sbox::basis::BasisShell shell;
    shell.atom_index = 0;
    shell.angular_momentum = 3;
    shell.primitives.push_back({0.0, 1.0});

    std::vector<double> values(7, 0.0);
    const double count = sbox::basis::evaluate_shell(
        shell, Eigen::Vector3d::Zero(), Eigen::Vector3d(0.0, 0.0, 1.0), 0, values);

    EXPECT_DOUBLE_EQ(count, 7.0);
    EXPECT_NEAR(values[0], 1.0, 1e-12);
    EXPECT_NEAR(values[1], 0.0, 1e-12);
    EXPECT_NEAR(values[2], 0.0, 1e-12);
    EXPECT_NEAR(values[3], 0.0, 1e-12);
    EXPECT_NEAR(values[4], 0.0, 1e-12);
    EXPECT_NEAR(values[5], 0.0, 1e-12);
    EXPECT_NEAR(values[6], 0.0, 1e-12);
}

TEST(GaussianEvalTest, BasisEvaluationSupportsSphericalDAndF) {
    sbox::basis::MOData mo_data;
    mo_data.atom_positions.push_back(Eigen::Vector3d::Zero());
    mo_data.basis.spherical = true;

    sbox::basis::BasisShell d_shell;
    d_shell.atom_index = 0;
    d_shell.angular_momentum = 2;
    d_shell.primitives.push_back({0.0, 1.0});
    mo_data.basis.shells.push_back(d_shell);

    sbox::basis::BasisShell f_shell;
    f_shell.atom_index = 0;
    f_shell.angular_momentum = 3;
    f_shell.primitives.push_back({0.0, 1.0});
    mo_data.basis.shells.push_back(f_shell);

    mo_data.coefficients = Eigen::MatrixXd::Zero(12, 1);
    Eigen::VectorXd basis_values;
    sbox::basis::evaluate_basis_at_point(mo_data, Eigen::Vector3d(0.0, 0.0, 1.0), basis_values);

    ASSERT_EQ(basis_values.size(), 12);
    EXPECT_NEAR(basis_values(0), 1.0, 1e-12);
    EXPECT_NEAR(basis_values(5), 1.0, 1e-12);
}
