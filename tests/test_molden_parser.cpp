#include "core/molden_parser.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <string>

namespace {

std::string h2_sto3g_path() {
    const std::filesystem::path here(__FILE__);
    return (here.parent_path() / "data" / "h2_sto3g.molden").string();
}

}  // namespace

TEST(MoldenParserTest, ParsesH2Sto3gFromFile) {
    const sbox::basis::MOData data = sbox::molden::parse_molden_file(h2_sto3g_path());

    ASSERT_EQ(data.atom_positions.size(), 2U);
    EXPECT_NEAR(data.atom_positions[0].z(), -0.7, 1e-10);
    EXPECT_NEAR(data.atom_positions[1].z(), 0.7, 1e-10);

    ASSERT_EQ(data.basis.shells.size(), 2U);
    EXPECT_TRUE(data.basis.spherical);
    EXPECT_EQ(data.basis.num_basis_functions(), 2);

    EXPECT_EQ(data.basis.shells[0].angular_momentum, 0);
    EXPECT_EQ(data.basis.shells[1].angular_momentum, 0);
    ASSERT_EQ(data.basis.shells[0].primitives.size(), 3U);
    ASSERT_EQ(data.basis.shells[1].primitives.size(), 3U);
    EXPECT_NEAR(data.basis.shells[0].primitives[0].exponent, 3.42525091, 1e-9);
    EXPECT_NEAR(data.basis.shells[0].primitives[1].coefficient, 0.53532814, 1e-9);

    ASSERT_EQ(data.coefficients.rows(), 2);
    ASSERT_EQ(data.coefficients.cols(), 2);
    EXPECT_NEAR(data.coefficients(0, 0), 0.54893404, 1e-8);
    EXPECT_NEAR(data.coefficients(1, 0), 0.54893404, 1e-8);
    EXPECT_NEAR(data.coefficients(0, 1), 1.21146408, 1e-8);
    EXPECT_NEAR(data.coefficients(1, 1), -1.21146408, 1e-8);

    ASSERT_EQ(data.energies.size(), 2);
    ASSERT_EQ(data.occupations.size(), 2);
    EXPECT_NEAR(data.energies(0), -0.5782029, 1e-8);
    EXPECT_NEAR(data.energies(1), 0.6702678, 1e-8);
    EXPECT_NEAR(data.occupations(0), 2.0, 1e-12);
    EXPECT_NEAR(data.occupations(1), 0.0, 1e-12);
}

TEST(MoldenParserTest, CanRenormalizeContractionCoefficients) {
    const sbox::basis::MOData default_data = sbox::molden::parse_molden_file(h2_sto3g_path());

    sbox::molden::ParseOptions options;
    options.contraction_coefficients_include_shell_normalization = false;
    const sbox::basis::MOData renorm_data = sbox::molden::parse_molden_file(h2_sto3g_path(), options);

    ASSERT_FALSE(default_data.basis.shells.empty());
    ASSERT_FALSE(renorm_data.basis.shells.empty());
    const double default_c0 = default_data.basis.shells[0].primitives[0].coefficient;
    const double renorm_c0 = renorm_data.basis.shells[0].primitives[0].coefficient;

    // STO-3G H-shell coefficients become larger after per-shell renormalization.
    EXPECT_GT(renorm_c0, default_c0);
}
