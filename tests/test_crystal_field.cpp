#include "analysis/crystal_field.h"

#include <gtest/gtest.h>

namespace {

using sbox::analysis::DOrbitalEnergies;
using sbox::chem::CoordinationGeometry;

}  // namespace

TEST(CrystalFieldTest, OctahedralGroupingAndDelta) {
    DOrbitalEnergies d;
    d.dxy = -0.4;
    d.dxz = -0.4;
    d.dyz = -0.4;
    d.dz2 = 0.6;
    d.dx2y2 = 0.6;

    sbox::analysis::identify_splitting(d, CoordinationGeometry::Octahedral);
    ASSERT_EQ(d.groups.size(), 2u);
    EXPECT_EQ(d.groups[0].label, "t2g");
    EXPECT_EQ(d.groups[1].label, "e_g");
    EXPECT_NEAR(d.groups[0].average_energy, -0.4, 1.0e-9);
    EXPECT_NEAR(d.groups[1].average_energy, 0.6, 1.0e-9);
    EXPECT_NEAR(d.delta_oct(), 1.0, 1.0e-9);
}

TEST(CrystalFieldTest, TetrahedralGroupingAndDelta) {
    DOrbitalEnergies d;
    d.dxy = 0.2;
    d.dxz = 0.2;
    d.dyz = 0.2;
    d.dz2 = -0.3;
    d.dx2y2 = -0.3;

    sbox::analysis::identify_splitting(d, CoordinationGeometry::Tetrahedral);
    ASSERT_EQ(d.groups.size(), 2u);
    EXPECT_EQ(d.groups[0].label, "e");
    EXPECT_EQ(d.groups[1].label, "t2");
    EXPECT_NEAR(d.delta_tet(), 0.5, 1.0e-9);
}

TEST(CrystalFieldTest, DElectronCountForIronTwoPlus) {
    EXPECT_EQ(sbox::analysis::d_electron_count(26, 2), 6);
}

TEST(CrystalFieldTest, OctahedralCfseForHighSpinD6) {
    EXPECT_NEAR(sbox::analysis::octahedral_cfse_dq(6, true), -4.0, 1.0e-9);
}
