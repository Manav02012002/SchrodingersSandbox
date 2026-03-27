#include "core/molecular_system.h"

#include <Eigen/Core>

#include <gtest/gtest.h>

#include <cmath>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

}  // namespace

TEST(BondPerceptionTest, PerceivesWaterBondsOnly) {
    sbox::chem::MolecularSystem system;
    const double bond_length = 1.81;
    const double angle = 104.5 * kPi / 180.0;

    system.add_atom({8, Eigen::Vector3d(0.0, 0.0, 0.0), "O"});
    system.add_atom({1, Eigen::Vector3d(bond_length, 0.0, 0.0), "H1"});
    system.add_atom({1, Eigen::Vector3d(bond_length * std::cos(angle), bond_length * std::sin(angle), 0.0), "H2"});

    system.perceive_bonds();

    EXPECT_EQ(system.num_bonds(), 2);
    EXPECT_TRUE(system.has_bond(0, 1));
    EXPECT_TRUE(system.has_bond(0, 2));
    EXPECT_FALSE(system.has_bond(1, 2));
}

TEST(BondPerceptionTest, PerceivesMethaneCarbonHydrogenBonds) {
    sbox::chem::MolecularSystem system;
    const double bond_length = 2.06;
    const double scale = bond_length / std::sqrt(3.0);

    system.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "C"});
    system.add_atom({1, scale * Eigen::Vector3d(1.0, 1.0, 1.0), "H1"});
    system.add_atom({1, scale * Eigen::Vector3d(1.0, -1.0, -1.0), "H2"});
    system.add_atom({1, scale * Eigen::Vector3d(-1.0, 1.0, -1.0), "H3"});
    system.add_atom({1, scale * Eigen::Vector3d(-1.0, -1.0, 1.0), "H4"});

    system.perceive_bonds();

    EXPECT_EQ(system.num_bonds(), 4);
    for (int i = 1; i <= 4; ++i) {
        EXPECT_TRUE(system.has_bond(0, i));
    }
    EXPECT_FALSE(system.has_bond(1, 2));
    EXPECT_FALSE(system.has_bond(1, 3));
    EXPECT_FALSE(system.has_bond(1, 4));
}

TEST(BondPerceptionTest, PerceivesNitrogenBond) {
    sbox::chem::MolecularSystem system;
    system.add_atom({7, Eigen::Vector3d(0.0, 0.0, 0.0), "N1"});
    system.add_atom({7, Eigen::Vector3d(2.07, 0.0, 0.0), "N2"});

    system.perceive_bonds();

    EXPECT_EQ(system.num_bonds(), 1);
    EXPECT_TRUE(system.has_bond(0, 1));
}

TEST(BondPerceptionTest, PerceivesTwoSeparatedHydrogenMolecules) {
    sbox::chem::MolecularSystem system;
    system.add_atom({1, Eigen::Vector3d(0.0, 0.0, 0.0), "H1"});
    system.add_atom({1, Eigen::Vector3d(1.1, 0.0, 0.0), "H2"});
    system.add_atom({1, Eigen::Vector3d(20.0, 0.0, 0.0), "H3"});
    system.add_atom({1, Eigen::Vector3d(21.1, 0.0, 0.0), "H4"});

    system.perceive_bonds();

    EXPECT_EQ(system.num_bonds(), 2);
    EXPECT_TRUE(system.has_bond(0, 1));
    EXPECT_TRUE(system.has_bond(2, 3));
    EXPECT_FALSE(system.has_bond(0, 2));
    EXPECT_FALSE(system.has_bond(0, 3));
    EXPECT_FALSE(system.has_bond(1, 2));
    EXPECT_FALSE(system.has_bond(1, 3));
}

TEST(BondPerceptionTest, HasBondReflectsPerceivedConnectivity) {
    sbox::chem::MolecularSystem system;
    system.add_atom({8, Eigen::Vector3d(0.0, 0.0, 0.0), "O"});
    system.add_atom({1, Eigen::Vector3d(1.81, 0.0, 0.0), "H1"});
    system.add_atom({1, Eigen::Vector3d(-1.81, 0.0, 0.0), "H2"});

    system.perceive_bonds();

    EXPECT_TRUE(system.has_bond(0, 1));
    EXPECT_TRUE(system.has_bond(1, 0));
    EXPECT_TRUE(system.has_bond(0, 2));
    EXPECT_FALSE(system.has_bond(1, 2));
}
