#include "core/valence.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

#include <cmath>

namespace {

constexpr double kDeg = 180.0 / 3.14159265358979323846;

}  // namespace

TEST(ValenceTest, DefaultValenceValues) {
    EXPECT_EQ(sbox::chem::default_valence(6), 4);
    EXPECT_EQ(sbox::chem::default_valence(7), 3);
    EXPECT_EQ(sbox::chem::default_valence(8), 2);
    EXPECT_EQ(sbox::chem::default_valence(1), 1);
    EXPECT_EQ(sbox::chem::default_valence(26), -1);
}

TEST(ValenceTest, MissingHydrogensCounts) {
    sbox::chem::MolecularSystem methane_center;
    methane_center.add_atom({6, Eigen::Vector3d::Zero(), "", 0});
    EXPECT_EQ(sbox::chem::missing_hydrogens(methane_center, 0), 4);

    sbox::chem::MolecularSystem carbon_two_bonds;
    carbon_two_bonds.add_atom({6, Eigen::Vector3d::Zero(), "", 0});
    carbon_two_bonds.add_atom({6, Eigen::Vector3d(1.0, 0.0, 0.0), "", 0});
    carbon_two_bonds.add_atom({6, Eigen::Vector3d(-1.0, 0.0, 0.0), "", 0});
    carbon_two_bonds.add_bond(0, 1, sbox::chem::BondOrder::Single);
    carbon_two_bonds.add_bond(0, 2, sbox::chem::BondOrder::Single);
    EXPECT_EQ(sbox::chem::missing_hydrogens(carbon_two_bonds, 0), 2);

    sbox::chem::MolecularSystem oxygen_one_bond;
    oxygen_one_bond.add_atom({8, Eigen::Vector3d::Zero(), "", 0});
    oxygen_one_bond.add_atom({6, Eigen::Vector3d(1.0, 0.0, 0.0), "", 0});
    oxygen_one_bond.add_bond(0, 1, sbox::chem::BondOrder::Single);
    EXPECT_EQ(sbox::chem::missing_hydrogens(oxygen_one_bond, 0), 1);
}

TEST(ValenceTest, AddHydrogensToBareCarbonGivesMethaneGeometry) {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({6, Eigen::Vector3d::Zero(), "", 0});

    sbox::chem::add_hydrogens(mol);
    ASSERT_EQ(mol.num_atoms(), 5);
    ASSERT_EQ(mol.num_bonds(), 4);

    for (int i = 1; i < mol.num_atoms(); ++i) {
        EXPECT_NEAR(mol.distance(0, i), 2.06, 1e-2);
    }

    for (int i = 1; i < mol.num_atoms(); ++i) {
        for (int j = i + 1; j < mol.num_atoms(); ++j) {
            EXPECT_NEAR(mol.angle(i, 0, j) * kDeg, 109.47, 2.0);
        }
    }
}

TEST(ValenceTest, AddHydrogensToWaterLikeOxygenAddsTwoHydrogens) {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({8, Eigen::Vector3d::Zero(), "", 0});

    sbox::chem::add_hydrogens(mol);
    EXPECT_EQ(mol.num_atoms(), 3);
    EXPECT_EQ(mol.num_bonds(), 2);
}

TEST(ValenceTest, AddHydrogensToEthyleneRespectsSp2Count) {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({6, Eigen::Vector3d(-1.27, 0.0, 0.0), "", 0});
    mol.add_atom({6, Eigen::Vector3d(1.27, 0.0, 0.0), "", 0});
    mol.add_bond(0, 1, sbox::chem::BondOrder::Double);

    sbox::chem::add_hydrogens(mol);
    EXPECT_EQ(mol.num_atoms(), 6);
    EXPECT_EQ(mol.num_bonds(), 5);

    int carbon0_h = 0;
    int carbon1_h = 0;
    for (const auto& bond : mol.bonds()) {
        if ((bond.atom_i == 0 && mol.atom(bond.atom_j).Z == 1) || (bond.atom_j == 0 && mol.atom(bond.atom_i).Z == 1)) {
            ++carbon0_h;
        }
        if ((bond.atom_i == 1 && mol.atom(bond.atom_j).Z == 1) || (bond.atom_j == 1 && mol.atom(bond.atom_i).Z == 1)) {
            ++carbon1_h;
        }
    }
    EXPECT_EQ(carbon0_h, 2);
    EXPECT_EQ(carbon1_h, 2);

    std::vector<int> neighbors0;
    for (int n : mol.neighbors(0)) {
        if (n != 1) {
            neighbors0.push_back(n);
        }
    }
    ASSERT_EQ(neighbors0.size(), 2u);
    EXPECT_NEAR(mol.angle(neighbors0[0], 0, neighbors0[1]) * kDeg, 120.0, 2.0);
}

TEST(ValenceTest, RemoveHydrogensLeavesHeavyAtoms) {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({6, Eigen::Vector3d::Zero(), "", 0});
    sbox::chem::add_hydrogens(mol);
    ASSERT_EQ(mol.num_atoms(), 5);

    sbox::chem::remove_hydrogens(mol);
    EXPECT_EQ(mol.num_atoms(), 1);
    EXPECT_EQ(mol.atom(0).Z, 6);
    EXPECT_EQ(mol.num_bonds(), 0);
}
