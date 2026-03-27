#include "core/elements.h"
#include "core/molecular_system.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <vector>

#include <gtest/gtest.h>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

sbox::chem::MolecularSystem make_water_geometry_for_com() {
    sbox::chem::MolecularSystem system;
    system.add_atom({8, Eigen::Vector3d(0.0, 0.0, 0.0), "O"});
    system.add_atom({1, Eigen::Vector3d(1.0, 0.0, 0.0), "H1"});
    system.add_atom({1, Eigen::Vector3d(0.0, 1.0, 0.0), "H2"});
    return system;
}

}  // namespace

TEST(MolecularSystemTest, AddAtomAndNumAtoms) {
    sbox::chem::MolecularSystem system;

    EXPECT_EQ(system.num_atoms(), 0);

    const int oxygen = system.add_atom({8, Eigen::Vector3d(0.0, 0.0, 0.0), "O"});
    const int hydrogen = system.add_atom({1, Eigen::Vector3d(1.8, 0.0, 0.0), "H"});

    EXPECT_EQ(oxygen, 0);
    EXPECT_EQ(hydrogen, 1);
    EXPECT_EQ(system.num_atoms(), 2);
    EXPECT_EQ(system.atom(0).Z, 8);
    EXPECT_EQ(system.atom(1).label, "H");
}

TEST(MolecularSystemTest, AddBondAndNeighbors) {
    sbox::chem::MolecularSystem system;
    system.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "C1"});
    system.add_atom({6, Eigen::Vector3d(1.5, 0.0, 0.0), "C2"});
    system.add_atom({1, Eigen::Vector3d(-1.0, 0.0, 0.0), "H"});

    const int bond0 = system.add_bond(0, 1, sbox::chem::BondOrder::Single);
    const int bond1 = system.add_bond(0, 2, sbox::chem::BondOrder::Single);

    EXPECT_EQ(bond0, 0);
    EXPECT_EQ(bond1, 1);
    EXPECT_EQ(system.num_bonds(), 2);

    std::vector<int> neighbors = system.neighbors(0);
    std::sort(neighbors.begin(), neighbors.end());
    EXPECT_EQ(neighbors, (std::vector<int>{1, 2}));
    EXPECT_EQ(system.coordination_number(0), 2);
}

TEST(MolecularSystemTest, RemoveAtomRemovesAssociatedBondsAndRenumbers) {
    sbox::chem::MolecularSystem system;
    system.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "A"});
    system.add_atom({6, Eigen::Vector3d(1.0, 0.0, 0.0), "B"});
    system.add_atom({6, Eigen::Vector3d(2.0, 0.0, 0.0), "C"});
    system.add_atom({6, Eigen::Vector3d(3.0, 0.0, 0.0), "D"});

    system.add_bond(0, 1, sbox::chem::BondOrder::Single);
    system.add_bond(1, 2, sbox::chem::BondOrder::Double);
    system.add_bond(2, 3, sbox::chem::BondOrder::Triple);

    system.remove_atom(1);

    EXPECT_EQ(system.num_atoms(), 3);
    EXPECT_EQ(system.num_bonds(), 1);
    EXPECT_EQ(system.bond(0).atom_i, 1);
    EXPECT_EQ(system.bond(0).atom_j, 2);
    EXPECT_EQ(system.bond(0).order, sbox::chem::BondOrder::Triple);
}

TEST(MolecularSystemTest, DistanceMatchesKnownPositions) {
    sbox::chem::MolecularSystem system;
    system.add_atom({1, Eigen::Vector3d(0.0, 0.0, 0.0), "H1"});
    system.add_atom({1, Eigen::Vector3d(3.0, 4.0, 0.0), "H2"});

    EXPECT_NEAR(system.distance(0, 1), 5.0, 1e-12);
}

TEST(MolecularSystemTest, AngleMatchesWaterLikeGeometry) {
    sbox::chem::MolecularSystem system;
    const double theta = 104.5 * kPi / 180.0;
    system.add_atom({1, Eigen::Vector3d(1.0, 0.0, 0.0), "H1"});
    system.add_atom({8, Eigen::Vector3d(0.0, 0.0, 0.0), "O"});
    system.add_atom({1, Eigen::Vector3d(std::cos(theta), std::sin(theta), 0.0), "H2"});

    EXPECT_NEAR(system.angle(0, 1, 2), theta, 1e-12);
}

TEST(MolecularSystemTest, DihedralMatchesKnownGeometry) {
    sbox::chem::MolecularSystem system;
    system.add_atom({6, Eigen::Vector3d(0.0, 1.0, 0.0), "A"});
    system.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "B"});
    system.add_atom({6, Eigen::Vector3d(1.0, 0.0, 0.0), "C"});
    system.add_atom({6, Eigen::Vector3d(1.0, 0.5, std::sqrt(3.0) * 0.5), "D"});

    EXPECT_NEAR(system.dihedral(0, 1, 2, 3), kPi / 3.0, 1e-12);
}

TEST(MolecularSystemTest, CenterOfMassUsesElementMasses) {
    const sbox::chem::MolecularSystem system = make_water_geometry_for_com();
    const double h_mass = sbox::elements::get_element(1).atomic_mass;
    const double o_mass = sbox::elements::get_element(8).atomic_mass;
    const double total_mass = o_mass + 2.0 * h_mass;

    const Eigen::Vector3d expected(h_mass / total_mass, h_mass / total_mass, 0.0);
    const Eigen::Vector3d actual = system.center_of_mass();

    EXPECT_NEAR(actual.x(), expected.x(), 1e-12);
    EXPECT_NEAR(actual.y(), expected.y(), 1e-12);
    EXPECT_NEAR(actual.z(), expected.z(), 1e-12);
}

TEST(MolecularSystemTest, CenterMovesCenterOfMassToOrigin) {
    sbox::chem::MolecularSystem system = make_water_geometry_for_com();

    system.center();
    const Eigen::Vector3d centered = system.center_of_mass();

    EXPECT_NEAR(centered.x(), 0.0, 1e-12);
    EXPECT_NEAR(centered.y(), 0.0, 1e-12);
    EXPECT_NEAR(centered.z(), 0.0, 1e-12);
}
