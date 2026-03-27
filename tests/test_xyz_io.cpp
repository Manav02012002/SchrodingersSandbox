#include "core/molecular_system.h"
#include "io/xyz_io.h"

#include <Eigen/Core>

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <string>

namespace {

sbox::chem::MolecularSystem make_methane() {
    sbox::chem::MolecularSystem mol;
    mol.set_name("methane");

    const double bond_length = 2.06;
    const double scale = bond_length / std::sqrt(3.0);

    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "C"});
    mol.add_atom({1, scale * Eigen::Vector3d(1.0, 1.0, 1.0), "H1"});
    mol.add_atom({1, scale * Eigen::Vector3d(1.0, -1.0, -1.0), "H2"});
    mol.add_atom({1, scale * Eigen::Vector3d(-1.0, 1.0, -1.0), "H3"});
    mol.add_atom({1, scale * Eigen::Vector3d(-1.0, -1.0, 1.0), "H4"});

    return mol;
}

std::filesystem::path temp_xyz_path() {
    return std::filesystem::temp_directory_path() / "schrodingerssandbox_test_methane.xyz";
}

}  // namespace

TEST(XyzIoTest, ReadXyzStringParsesWaterAndPerceivesBonds) {
    const std::string xyz =
        "3\n"
        "water\n"
        "O  0.000000  0.000000  0.117300\n"
        "H  0.000000  0.756950 -0.469200\n"
        "H  0.000000 -0.756950 -0.469200\n";

    const sbox::chem::MolecularSystem mol = sbox::io::read_xyz_string(xyz);

    ASSERT_EQ(mol.num_atoms(), 3);
    EXPECT_EQ(mol.atom(0).Z, 8);
    EXPECT_EQ(mol.atom(1).Z, 1);
    EXPECT_EQ(mol.atom(2).Z, 1);
    EXPECT_EQ(mol.name(), "water");
    EXPECT_EQ(mol.num_bonds(), 2);
}

TEST(XyzIoTest, RoundTripPreservesAtomsAndPositions) {
    const sbox::chem::MolecularSystem original = make_methane();
    const std::filesystem::path path = temp_xyz_path();

    sbox::io::write_xyz(path.string(), original);
    const sbox::chem::MolecularSystem round_tripped = sbox::io::read_xyz(path.string());
    std::filesystem::remove(path);

    ASSERT_EQ(round_tripped.num_atoms(), original.num_atoms());
    for (int i = 0; i < original.num_atoms(); ++i) {
        EXPECT_EQ(round_tripped.atom(i).Z, original.atom(i).Z);
        EXPECT_NEAR(round_tripped.atom(i).position.x(), original.atom(i).position.x(), 1e-6);
        EXPECT_NEAR(round_tripped.atom(i).position.y(), original.atom(i).position.y(), 1e-6);
        EXPECT_NEAR(round_tripped.atom(i).position.z(), original.atom(i).position.z(), 1e-6);
    }
}

TEST(XyzIoTest, ReadSingleAtomXyzString) {
    const sbox::chem::MolecularSystem mol = sbox::io::read_xyz_string("1\nhelium\nHe 0.0 0.0 0.0\n");

    ASSERT_EQ(mol.num_atoms(), 1);
    EXPECT_EQ(mol.atom(0).Z, 2);
    EXPECT_EQ(mol.name(), "helium");
    EXPECT_EQ(mol.num_bonds(), 0);
}

TEST(XyzIoTest, ReadXyzStringThrowsOnMalformedInput) {
    EXPECT_THROW(sbox::io::read_xyz_string("2\nbad\nH 0.0 0.0 0.0\n"), std::runtime_error);
    EXPECT_THROW(sbox::io::read_xyz_string("1\nbad\nH 0.0 0.0\n"), std::runtime_error);
}
