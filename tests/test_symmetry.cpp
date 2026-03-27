#include "core/symmetry.h"

#include <Eigen/Core>

#include <gtest/gtest.h>

#include <cmath>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

double deg(double value) {
    return value * kPi / 180.0;
}

sbox::chem::MolecularSystem make_water() {
    sbox::chem::MolecularSystem mol;
    const double r = 1.81;
    const double theta = deg(52.25);
    mol.add_atom({8, Eigen::Vector3d(0.0, 0.0, 0.0), "O"});
    mol.add_atom({1, Eigen::Vector3d(r * std::sin(theta), 0.0, r * std::cos(theta)), "H1"});
    mol.add_atom({1, Eigen::Vector3d(-r * std::sin(theta), 0.0, r * std::cos(theta)), "H2"});
    return mol;
}

sbox::chem::MolecularSystem make_ammonia() {
    sbox::chem::MolecularSystem mol;
    const double r = 1.90;
    const double z = 0.7;
    const double rho = std::sqrt(r * r - z * z);
    mol.add_atom({7, Eigen::Vector3d(0.0, 0.0, 0.4), "N"});
    for (int i = 0; i < 3; ++i) {
        const double phi = 2.0 * kPi * static_cast<double>(i) / 3.0;
        mol.add_atom({1, Eigen::Vector3d(rho * std::cos(phi), rho * std::sin(phi), -0.3), "H"});
    }
    return mol;
}

sbox::chem::MolecularSystem make_methane() {
    sbox::chem::MolecularSystem mol;
    const double scale = 2.06 / std::sqrt(3.0);
    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "C"});
    mol.add_atom({1, scale * Eigen::Vector3d(1.0, 1.0, 1.0), "H1"});
    mol.add_atom({1, scale * Eigen::Vector3d(1.0, -1.0, -1.0), "H2"});
    mol.add_atom({1, scale * Eigen::Vector3d(-1.0, 1.0, -1.0), "H3"});
    mol.add_atom({1, scale * Eigen::Vector3d(-1.0, -1.0, 1.0), "H4"});
    return mol;
}

sbox::chem::MolecularSystem make_benzene() {
    sbox::chem::MolecularSystem mol;
    const double rc = 2.64;
    const double rh = 4.70;
    for (int i = 0; i < 6; ++i) {
        const double phi = 2.0 * kPi * static_cast<double>(i) / 6.0;
        mol.add_atom({6, Eigen::Vector3d(rc * std::cos(phi), rc * std::sin(phi), 0.0), "C"});
    }
    for (int i = 0; i < 6; ++i) {
        const double phi = 2.0 * kPi * static_cast<double>(i) / 6.0;
        mol.add_atom({1, Eigen::Vector3d(rh * std::cos(phi), rh * std::sin(phi), 0.0), "H"});
    }
    return mol;
}

sbox::chem::MolecularSystem make_h2() {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({1, Eigen::Vector3d(0.0, 0.0, -0.7), "H1"});
    mol.add_atom({1, Eigen::Vector3d(0.0, 0.0, 0.7), "H2"});
    return mol;
}

sbox::chem::MolecularSystem make_hcl() {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({1, Eigen::Vector3d(0.0, 0.0, -1.2), "H"});
    mol.add_atom({17, Eigen::Vector3d(0.0, 0.0, 1.0), "Cl"});
    return mol;
}

sbox::chem::MolecularSystem make_chfclbr() {
    sbox::chem::MolecularSystem mol;
    const double scale = 2.1 / std::sqrt(3.0);
    mol.add_atom({6, Eigen::Vector3d(0.1, -0.05, 0.0), "C"});
    mol.add_atom({1, Eigen::Vector3d(0.1, -0.05, 0.0) + scale * Eigen::Vector3d(1.0, 1.0, 1.0), "H"});
    mol.add_atom({9, Eigen::Vector3d(0.1, -0.05, 0.0) + scale * Eigen::Vector3d(1.0, -1.0, -1.0), "F"});
    mol.add_atom({17, Eigen::Vector3d(0.1, -0.05, 0.0) + scale * Eigen::Vector3d(-1.0, 1.0, -1.0), "Cl"});
    mol.add_atom({35, Eigen::Vector3d(0.1, -0.05, 0.0) + scale * Eigen::Vector3d(-1.0, -1.0, 1.0), "Br"});
    return mol;
}

sbox::chem::MolecularSystem make_trans_dichloroethene() {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({6, Eigen::Vector3d(-1.26, 0.0, 0.0), "C1"});
    mol.add_atom({6, Eigen::Vector3d(1.26, 0.0, 0.0), "C2"});
    mol.add_atom({1, Eigen::Vector3d(-2.20, 0.0, 1.70), "H1"});
    mol.add_atom({17, Eigen::Vector3d(-2.80, 0.0, -1.70), "Cl1"});
    mol.add_atom({17, Eigen::Vector3d(2.80, 0.0, 1.70), "Cl2"});
    mol.add_atom({1, Eigen::Vector3d(2.20, 0.0, -1.70), "H2"});
    return mol;
}

sbox::chem::MolecularSystem make_allene() {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "C2"});
    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, 2.5), "C3"});
    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, -2.5), "C1"});
    mol.add_atom({1, Eigen::Vector3d(1.8, 0.0, 3.4), "H1"});
    mol.add_atom({1, Eigen::Vector3d(-1.8, 0.0, 3.4), "H2"});
    mol.add_atom({1, Eigen::Vector3d(0.0, 1.8, -3.4), "H3"});
    mol.add_atom({1, Eigen::Vector3d(0.0, -1.8, -3.4), "H4"});
    return mol;
}

sbox::chem::MolecularSystem make_ethane_staggered() {
    sbox::chem::MolecularSystem mol;
    const double c = 1.45;
    const double rho = 1.78;
    const double hz = 1.02;
    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, c), "C1"});
    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, -c), "C2"});
    for (int i = 0; i < 3; ++i) {
        const double phi = 2.0 * kPi * static_cast<double>(i) / 3.0;
        mol.add_atom({1, Eigen::Vector3d(rho * std::cos(phi), rho * std::sin(phi), c + hz), "H"});
    }
    for (int i = 0; i < 3; ++i) {
        const double phi = 2.0 * kPi * static_cast<double>(i) / 3.0 + kPi / 3.0;
        mol.add_atom({1, Eigen::Vector3d(rho * std::cos(phi), rho * std::sin(phi), -c - hz), "H"});
    }
    return mol;
}

void expect_pg(const sbox::chem::MolecularSystem& mol, sbox::chem::PointGroup pg) {
    EXPECT_EQ(sbox::chem::detect_point_group(mol, 0.3), pg);
}

}  // namespace

TEST(SymmetryTest, WaterIsC2v) { expect_pg(make_water(), sbox::chem::PointGroup::C2v); }
TEST(SymmetryTest, AmmoniaIsC3v) { expect_pg(make_ammonia(), sbox::chem::PointGroup::C3v); }
TEST(SymmetryTest, MethaneIsTd) { expect_pg(make_methane(), sbox::chem::PointGroup::Td); }
TEST(SymmetryTest, BenzeneIsD6h) { expect_pg(make_benzene(), sbox::chem::PointGroup::D6h); }
TEST(SymmetryTest, H2IsDinfh) { expect_pg(make_h2(), sbox::chem::PointGroup::Dinfh); }
TEST(SymmetryTest, HClIsCinfv) { expect_pg(make_hcl(), sbox::chem::PointGroup::Cinfv); }
TEST(SymmetryTest, ChiralCarbonIsC1) { expect_pg(make_chfclbr(), sbox::chem::PointGroup::C1); }
TEST(SymmetryTest, TransDichloroetheneIsC2h) { expect_pg(make_trans_dichloroethene(), sbox::chem::PointGroup::C2h); }
TEST(SymmetryTest, AlleneIsD2d) { expect_pg(make_allene(), sbox::chem::PointGroup::D2d); }
TEST(SymmetryTest, EthaneStaggeredIsD3d) { expect_pg(make_ethane_staggered(), sbox::chem::PointGroup::D3d); }
