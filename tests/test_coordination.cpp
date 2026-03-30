#include "chem/coordination.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

namespace {

using sbox::chem::Atom;
using sbox::chem::BondOrder;
using sbox::chem::CoordinationGeometry;
using sbox::chem::LigandSpec;
using sbox::chem::MolecularSystem;

constexpr double kPi = 3.14159265358979323846;

double degrees(double radians) {
    return radians * 180.0 / kPi;
}

MolecularSystem make_water_ligand() {
    MolecularSystem ligand;
    const double oh = 1.81;
    const double half_angle = 52.25 * kPi / 180.0;
    ligand.add_atom({8, Eigen::Vector3d(0.0, 0.0, 0.0), "O", 0});
    ligand.add_atom({1, Eigen::Vector3d(oh * std::sin(half_angle), 0.0, oh * std::cos(half_angle)), "H", 0});
    ligand.add_atom({1, Eigen::Vector3d(-oh * std::sin(half_angle), 0.0, oh * std::cos(half_angle)), "H", 0});
    ligand.add_bond(0, 1, BondOrder::Single);
    ligand.add_bond(0, 2, BondOrder::Single);
    return ligand;
}

MolecularSystem make_ammonia_ligand() {
    MolecularSystem ligand;
    const double nh = 1.91;
    const double z = nh * 0.40;
    const double r = std::sqrt(nh * nh - z * z);
    ligand.add_atom({7, Eigen::Vector3d(0.0, 0.0, 0.0), "N", 0});
    for (int i = 0; i < 3; ++i) {
        const double a = 2.0 * kPi * static_cast<double>(i) / 3.0;
        ligand.add_atom({1, Eigen::Vector3d(r * std::cos(a), r * std::sin(a), z), "H", 0});
        ligand.add_bond(0, i + 1, BondOrder::Single);
    }
    return ligand;
}

}  // namespace

TEST(CoordinationTest, OctahedralTemplateHasExpectedDirections) {
    const auto& tmpl = sbox::chem::get_template(CoordinationGeometry::Octahedral);
    ASSERT_EQ(tmpl.coordination_number, 6);
    ASSERT_EQ(tmpl.directions.size(), 6u);

    for (const auto& dir : tmpl.directions) {
        EXPECT_NEAR(dir.norm(), 1.0, 1.0e-9);
    }

    int opposite_pairs = 0;
    for (std::size_t i = 0; i < tmpl.directions.size(); ++i) {
        for (std::size_t j = i + 1; j < tmpl.directions.size(); ++j) {
            const double dot = tmpl.directions[i].dot(tmpl.directions[j]);
            if (std::abs(dot + 1.0) < 1.0e-9) {
                ++opposite_pairs;
            } else {
                EXPECT_NEAR(dot, 0.0, 1.0e-9);
            }
        }
    }
    EXPECT_EQ(opposite_pairs, 3);
}

TEST(CoordinationTest, TetrahedralTemplateAnglesAreCorrect) {
    const auto& tmpl = sbox::chem::get_template(CoordinationGeometry::Tetrahedral);
    ASSERT_EQ(tmpl.directions.size(), 4u);

    for (std::size_t i = 0; i < tmpl.directions.size(); ++i) {
        for (std::size_t j = i + 1; j < tmpl.directions.size(); ++j) {
            const double angle_deg = degrees(std::acos(std::clamp(tmpl.directions[i].dot(tmpl.directions[j]), -1.0, 1.0)));
            EXPECT_NEAR(angle_deg, 109.4712, 1.0e-3);
        }
    }
}

TEST(CoordinationTest, AssembleFeHexaaquaComplex) {
    LigandSpec water;
    water.donor_Z = 8;
    water.ligand = make_water_ligand();
    water.name = "H2O";

    std::vector<LigandSpec> ligands(6, water);
    MolecularSystem complex = sbox::chem::assemble_complex(26, 2, CoordinationGeometry::Octahedral, ligands);

    EXPECT_EQ(complex.num_atoms(), 19);
    EXPECT_EQ(complex.num_bonds(), 18);
    EXPECT_EQ(complex.charge(), 2);

    std::vector<double> distances;
    for (int i = 1; i < complex.num_atoms(); i += 3) {
        distances.push_back((complex.atom(i).position - complex.atom(0).position).norm());
    }
    ASSERT_EQ(distances.size(), 6u);
    for (double d : distances) {
        EXPECT_NEAR(d, distances.front(), 1.0e-2);
        EXPECT_NEAR(d, 3.97, 1.0e-2);
    }
}

TEST(CoordinationTest, AssembleNickelTetrammineComplex) {
    LigandSpec ammonia;
    ammonia.donor_Z = 7;
    ammonia.ligand = make_ammonia_ligand();
    ammonia.name = "NH3";

    std::vector<LigandSpec> ligands(4, ammonia);
    MolecularSystem complex = sbox::chem::assemble_complex(28, 2, CoordinationGeometry::Tetrahedral, ligands);

    EXPECT_EQ(complex.num_atoms(), 17);
    EXPECT_EQ(complex.num_bonds(), 16);
    EXPECT_EQ(complex.charge(), 2);

    for (int i = 1; i < complex.num_atoms(); i += 4) {
        EXPECT_NEAR((complex.atom(i).position - complex.atom(0).position).norm(), 3.59, 1.0e-2);
    }
}
