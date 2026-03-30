#include "analysis/orbital_composition.h"
#include "core/basis_set.h"
#include "core/molecular_system.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

#include <cmath>

namespace {

sbox::basis::MOData make_h2_sto3g() {
    sbox::basis::MOData mo_data;
    mo_data.basis.spherical = true;
    mo_data.atomic_numbers = {1, 1};
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

sbox::chem::MolecularSystem make_h2_molecule() {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({1, Eigen::Vector3d(0.0, 0.0, -0.7), "", 0});
    mol.add_atom({1, Eigen::Vector3d(0.0, 0.0, 0.7), "", 0});
    mol.add_bond(0, 1, sbox::chem::BondOrder::Single);
    return mol;
}

}  // namespace

TEST(OrbitalCompositionTest, BondingMoIsSplitEvenlyAcrossAtoms) {
    const auto mo = make_h2_sto3g();
    const auto mol = make_h2_molecule();
    const auto composition = sbox::analysis::analyze_orbital_composition(mo, mol, 0);

    ASSERT_EQ(composition.atom_contributions.size(), 2U);
    EXPECT_NEAR(composition.atom_contributions[0].total_weight, 0.5, 1e-8);
    EXPECT_NEAR(composition.atom_contributions[1].total_weight, 0.5, 1e-8);
}

TEST(OrbitalCompositionTest, AntibondingMoIsAlsoSplitEvenlyAcrossAtoms) {
    const auto mo = make_h2_sto3g();
    const auto mol = make_h2_molecule();
    const auto composition = sbox::analysis::analyze_orbital_composition(mo, mol, 1);

    ASSERT_EQ(composition.atom_contributions.size(), 2U);
    EXPECT_NEAR(composition.atom_contributions[0].total_weight, 0.5, 1e-8);
    EXPECT_NEAR(composition.atom_contributions[1].total_weight, 0.5, 1e-8);
}

TEST(OrbitalCompositionTest, SummaryContainsHydrogenAndContributionsAreSorted) {
    const auto mo = make_h2_sto3g();
    const auto mol = make_h2_molecule();
    const auto composition = sbox::analysis::analyze_orbital_composition(mo, mol, 0);

    EXPECT_NE(composition.summary.find("H"), std::string::npos);
    ASSERT_GE(composition.atom_contributions.size(), 2U);
    EXPECT_GE(composition.atom_contributions[0].total_weight, composition.atom_contributions[1].total_weight);
}
