#include "core/zmatrix.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

double degrees_to_radians(double degrees) {
    return degrees * kPi / 180.0;
}

}  // namespace

TEST(ZMatrixTest, BuildsH2) {
    const std::vector<sbox::chem::ZMatrixEntry> zmat = {
        {1},
        {1, 0, 1.4},
    };

    const sbox::chem::MolecularSystem mol = sbox::chem::zmatrix_to_cartesian(zmat);

    ASSERT_EQ(mol.num_atoms(), 2);
    EXPECT_NEAR(mol.distance(0, 1), 1.4, 1e-12);
}

TEST(ZMatrixTest, BuildsWater) {
    const std::vector<sbox::chem::ZMatrixEntry> zmat = {
        {8},
        {1, 0, 1.81},
        {1, 0, 1.81, 1, degrees_to_radians(104.5)},
    };

    const sbox::chem::MolecularSystem mol = sbox::chem::zmatrix_to_cartesian(zmat);

    ASSERT_EQ(mol.num_atoms(), 3);
    EXPECT_NEAR(mol.distance(0, 1), 1.81, 1e-12);
    EXPECT_NEAR(mol.distance(0, 2), 1.81, 1e-12);
    EXPECT_NEAR(mol.angle(1, 0, 2), degrees_to_radians(104.5), 1e-6);
}

TEST(ZMatrixTest, BuildsHydrogenPeroxideDihedral) {
    const double hoo_angle = degrees_to_radians(101.9);
    const double dihedral = degrees_to_radians(111.5);
    const std::vector<sbox::chem::ZMatrixEntry> zmat = {
        {1},
        {8, 0, 1.81},
        {8, 1, 2.80, 0, hoo_angle},
        {1, 2, 1.81, 1, hoo_angle, 0, dihedral},
    };

    const sbox::chem::MolecularSystem mol = sbox::chem::zmatrix_to_cartesian(zmat);

    ASSERT_EQ(mol.num_atoms(), 4);
    EXPECT_NEAR(std::abs(mol.dihedral(0, 1, 2, 3)), dihedral, 1e-4);
}

TEST(ZMatrixTest, BuildsMethaneWithTetrahedralAngles) {
    const double tetra = degrees_to_radians(109.47);
    const std::vector<sbox::chem::ZMatrixEntry> zmat = {
        {6},
        {1, 0, 2.06},
        {1, 0, 2.06, 1, tetra},
        {1, 0, 2.06, 1, tetra, 2, degrees_to_radians(120.0)},
        {1, 0, 2.06, 1, tetra, 2, degrees_to_radians(-120.0)},
    };

    const sbox::chem::MolecularSystem mol = sbox::chem::zmatrix_to_cartesian(zmat);

    ASSERT_EQ(mol.num_atoms(), 5);
    for (int i = 1; i <= 4; ++i) {
        EXPECT_NEAR(mol.distance(0, i), 2.06, 1e-10);
    }

    for (int i = 1; i <= 4; ++i) {
        for (int j = i + 1; j <= 4; ++j) {
            EXPECT_NEAR(mol.angle(i, 0, j), tetra, 1e-4);
        }
    }
}
