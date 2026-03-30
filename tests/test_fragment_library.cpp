#include "editor/fragment_library.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

TEST(FragmentLibraryTest, HasAtLeastTwentyFragments) {
    sbox::editor::FragmentLibrary library;
    EXPECT_GE(static_cast<int>(library.all().size()), 20);
}

TEST(FragmentLibraryTest, HasMultipleCategories) {
    sbox::editor::FragmentLibrary library;
    EXPECT_GE(static_cast<int>(library.categories().size()), 4);
}

TEST(FragmentLibraryTest, BenzeneHasExpectedCounts) {
    sbox::editor::FragmentLibrary library;
    const sbox::editor::Fragment* benzene = library.find("benzene");
    ASSERT_NE(benzene, nullptr);
    EXPECT_EQ(benzene->molecule.num_atoms(), 12);
    EXPECT_EQ(benzene->molecule.num_bonds(), 12);
}

TEST(FragmentLibraryTest, MethylHasAttachmentAtom) {
    sbox::editor::FragmentLibrary library;
    const sbox::editor::Fragment* methyl = library.find("methyl");
    ASSERT_NE(methyl, nullptr);
    EXPECT_GE(methyl->attachment_atom, 0);
}

TEST(FragmentLibraryTest, PlaceMovesFragmentAnchorToPosition) {
    sbox::editor::FragmentLibrary library;
    const sbox::editor::Fragment* methyl = library.find("methyl");
    ASSERT_NE(methyl, nullptr);
    const Eigen::Vector3d target(4.0, -2.0, 1.0);
    const sbox::chem::MolecularSystem placed = library.place(*methyl, target);
    ASSERT_GE(methyl->attachment_atom, 0);
    EXPECT_NEAR((placed.atom(methyl->attachment_atom).position - target).norm(), 0.0, 1e-9);
}

TEST(FragmentLibraryTest, AllFragmentsHaveValidMolecules) {
    sbox::editor::FragmentLibrary library;
    for (const auto& fragment : library.all()) {
        EXPECT_GT(fragment.molecule.num_atoms(), 0) << fragment.name;
        for (const auto& bond : fragment.molecule.bonds()) {
            EXPECT_GE(bond.atom_i, 0) << fragment.name;
            EXPECT_GE(bond.atom_j, 0) << fragment.name;
            EXPECT_LT(bond.atom_i, fragment.molecule.num_atoms()) << fragment.name;
            EXPECT_LT(bond.atom_j, fragment.molecule.num_atoms()) << fragment.name;
        }
    }
}
