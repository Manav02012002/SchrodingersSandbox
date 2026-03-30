#include "chem/ligand_library.h"

#include <gtest/gtest.h>

namespace {

using sbox::chem::LigandDenticity;
using sbox::chem::LigandLibrary;

}  // namespace

TEST(LigandLibraryTest, HasAtLeastTwentyLigands) {
    LigandLibrary library;
    EXPECT_GE(static_cast<int>(library.all().size()), 20);
}

TEST(LigandLibraryTest, WaterPropertiesAreCorrect) {
    LigandLibrary library;
    const auto* water = library.find("water");
    ASSERT_NE(water, nullptr);
    EXPECT_EQ(water->molecule.num_atoms(), 3);
    ASSERT_EQ(water->donor_atoms.size(), 1u);
    EXPECT_EQ(water->charge, 0);
    EXPECT_EQ(water->denticity, LigandDenticity::Monodentate);
}

TEST(LigandLibraryTest, EnIsBidentate) {
    LigandLibrary library;
    const auto* en = library.find("en");
    ASSERT_NE(en, nullptr);
    EXPECT_EQ(en->denticity, LigandDenticity::Bidentate);
    EXPECT_EQ(en->donor_atoms.size(), 2u);
}

TEST(LigandLibraryTest, EdtaIsHexadentate) {
    LigandLibrary library;
    const auto* edta = library.find("edta");
    ASSERT_NE(edta, nullptr);
    EXPECT_EQ(edta->denticity, LigandDenticity::Hexadentate);
    EXPECT_EQ(edta->donor_atoms.size(), 6u);
}

TEST(LigandLibraryTest, CyanidePropertiesAreCorrect) {
    LigandLibrary library;
    const auto* cn = library.find("cn");
    ASSERT_NE(cn, nullptr);
    EXPECT_EQ(cn->charge, -1);
    EXPECT_EQ(cn->field_strength, "very strong");
}

TEST(LigandLibraryTest, AllLigandsHaveValidDonorAtoms) {
    LigandLibrary library;
    for (const auto& ligand : library.all()) {
        EXPECT_GT(ligand.molecule.num_atoms(), 0) << ligand.name;
        for (int donor : ligand.donor_atoms) {
            EXPECT_GE(donor, 0) << ligand.name;
            EXPECT_LT(donor, ligand.molecule.num_atoms()) << ligand.name;
        }
    }
}

TEST(LigandLibraryTest, ToLigandSpecUsesSelectedDonor) {
    LigandLibrary library;
    const auto* en = library.find("en");
    ASSERT_NE(en, nullptr);
    const auto spec = library.to_ligand_spec(*en, 1);
    EXPECT_EQ(spec.donor_Z, 7);
    ASSERT_GT(spec.ligand.num_atoms(), 0);
    EXPECT_EQ(spec.ligand.atom(0).Z, 7);
}
