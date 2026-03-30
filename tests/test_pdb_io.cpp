#include "io/pdb_io.h"

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

namespace {

std::filesystem::path temp_pdb_path(const char* name) {
    return std::filesystem::temp_directory_path() / name;
}

std::string minimal_peptide_pdb() {
    return
        "TITLE     Minimal peptide\n"
        "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       3.500   3.000   3.000  1.00  0.00           O\n"
        "ATOM      5  N   GLY A   2       4.000   2.000   3.000  1.00  0.00           N\n"
        "ATOM      6  CA  GLY A   2       5.000   2.000   3.000  1.00  0.00           C\n"
        "ATOM      7  C   GLY A   2       6.000   2.000   3.000  1.00  0.00           C\n"
        "ATOM      8  O   GLY A   2       6.500   3.000   3.000  1.00  0.00           O\n"
        "END\n";
}

void write_text_file(const std::filesystem::path& path, const std::string& content) {
    std::ofstream out(path);
    ASSERT_TRUE(static_cast<bool>(out));
    out << content;
}

}  // namespace

TEST(PdbIoTest, ReadsMinimalPeptideAndGroupsResidues) {
    const std::filesystem::path path = temp_pdb_path("schrodingerssandbox_test_peptide.pdb");
    write_text_file(path, minimal_peptide_pdb());

    const sbox::io::PDBData data = sbox::io::read_pdb(path.string());
    std::filesystem::remove(path);

    ASSERT_EQ(data.atoms.size(), 8u);
    ASSERT_EQ(data.residues.size(), 2u);
    ASSERT_EQ(data.chains.size(), 1u);
    EXPECT_EQ(data.atoms[0].Z, 7);
    EXPECT_EQ(data.atoms[1].Z, 6);
    EXPECT_EQ(data.residues[0].name, "ALA");
    EXPECT_EQ(data.residues[1].name, "GLY");
    EXPECT_EQ(data.chains[0].id, "A");
}

TEST(PdbIoTest, ToMolecularSystemPerceivesBondsWithoutConect) {
    const std::filesystem::path path = temp_pdb_path("schrodingerssandbox_test_peptide_bonds.pdb");
    write_text_file(path, minimal_peptide_pdb());

    const sbox::io::PDBData data = sbox::io::read_pdb(path.string());
    std::filesystem::remove(path);

    const sbox::chem::MolecularSystem mol = data.to_molecular_system();
    ASSERT_EQ(mol.num_atoms(), 8);
    EXPECT_GT(mol.num_bonds(), 0);
}

TEST(PdbIoTest, ParsesHetatmWater) {
    const std::filesystem::path path = temp_pdb_path("schrodingerssandbox_test_water.pdb");
    write_text_file(path,
                    "TITLE     Water\n"
                    "HETATM    1  O   HOH B   7       0.000   0.000   0.000  1.00 10.00           O\n"
                    "HETATM    2  H1  HOH B   7       0.758   0.000   0.504  1.00 10.00           H\n"
                    "HETATM    3  H2  HOH B   7      -0.758   0.000   0.504  1.00 10.00           H\n"
                    "END\n");

    const sbox::io::PDBData data = sbox::io::read_pdb(path.string());
    std::filesystem::remove(path);

    ASSERT_EQ(data.atoms.size(), 3u);
    EXPECT_EQ(data.residues.size(), 1u);
    EXPECT_EQ(data.residues[0].name, "HOH");
    EXPECT_EQ(data.atoms[0].Z, 8);
    EXPECT_EQ(data.atoms[1].Z, 1);
}

TEST(PdbIoTest, UsesConectRecordsWhenPresent) {
    const std::filesystem::path path = temp_pdb_path("schrodingerssandbox_test_conect.pdb");
    write_text_file(path,
                    "TITLE     Connected fragment\n"
                    "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N\n"
                    "ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C\n"
                    "ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C\n"
                    "CONECT    1    2\n"
                    "CONECT    2    1    3\n"
                    "CONECT    3    2\n"
                    "END\n");

    const sbox::io::PDBData data = sbox::io::read_pdb(path.string());
    std::filesystem::remove(path);

    ASSERT_EQ(data.conect_bonds.size(), 4u);
    const sbox::chem::MolecularSystem mol = data.to_molecular_system();
    ASSERT_EQ(mol.num_atoms(), 3);
    EXPECT_EQ(mol.num_bonds(), 2);
    EXPECT_TRUE(mol.has_bond(0, 1));
    EXPECT_TRUE(mol.has_bond(1, 2));
}
