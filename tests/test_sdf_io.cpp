#include "core/molecular_system.h"
#include "io/sdf_io.h"

#include <Eigen/Core>

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <tuple>
#include <vector>

namespace {

using BondSignature = std::tuple<int, int, int>;

std::vector<BondSignature> bond_signatures(const sbox::chem::MolecularSystem& mol) {
    std::vector<BondSignature> out;
    for (const sbox::chem::Bond& bond : mol.bonds()) {
        const int i = std::min(bond.atom_i, bond.atom_j);
        const int j = std::max(bond.atom_i, bond.atom_j);
        out.emplace_back(i, j, static_cast<int>(bond.order));
    }
    std::sort(out.begin(), out.end());
    return out;
}

sbox::chem::MolecularSystem make_formaldehyde() {
    sbox::chem::MolecularSystem mol;
    mol.set_name("formaldehyde");
    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "C"});
    mol.add_atom({8, Eigen::Vector3d(2.30, 0.0, 0.0), "O"});
    mol.add_atom({1, Eigen::Vector3d(-1.00, 1.60, 0.0), "H1"});
    mol.add_atom({1, Eigen::Vector3d(-1.00, -1.60, 0.0), "H2"});
    mol.add_bond(0, 1, sbox::chem::BondOrder::Double);
    mol.add_bond(0, 2, sbox::chem::BondOrder::Single);
    mol.add_bond(0, 3, sbox::chem::BondOrder::Single);
    return mol;
}

std::filesystem::path temp_sdf_path() {
    return std::filesystem::temp_directory_path() / "schrodingerssandbox_formaldehyde.sdf";
}

}  // namespace

TEST(SdfIoTest, ReadSdfStringParsesEthanol) {
    const std::string sdf =
        "ethanol\n"
        "Codex\n"
        "\n"
        "  9  8  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.5200    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    2.8800    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.5400    0.9300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.5400   -0.9300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.0000    0.0000    1.0900 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.5200    0.0000   -1.0900 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.9800    0.9300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    3.2400    0.7600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  2  3  1  0  0  0  0\n"
        "  1  4  1  0  0  0  0\n"
        "  1  5  1  0  0  0  0\n"
        "  1  6  1  0  0  0  0\n"
        "  2  7  1  0  0  0  0\n"
        "  2  8  1  0  0  0  0\n"
        "  3  9  1  0  0  0  0\n"
        "M  END\n";

    const sbox::chem::MolecularSystem mol = sbox::io::read_sdf_string(sdf);

    ASSERT_EQ(mol.num_atoms(), 9);
    ASSERT_EQ(mol.num_bonds(), 8);
    EXPECT_EQ(mol.atom(0).Z, 6);
    EXPECT_EQ(mol.atom(1).Z, 6);
    EXPECT_EQ(mol.atom(2).Z, 8);
    for (int i = 3; i < 9; ++i) {
        EXPECT_EQ(mol.atom(i).Z, 1);
    }
    for (const sbox::chem::Bond& bond : mol.bonds()) {
        EXPECT_EQ(bond.order, sbox::chem::BondOrder::Single);
    }
}

TEST(SdfIoTest, RoundTripPreservesAtomsAndBonds) {
    const sbox::chem::MolecularSystem original = make_formaldehyde();
    const std::filesystem::path path = temp_sdf_path();

    sbox::io::write_sdf(path.string(), original);
    const sbox::chem::MolecularSystem round_tripped = sbox::io::read_sdf(path.string());
    std::filesystem::remove(path);

    ASSERT_EQ(round_tripped.num_atoms(), original.num_atoms());
    ASSERT_EQ(round_tripped.num_bonds(), original.num_bonds());
    for (int i = 0; i < original.num_atoms(); ++i) {
        EXPECT_EQ(round_tripped.atom(i).Z, original.atom(i).Z);
        EXPECT_NEAR(round_tripped.atom(i).position.x(), original.atom(i).position.x(), 1e-4);
        EXPECT_NEAR(round_tripped.atom(i).position.y(), original.atom(i).position.y(), 1e-4);
        EXPECT_NEAR(round_tripped.atom(i).position.z(), original.atom(i).position.z(), 1e-4);
    }
    EXPECT_EQ(bond_signatures(round_tripped), bond_signatures(original));
}

TEST(SdfIoTest, BondOrdersSurviveRoundTrip) {
    const sbox::chem::MolecularSystem original = make_formaldehyde();
    const std::filesystem::path path = temp_sdf_path();

    sbox::io::write_sdf(path.string(), original);
    const sbox::chem::MolecularSystem round_tripped = sbox::io::read_sdf(path.string());
    std::filesystem::remove(path);

    ASSERT_EQ(round_tripped.num_bonds(), 3);
    bool saw_double = false;
    for (const sbox::chem::Bond& bond : round_tripped.bonds()) {
        if ((bond.atom_i == 0 && bond.atom_j == 1) || (bond.atom_i == 1 && bond.atom_j == 0)) {
            EXPECT_EQ(bond.order, sbox::chem::BondOrder::Double);
            saw_double = true;
        }
    }
    EXPECT_TRUE(saw_double);
}

TEST(SdfIoTest, ReadSdfStringParsesAromaticBenzene) {
    const std::string sdf =
        "benzene\n"
        "Codex\n"
        "\n"
        "  6  6  0  0  0  0  0  0  0999 V2000\n"
        "    1.3960    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.6980    1.2090    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.6980    1.2090    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -1.3960    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.6980   -1.2090    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.6980   -1.2090    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  4  0  0  0  0\n"
        "  2  3  4  0  0  0  0\n"
        "  3  4  4  0  0  0  0\n"
        "  4  5  4  0  0  0  0\n"
        "  5  6  4  0  0  0  0\n"
        "  6  1  4  0  0  0  0\n"
        "M  END\n";

    const sbox::chem::MolecularSystem mol = sbox::io::read_sdf_string(sdf);

    ASSERT_EQ(mol.num_bonds(), 6);
    for (const sbox::chem::Bond& bond : mol.bonds()) {
        EXPECT_EQ(bond.order, sbox::chem::BondOrder::Aromatic);
    }
}
