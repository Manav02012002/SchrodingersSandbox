#include "io/fchk_io.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

namespace {

std::string minimal_h2_fchk() {
    return
        "H2 minimal test\n"
        "RHF STO-3G\n"
        "Number of atoms                            I              2\n"
        "Charge                                     I              0\n"
        "Multiplicity                               I              1\n"
        "Number of basis functions                  I              2\n"
        "Atomic numbers                             I   N=2\n"
        " 1 1\n"
        "Current cartesian coordinates              R   N=6\n"
        " 0.000000000000E+00 0.000000000000E+00 -7.000000000000E-01 0.000000000000E+00 0.000000000000E+00 7.000000000000E-01\n"
        "Total Energy                               R    -1.100000000000E+00\n"
        "Alpha Orbital Energies                     R   N=2\n"
        " -5.000000000000E-01 2.000000000000E-01\n"
        "Alpha MO coefficients                      R   N=4\n"
        " 7.071067810000E-01 7.071067810000E-01 7.071067810000E-01 -7.071067810000E-01\n";
}

std::string minimal_h2_fchk_with_basis_metadata() {
    return
        "H2 minimal with basis metadata\n"
        "RHF STO-3G\n"
        "Number of atoms                            I              2\n"
        "Charge                                     I              0\n"
        "Multiplicity                               I              1\n"
        "Number of basis functions                  I              2\n"
        "Atomic numbers                             I   N=2\n"
        " 1 1\n"
        "Shell types                                I   N=2\n"
        " 0 0\n"
        "Number of primitives per shell             I   N=2\n"
        " 3 3\n"
        "Shell to atom map                          I   N=2\n"
        " 1 2\n"
        "Primitive exponents                        R   N=6\n"
        " 3.425250910000E+00 6.239137300000E-01 1.688554000000E-01 3.425250910000E+00 6.239137300000E-01 1.688554000000E-01\n"
        "Contraction coefficients                   R   N=6\n"
        " 1.543289700000E-01 5.353281400000E-01 4.446345400000E-01 1.543289700000E-01 5.353281400000E-01 4.446345400000E-01\n"
        "Current cartesian coordinates              R   N=6\n"
        " 0.000000000000E+00 0.000000000000E+00 -7.000000000000E-01 0.000000000000E+00 0.000000000000E+00 7.000000000000E-01\n"
        "Total Energy                               R    -1.100000000000E+00\n"
        "Alpha Orbital Energies                     R   N=2\n"
        " -5.000000000000E-01 2.000000000000E-01\n"
        "Alpha MO coefficients                      R   N=4\n"
        " 7.071067810000E-01 7.071067810000E-01 7.071067810000E-01 -7.071067810000E-01\n";
}

std::filesystem::path temp_fchk_path(const char* stem) {
    return std::filesystem::temp_directory_path()
        / (std::string("schrodingerssandbox_") + stem + ".fchk");
}

}  // namespace

TEST(FchkIoTest, ParsesMinimalH2File) {
    const std::filesystem::path path = temp_fchk_path("fchk_minimal");
    {
        std::ofstream out(path);
        out << minimal_h2_fchk();
    }

    const sbox::io::FchkData data = sbox::io::read_fchk(path.string());
    std::filesystem::remove(path);

    EXPECT_EQ(data.num_atoms, 2);
    ASSERT_EQ(data.atomic_numbers.size(), 2U);
    EXPECT_EQ(data.atomic_numbers[0], 1);
    EXPECT_EQ(data.atomic_numbers[1], 1);
    EXPECT_NEAR(data.total_energy, -1.1, 1e-12);
    ASSERT_EQ(data.mo_energies.size(), 2);
    ASSERT_EQ(data.mo_coefficients.rows(), 2);
    ASSERT_EQ(data.mo_coefficients.cols(), 2);
}

TEST(FchkIoTest, MissingOptionalKeysDoNotThrow) {
    const std::filesystem::path path = temp_fchk_path("fchk_optional");
    {
        std::ofstream out(path);
        out << minimal_h2_fchk();
    }

    const sbox::io::FchkData data = sbox::io::read_fchk(path.string());
    std::filesystem::remove(path);

    EXPECT_TRUE(data.mulliken_charges.empty());
    EXPECT_NEAR(data.dipole_moment.norm(), 0.0, 1e-12);
}

TEST(FchkIoTest, FillsOccupationsForClosedShellH2) {
    const std::filesystem::path path = temp_fchk_path("fchk_occupations");
    {
        std::ofstream out(path);
        out << minimal_h2_fchk();
    }

    const sbox::io::FchkData data = sbox::io::read_fchk(path.string());
    std::filesystem::remove(path);

    ASSERT_EQ(data.occupations.size(), 2);
    EXPECT_NEAR(data.occupations(0), 2.0, 1e-12);
    EXPECT_NEAR(data.occupations(1), 0.0, 1e-12);
}

TEST(FchkIoTest, ParsesBasisMetadataWhenPresent) {
    const std::filesystem::path path = temp_fchk_path("fchk_basis_metadata");
    {
        std::ofstream out(path);
        out << minimal_h2_fchk_with_basis_metadata();
    }

    const sbox::io::FchkData data = sbox::io::read_fchk(path.string());
    std::filesystem::remove(path);

    ASSERT_EQ(data.shell_types.size(), 2U);
    ASSERT_EQ(data.primitives_per_shell.size(), 2U);
    ASSERT_EQ(data.shell_to_atom_map.size(), 2U);
    ASSERT_EQ(data.primitive_exponents.size(), 6U);
    ASSERT_EQ(data.contraction_coefficients.size(), 6U);
    EXPECT_EQ(data.shell_types[0], 0);
    EXPECT_EQ(data.shell_to_atom_map[1], 2);
    EXPECT_EQ(data.primitives_per_shell[0], 3);
    EXPECT_NEAR(data.primitive_exponents[0], 3.42525091, 1e-12);
    EXPECT_NEAR(data.contraction_coefficients[1], 0.53532814, 1e-12);
}
