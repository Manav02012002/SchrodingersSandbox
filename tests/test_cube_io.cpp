#include "io/cube_io.h"

#include <Eigen/Core>

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

namespace {

sbox::io::CubeData make_small_cube() {
    sbox::io::CubeData cube;
    cube.comment1 = "small cube";
    cube.comment2 = "round trip";
    cube.atom_Z = {1};
    cube.atom_pos = {Eigen::Vector3d(0.0, 0.0, 0.0)};
    cube.origin = Eigen::Vector3d(-1.0, -1.0, -1.0);
    cube.step_x = Eigen::Vector3d(0.25, 0.0, 0.0);
    cube.step_y = Eigen::Vector3d(0.0, 0.25, 0.0);
    cube.step_z = Eigen::Vector3d(0.0, 0.0, 0.25);
    cube.nx = 4;
    cube.ny = 4;
    cube.nz = 4;
    cube.data.resize(64);
    for (std::size_t i = 0; i < cube.data.size(); ++i) {
        cube.data[i] = static_cast<float>(i) * 0.001f;
    }
    return cube;
}

std::filesystem::path small_cube_path() {
    return std::filesystem::temp_directory_path() / "schrodingerssandbox_small.cube";
}

std::filesystem::path h2_cube_path() {
    return std::filesystem::temp_directory_path() / "schrodingerssandbox_h2.cube";
}

}  // namespace

TEST(CubeIoTest, RoundTripPreservesGridAndData) {
    const sbox::io::CubeData original = make_small_cube();
    const std::filesystem::path path = small_cube_path();

    sbox::io::write_cube(path.string(), original);
    const sbox::io::CubeData round_tripped = sbox::io::read_cube(path.string());
    std::filesystem::remove(path);

    EXPECT_EQ(round_tripped.nx, original.nx);
    EXPECT_EQ(round_tripped.ny, original.ny);
    EXPECT_EQ(round_tripped.nz, original.nz);
    EXPECT_NEAR(round_tripped.origin.x(), original.origin.x(), 1e-5);
    EXPECT_NEAR(round_tripped.origin.y(), original.origin.y(), 1e-5);
    EXPECT_NEAR(round_tripped.origin.z(), original.origin.z(), 1e-5);
    EXPECT_NEAR(round_tripped.step_x.x(), original.step_x.x(), 1e-5);
    EXPECT_NEAR(round_tripped.step_y.y(), original.step_y.y(), 1e-5);
    EXPECT_NEAR(round_tripped.step_z.z(), original.step_z.z(), 1e-5);
    ASSERT_EQ(round_tripped.atom_Z.size(), 1U);
    EXPECT_EQ(round_tripped.atom_Z[0], 1);
    EXPECT_NEAR(round_tripped.atom_pos[0].norm(), 0.0, 1e-5);
    ASSERT_EQ(round_tripped.data.size(), original.data.size());
    for (std::size_t i = 0; i < original.data.size(); ++i) {
        EXPECT_NEAR(round_tripped.data[i], original.data[i], 1e-5);
    }
}

TEST(CubeIoTest, RoundTripPreservesTwoAtomPositions) {
    sbox::io::CubeData cube;
    cube.comment1 = "h2";
    cube.comment2 = "two atom test";
    cube.atom_Z = {1, 1};
    cube.atom_pos = {
        Eigen::Vector3d(-0.7, 0.0, 0.0),
        Eigen::Vector3d(0.7, 0.0, 0.0),
    };
    cube.origin = Eigen::Vector3d(-2.0, -2.0, -2.0);
    cube.step_x = Eigen::Vector3d(0.5, 0.0, 0.0);
    cube.step_y = Eigen::Vector3d(0.0, 0.5, 0.0);
    cube.step_z = Eigen::Vector3d(0.0, 0.0, 0.5);
    cube.nx = 8;
    cube.ny = 8;
    cube.nz = 8;
    cube.data.assign(static_cast<std::size_t>(cube.nx * cube.ny * cube.nz), 0.0f);

    const std::filesystem::path path = h2_cube_path();
    sbox::io::write_cube(path.string(), cube);
    const sbox::io::CubeData round_tripped = sbox::io::read_cube(path.string());
    std::filesystem::remove(path);

    ASSERT_EQ(round_tripped.atom_pos.size(), 2U);
    EXPECT_NEAR(round_tripped.atom_pos[0].x(), -0.7, 1e-5);
    EXPECT_NEAR(round_tripped.atom_pos[1].x(), 0.7, 1e-5);
    EXPECT_NEAR(round_tripped.atom_pos[0].y(), 0.0, 1e-5);
    EXPECT_NEAR(round_tripped.atom_pos[1].z(), 0.0, 1e-5);
}

TEST(CubeIoTest, ReadMalformedCubeThrows) {
    const std::filesystem::path path = std::filesystem::temp_directory_path() / "schrodingerssandbox_bad.cube";
    {
        std::ofstream out(path);
        out << "comment 1\n";
        out << "comment 2\n";
        out << "    1  0.00000e+00  0.00000e+00  0.00000e+00\n";
        out << "    2  1.00000e+00  0.00000e+00  0.00000e+00\n";
        out << "    2  0.00000e+00  1.00000e+00  0.00000e+00\n";
        out << "    2  0.00000e+00  0.00000e+00  1.00000e+00\n";
        out << "    1  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00\n";
        out << " 1.00000e-01 2.00000e-01\n";
    }

    EXPECT_THROW(sbox::io::read_cube(path.string()), std::runtime_error);
    std::filesystem::remove(path);
}
