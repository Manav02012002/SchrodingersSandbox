#include "io/trajectory_io.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {

constexpr double kBohrToAngstrom = 0.529177;

std::string temp_path(const char* stem) {
    static int counter = 0;
    const auto path = std::filesystem::temp_directory_path() /
                      ("sbox_traj_" + std::string(stem) + "_" + std::to_string(counter++) + ".xyz");
    return path.string();
}

void write_text_file(const std::string& path, const std::string& content) {
    std::ofstream out(path);
    ASSERT_TRUE(static_cast<bool>(out));
    out << content;
}

TEST(TrajectoryIO, ReadsThreeFrameTrajectory) {
    const std::string content =
        "2\n"
        "E = -1.117\n"
        "H 0.0 0.0 -0.37\n"
        "H 0.0 0.0 0.37\n"
        "2\n"
        "E = -1.100\n"
        "H 0.0 0.0 -0.40\n"
        "H 0.0 0.0 0.40\n"
        "2\n"
        "E = -1.050\n"
        "H 0.0 0.0 -0.50\n"
        "H 0.0 0.0 0.50\n";

    const std::string path = temp_path("traj_read");
    write_text_file(path, content);

    const auto traj = sbox::io::read_trajectory_xyz(path);
    EXPECT_EQ(traj.num_frames(), 3);
    EXPECT_NEAR(traj.frames[0].energy, -1.117, 1.0e-12);
    EXPECT_NEAR(traj.frames[1].energy, -1.100, 1.0e-12);
    EXPECT_NEAR(traj.frames[2].energy, -1.050, 1.0e-12);
    EXPECT_EQ(traj.frames[0].geometry.num_atoms(), 2);
}

TEST(TrajectoryIO, RoundTripsTrajectory) {
    sbox::io::Trajectory traj;
    for (int i = 0; i < 2; ++i) {
        sbox::io::TrajectoryFrame frame;
        frame.frame_index = i;
        frame.energy = -1.0 - 0.1 * i;
        frame.comment = "Frame";
        frame.geometry.add_atom({1, Eigen::Vector3d(0.0, 0.0, (-0.3 - 0.05 * i) / kBohrToAngstrom), "", 0});
        frame.geometry.add_atom({1, Eigen::Vector3d(0.0, 0.0, (0.3 + 0.05 * i) / kBohrToAngstrom), "", 0});
        frame.geometry.perceive_bonds();
        traj.frames.push_back(frame);
    }

    const std::string path = temp_path("traj_roundtrip");
    sbox::io::write_trajectory_xyz(path, traj);
    const auto reread = sbox::io::read_trajectory_xyz(path);
    ASSERT_EQ(reread.num_frames(), 2);
    EXPECT_NEAR(reread.frames[1].geometry.atom(1).position.z(), traj.frames[1].geometry.atom(1).position.z(), 1.0e-9);
}

TEST(TrajectoryIO, InterpolatesMidpoint) {
    const std::string content =
        "2\n"
        "E = -1.117\n"
        "H 0.0 0.0 -0.37\n"
        "H 0.0 0.0 0.37\n"
        "2\n"
        "E = -1.100\n"
        "H 0.0 0.0 -0.40\n"
        "H 0.0 0.0 0.40\n";
    const std::string path = temp_path("traj_interp");
    write_text_file(path, content);
    const auto traj = sbox::io::read_trajectory_xyz(path);
    const auto interp = traj.interpolate(0.5);
    EXPECT_NEAR(interp.atom(0).position.z(), -0.385 / kBohrToAngstrom, 1.0e-6);
    EXPECT_NEAR(interp.atom(1).position.z(), 0.385 / kBohrToAngstrom, 1.0e-6);
}

TEST(TrajectoryIO, InterpolatesExactFrameAtZero) {
    const std::string content =
        "2\n"
        "E = -1.117\n"
        "H 0.0 0.0 -0.37\n"
        "H 0.0 0.0 0.37\n"
        "2\n"
        "E = -1.100\n"
        "H 0.0 0.0 -0.40\n"
        "H 0.0 0.0 0.40\n";
    const std::string path = temp_path("traj_exact");
    write_text_file(path, content);
    const auto traj = sbox::io::read_trajectory_xyz(path);
    const auto interp = traj.interpolate(0.0);
    EXPECT_NEAR(interp.atom(0).position.z(), traj.frames[0].geometry.atom(0).position.z(), 1.0e-12);
    EXPECT_NEAR(interp.atom(1).position.z(), traj.frames[0].geometry.atom(1).position.z(), 1.0e-12);
}

TEST(TrajectoryIO, ReturnsEnergyVector) {
    const std::string content =
        "2\n"
        "E = -1.117\n"
        "H 0.0 0.0 -0.37\n"
        "H 0.0 0.0 0.37\n"
        "2\n"
        "E = -1.100\n"
        "H 0.0 0.0 -0.40\n"
        "H 0.0 0.0 0.40\n"
        "2\n"
        "E = -1.050\n"
        "H 0.0 0.0 -0.50\n"
        "H 0.0 0.0 0.50\n";
    const std::string path = temp_path("traj_energy");
    write_text_file(path, content);
    const auto traj = sbox::io::read_trajectory_xyz(path);
    const auto energies = traj.energies();
    ASSERT_EQ(energies.size(), 3);
    EXPECT_NEAR(energies[0], -1.117, 1.0e-12);
    EXPECT_NEAR(energies[1], -1.100, 1.0e-12);
    EXPECT_NEAR(energies[2], -1.050, 1.0e-12);
}

TEST(TrajectoryIO, MissingEnergyDefaultsToZero) {
    const std::string content =
        "2\n"
        "No energy here\n"
        "H 0.0 0.0 -0.37\n"
        "H 0.0 0.0 0.37\n";
    const std::string path = temp_path("traj_no_energy");
    write_text_file(path, content);
    const auto traj = sbox::io::read_trajectory_xyz(path);
    ASSERT_EQ(traj.num_frames(), 1);
    EXPECT_DOUBLE_EQ(traj.frames[0].energy, 0.0);
}

}  // namespace
