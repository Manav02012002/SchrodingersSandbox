#pragma once

#include <Eigen/Core>

#include <string>
#include <vector>

namespace sbox::io {

struct CubeData {
    std::string comment1;
    std::string comment2;
    std::vector<int> atom_Z;
    std::vector<Eigen::Vector3d> atom_pos;
    Eigen::Vector3d origin;
    Eigen::Vector3d step_x;
    Eigen::Vector3d step_y;
    Eigen::Vector3d step_z;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    std::vector<float> data;

    float at(int ix, int iy, int iz) const { return data[static_cast<std::size_t>(ix) * ny * nz + static_cast<std::size_t>(iy) * nz + static_cast<std::size_t>(iz)]; }
};

CubeData read_cube(const std::string& filepath);
void write_cube(const std::string& filepath, const CubeData& cube);

}  // namespace sbox::io
