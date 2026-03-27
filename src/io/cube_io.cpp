#include "io/cube_io.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace sbox::io {
namespace {

constexpr double kAngstromToBohr = 1.8897259886;

template <typename T>
T read_value(std::istream& input, const std::string& what) {
    T value{};
    if (!(input >> value)) {
        throw std::runtime_error("Malformed cube file: could not read " + what);
    }
    return value;
}

CubeData read_cube_stream(std::istream& input) {
    CubeData cube;

    if (!std::getline(input, cube.comment1)) {
        throw std::runtime_error("Malformed cube file: missing comment line 1");
    }
    if (!std::getline(input, cube.comment2)) {
        throw std::runtime_error("Malformed cube file: missing comment line 2");
    }

    int num_atoms_raw = read_value<int>(input, "atom count");
    cube.origin.x() = read_value<double>(input, "origin x");
    cube.origin.y() = read_value<double>(input, "origin y");
    cube.origin.z() = read_value<double>(input, "origin z");

    int nx_raw = read_value<int>(input, "nx");
    cube.step_x.x() = read_value<double>(input, "step_x x");
    cube.step_x.y() = read_value<double>(input, "step_x y");
    cube.step_x.z() = read_value<double>(input, "step_x z");

    int ny_raw = read_value<int>(input, "ny");
    cube.step_y.x() = read_value<double>(input, "step_y x");
    cube.step_y.y() = read_value<double>(input, "step_y y");
    cube.step_y.z() = read_value<double>(input, "step_y z");

    int nz_raw = read_value<int>(input, "nz");
    cube.step_z.x() = read_value<double>(input, "step_z x");
    cube.step_z.y() = read_value<double>(input, "step_z y");
    cube.step_z.z() = read_value<double>(input, "step_z z");

    const bool angstrom_units = (nx_raw < 0) || (ny_raw < 0) || (nz_raw < 0);
    cube.nx = std::abs(nx_raw);
    cube.ny = std::abs(ny_raw);
    cube.nz = std::abs(nz_raw);

    if (angstrom_units) {
        cube.origin *= kAngstromToBohr;
        cube.step_x *= kAngstromToBohr;
        cube.step_y *= kAngstromToBohr;
        cube.step_z *= kAngstromToBohr;
    }

    const int num_atoms = std::abs(num_atoms_raw);
    cube.atom_Z.reserve(static_cast<std::size_t>(num_atoms));
    cube.atom_pos.reserve(static_cast<std::size_t>(num_atoms));

    for (int i = 0; i < num_atoms; ++i) {
        const int z = read_value<int>(input, "atom Z");
        (void)read_value<double>(input, "atom charge");
        Eigen::Vector3d pos;
        pos.x() = read_value<double>(input, "atom x");
        pos.y() = read_value<double>(input, "atom y");
        pos.z() = read_value<double>(input, "atom z");
        if (angstrom_units) {
            pos *= kAngstromToBohr;
        }
        cube.atom_Z.push_back(z);
        cube.atom_pos.push_back(pos);
    }

    const std::size_t total_values =
        static_cast<std::size_t>(cube.nx) * static_cast<std::size_t>(cube.ny) * static_cast<std::size_t>(cube.nz);
    cube.data.reserve(total_values);

    float value = 0.0f;
    while (input >> value) {
        cube.data.push_back(value);
    }

    if (cube.data.size() != total_values) {
        throw std::runtime_error("Malformed cube file: expected " + std::to_string(total_values) +
                                 " volumetric values, got " + std::to_string(cube.data.size()));
    }

    return cube;
}

}  // namespace

CubeData read_cube(const std::string& filepath) {
    std::ifstream input(filepath);
    if (!input) {
        throw std::runtime_error("Could not open cube file: " + filepath);
    }
    return read_cube_stream(input);
}

void write_cube(const std::string& filepath, const CubeData& cube) {
    std::ofstream output(filepath);
    if (!output) {
        throw std::runtime_error("Could not open cube file for writing: " + filepath);
    }

    output << cube.comment1 << '\n';
    output << cube.comment2 << '\n';
    output << std::scientific << std::setprecision(5);
    output << std::setw(5) << cube.atom_Z.size()
           << std::setw(13) << cube.origin.x()
           << std::setw(13) << cube.origin.y()
           << std::setw(13) << cube.origin.z() << '\n';
    output << std::setw(5) << cube.nx
           << std::setw(13) << cube.step_x.x()
           << std::setw(13) << cube.step_x.y()
           << std::setw(13) << cube.step_x.z() << '\n';
    output << std::setw(5) << cube.ny
           << std::setw(13) << cube.step_y.x()
           << std::setw(13) << cube.step_y.y()
           << std::setw(13) << cube.step_y.z() << '\n';
    output << std::setw(5) << cube.nz
           << std::setw(13) << cube.step_z.x()
           << std::setw(13) << cube.step_z.y()
           << std::setw(13) << cube.step_z.z() << '\n';

    for (std::size_t i = 0; i < cube.atom_Z.size(); ++i) {
        output << std::setw(5) << cube.atom_Z[i]
               << std::setw(13) << 0.0
               << std::setw(13) << cube.atom_pos[i].x()
               << std::setw(13) << cube.atom_pos[i].y()
               << std::setw(13) << cube.atom_pos[i].z() << '\n';
    }

    for (std::size_t i = 0; i < cube.data.size(); ++i) {
        output << std::setw(13) << cube.data[i];
        if ((i + 1) % 6 == 0 || i + 1 == cube.data.size()) {
            output << '\n';
        }
    }
}

}  // namespace sbox::io
