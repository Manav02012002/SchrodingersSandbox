#include "io/trajectory_io.h"

#include "core/elements.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

namespace sbox::io {
namespace {

constexpr double kAngstromToBohr = 1.8897259886;
constexpr double kBohrToAngstrom = 1.0 / kAngstromToBohr;

int symbol_to_z(const std::string& symbol) {
    for (int z = 1; z <= 118; ++z) {
        if (symbol == sbox::elements::get_element(z).symbol) {
            return z;
        }
    }
    for (int z = 1; z <= 118; ++z) {
        std::string a = symbol;
        std::string b = sbox::elements::get_element(z).symbol;
        std::transform(a.begin(), a.end(), a.begin(), [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
        std::transform(b.begin(), b.end(), b.begin(), [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
        if (a == b) {
            return z;
        }
    }
    throw std::runtime_error("Unknown element symbol in trajectory XYZ: " + symbol);
}

double parse_tag_value(const std::string& comment, const std::regex& pattern) {
    std::smatch match;
    if (std::regex_search(comment, match, pattern) && match.size() >= 2) {
        return std::stod(match[1].str());
    }
    return 0.0;
}

TrajectoryFrame read_frame(std::istream& input, int frame_index) {
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty()) {
            break;
        }
    }
    if (!input || line.empty()) {
        throw std::runtime_error("No more trajectory frames");
    }

    std::istringstream count_stream(line);
    int atom_count = 0;
    if (!(count_stream >> atom_count) || atom_count < 0) {
        throw std::runtime_error("Malformed trajectory XYZ: invalid atom count");
    }

    std::string comment;
    if (!std::getline(input, comment)) {
        throw std::runtime_error("Malformed trajectory XYZ: missing comment line");
    }

    TrajectoryFrame frame;
    frame.comment = comment;
    frame.frame_index = frame_index;
    frame.energy = parse_tag_value(comment, std::regex(R"((?:\bE\b|\benergy\b)\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?))", std::regex::icase));
    frame.time_fs = parse_tag_value(comment, std::regex(R"((?:\bt\b|\btime\b)\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?))", std::regex::icase));
    frame.geometry.set_name(comment);

    for (int i = 0; i < atom_count; ++i) {
        if (!std::getline(input, line)) {
            throw std::runtime_error("Malformed trajectory XYZ: truncated atom block");
        }
        std::istringstream atom_stream(line);
        std::string symbol;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        if (!(atom_stream >> symbol >> x >> y >> z)) {
            throw std::runtime_error("Malformed trajectory XYZ: invalid atom line");
        }
        frame.geometry.add_atom({symbol_to_z(symbol), Eigen::Vector3d(x * kAngstromToBohr, y * kAngstromToBohr, z * kAngstromToBohr), "", 0});
    }
    frame.geometry.perceive_bonds();
    return frame;
}

}  // namespace

double Trajectory::min_energy() const {
    if (frames.empty()) {
        return 0.0;
    }
    double min_value = std::numeric_limits<double>::infinity();
    for (const auto& frame : frames) {
        min_value = std::min(min_value, frame.energy);
    }
    return std::isfinite(min_value) ? min_value : 0.0;
}

double Trajectory::max_energy() const {
    if (frames.empty()) {
        return 0.0;
    }
    double max_value = -std::numeric_limits<double>::infinity();
    for (const auto& frame : frames) {
        max_value = std::max(max_value, frame.energy);
    }
    return std::isfinite(max_value) ? max_value : 0.0;
}

std::vector<double> Trajectory::energies() const {
    std::vector<double> out;
    out.reserve(frames.size());
    for (const auto& frame : frames) {
        out.push_back(frame.energy);
    }
    return out;
}

sbox::chem::MolecularSystem Trajectory::interpolate(double t) const {
    if (frames.empty()) {
        throw std::runtime_error("Cannot interpolate an empty trajectory");
    }
    if (frames.size() == 1) {
        return frames.front().geometry;
    }

    const double clamped = std::clamp(t, 0.0, static_cast<double>(frames.size() - 1));
    const int frame_a = static_cast<int>(std::floor(clamped));
    const int frame_b = static_cast<int>(std::ceil(clamped));
    const double frac = clamped - static_cast<double>(frame_a);

    if (frame_a == frame_b) {
        return frames[static_cast<std::size_t>(frame_a)].geometry;
    }

    const auto& a = frames[static_cast<std::size_t>(frame_a)].geometry;
    const auto& b = frames[static_cast<std::size_t>(frame_b)].geometry;
    if (a.num_atoms() != b.num_atoms()) {
        throw std::runtime_error("Cannot interpolate trajectory frames with different atom counts");
    }

    sbox::chem::MolecularSystem mol = a;
    for (int i = 0; i < mol.num_atoms(); ++i) {
        mol.atom(i).position = a.atom(i).position * (1.0 - frac) + b.atom(i).position * frac;
    }
    mol.perceive_bonds();
    return mol;
}

Trajectory read_trajectory_xyz(const std::string& filepath) {
    std::ifstream input(filepath);
    if (!input) {
        throw std::runtime_error("Could not open trajectory XYZ file: " + filepath);
    }

    Trajectory traj;
    int frame_index = 0;
    while (true) {
        const auto pos = input.tellg();
        std::string probe;
        while (std::getline(input, probe)) {
            if (!probe.empty()) {
                break;
            }
        }
        if (!input) {
            break;
        }
        input.clear();
        input.seekg(pos);
        traj.frames.push_back(read_frame(input, frame_index++));
    }
    return traj;
}

void write_trajectory_xyz(const std::string& filepath, const Trajectory& traj) {
    std::ofstream output(filepath);
    if (!output) {
        throw std::runtime_error("Could not open trajectory XYZ file for writing: " + filepath);
    }

    output << std::fixed << std::setprecision(10);
    for (const auto& frame : traj.frames) {
        output << frame.geometry.num_atoms() << '\n';
        if (!frame.comment.empty()) {
            output << frame.comment << '\n';
        } else {
            output << "Frame " << frame.frame_index << ", E = " << frame.energy << " Hartree, t = " << frame.time_fs << " fs\n";
        }
        for (const auto& atom : frame.geometry.atoms()) {
            output << sbox::elements::get_element(atom.Z).symbol << ' '
                   << atom.position.x() * kBohrToAngstrom << ' '
                   << atom.position.y() * kBohrToAngstrom << ' '
                   << atom.position.z() * kBohrToAngstrom << '\n';
        }
    }
}

}  // namespace sbox::io
