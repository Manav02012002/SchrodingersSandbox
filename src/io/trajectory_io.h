#pragma once

#include "core/molecular_system.h"

#include <string>
#include <vector>

namespace sbox::io {

struct TrajectoryFrame {
    sbox::chem::MolecularSystem geometry;
    double energy = 0.0;
    double time_fs = 0.0;
    std::string comment;
    int frame_index = 0;
};

struct Trajectory {
    std::vector<TrajectoryFrame> frames;

    int num_frames() const { return static_cast<int>(frames.size()); }
    bool empty() const { return frames.empty(); }

    double min_energy() const;
    double max_energy() const;
    std::vector<double> energies() const;
    sbox::chem::MolecularSystem interpolate(double t) const;
};

Trajectory read_trajectory_xyz(const std::string& filepath);
void write_trajectory_xyz(const std::string& filepath, const Trajectory& traj);

}  // namespace sbox::io
