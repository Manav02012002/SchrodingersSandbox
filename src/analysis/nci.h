#pragma once

#include "io/cube_io.h"

#include <Eigen/Core>

#include <vector>

namespace sbox::analysis {

struct NCIGrid {
    Eigen::Vector3d origin = Eigen::Vector3d::Zero();
    Eigen::Vector3d step = Eigen::Vector3d::Ones();
    int nx = 0;
    int ny = 0;
    int nz = 0;
    std::vector<float> rdg;
    std::vector<float> sign_lambda2_rho;
    float rdg_cutoff = 0.5f;
    float rho_cutoff = 0.05f;
};

NCIGrid compute_nci(
    const sbox::io::CubeData& density_cube,
    float rdg_cutoff = 0.5f,
    float rho_cutoff = 0.05f);

}  // namespace sbox::analysis
