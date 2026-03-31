#include "analysis/nci.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace sbox::analysis {

namespace {

constexpr double kDensityEpsilon = 1.0e-10;

std::size_t index_3d(int ix, int iy, int iz, int ny, int nz) {
    return static_cast<std::size_t>(ix) * static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz) +
           static_cast<std::size_t>(iy) * static_cast<std::size_t>(nz) +
           static_cast<std::size_t>(iz);
}

int clamp_index(int value, int max_value) {
    return std::clamp(value, 0, max_value - 1);
}

float sample(const sbox::io::CubeData& cube, int ix, int iy, int iz) {
    ix = clamp_index(ix, cube.nx);
    iy = clamp_index(iy, cube.ny);
    iz = clamp_index(iz, cube.nz);
    return cube.data[index_3d(ix, iy, iz, cube.ny, cube.nz)];
}

double axis_step(const Eigen::Vector3d& v) {
    return v.norm();
}

}  // namespace

NCIGrid compute_nci(const sbox::io::CubeData& density_cube, float rdg_cutoff, float rho_cutoff) {
    if (density_cube.nx <= 2 || density_cube.ny <= 2 || density_cube.nz <= 2 || density_cube.data.empty()) {
        throw std::runtime_error("Density cube is too small for NCI analysis");
    }

    const double hx = axis_step(density_cube.step_x);
    const double hy = axis_step(density_cube.step_y);
    const double hz = axis_step(density_cube.step_z);
    if (hx <= 0.0 || hy <= 0.0 || hz <= 0.0) {
        throw std::runtime_error("Density cube has invalid grid spacing");
    }

    NCIGrid grid;
    grid.origin = density_cube.origin;
    grid.step = Eigen::Vector3d(hx, hy, hz);
    grid.nx = density_cube.nx;
    grid.ny = density_cube.ny;
    grid.nz = density_cube.nz;
    grid.rdg_cutoff = rdg_cutoff;
    grid.rho_cutoff = rho_cutoff;
    grid.rdg.resize(density_cube.data.size(), 10.0f);
    grid.sign_lambda2_rho.resize(density_cube.data.size(), 0.0f);

    const bool use_full_hessian = static_cast<long long>(grid.nx) * grid.ny * grid.nz < 100LL * 100LL * 100LL;
    const double rdg_denom_factor = 2.0 * std::cbrt(3.0 * M_PI * M_PI);

    for (int ix = 0; ix < grid.nx; ++ix) {
        for (int iy = 0; iy < grid.ny; ++iy) {
            for (int iz = 0; iz < grid.nz; ++iz) {
                const std::size_t idx = index_3d(ix, iy, iz, grid.ny, grid.nz);
                const double rho = std::max<double>(density_cube.data[idx], 0.0);
                if (rho < kDensityEpsilon) {
                    grid.rdg[idx] = 10.0f;
                    grid.sign_lambda2_rho[idx] = 0.0f;
                    continue;
                }

                const double drdx = (sample(density_cube, ix + 1, iy, iz) - sample(density_cube, ix - 1, iy, iz)) / (2.0 * hx);
                const double drdy = (sample(density_cube, ix, iy + 1, iz) - sample(density_cube, ix, iy - 1, iz)) / (2.0 * hy);
                const double drdz = (sample(density_cube, ix, iy, iz + 1) - sample(density_cube, ix, iy, iz - 1)) / (2.0 * hz);
                const double grad_norm = std::sqrt(drdx * drdx + drdy * drdy + drdz * drdz);
                const double rdg = grad_norm / (rdg_denom_factor * std::pow(rho, 4.0 / 3.0));

                double sign_term = 0.0;
                const double hxx = (sample(density_cube, ix + 1, iy, iz) - 2.0 * rho + sample(density_cube, ix - 1, iy, iz)) / (hx * hx);
                const double hyy = (sample(density_cube, ix, iy + 1, iz) - 2.0 * rho + sample(density_cube, ix, iy - 1, iz)) / (hy * hy);
                const double hzz = (sample(density_cube, ix, iy, iz + 1) - 2.0 * rho + sample(density_cube, ix, iy, iz - 1)) / (hz * hz);

                if (use_full_hessian) {
                    const double hxy = (sample(density_cube, ix + 1, iy + 1, iz) -
                                        sample(density_cube, ix + 1, iy - 1, iz) -
                                        sample(density_cube, ix - 1, iy + 1, iz) +
                                        sample(density_cube, ix - 1, iy - 1, iz)) / (4.0 * hx * hy);
                    const double hxz = (sample(density_cube, ix + 1, iy, iz + 1) -
                                        sample(density_cube, ix + 1, iy, iz - 1) -
                                        sample(density_cube, ix - 1, iy, iz + 1) +
                                        sample(density_cube, ix - 1, iy, iz - 1)) / (4.0 * hx * hz);
                    const double hyz = (sample(density_cube, ix, iy + 1, iz + 1) -
                                        sample(density_cube, ix, iy + 1, iz - 1) -
                                        sample(density_cube, ix, iy - 1, iz + 1) +
                                        sample(density_cube, ix, iy - 1, iz - 1)) / (4.0 * hy * hz);
                    Eigen::Matrix3d hessian;
                    hessian << hxx, hxy, hxz,
                               hxy, hyy, hyz,
                               hxz, hyz, hzz;
                    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(hessian, Eigen::EigenvaluesOnly);
                    const Eigen::Vector3d eigenvalues = solver.eigenvalues();
                    const double lambda2 = eigenvalues[1];
                    sign_term = (lambda2 >= 0.0 ? 1.0 : -1.0) * rho;
                } else {
                    const double laplacian = hxx + hyy + hzz;
                    sign_term = (laplacian >= 0.0 ? 1.0 : -1.0) * rho;
                }

                if (rho > rho_cutoff || rdg > rdg_cutoff) {
                    grid.rdg[idx] = 10.0f;
                    grid.sign_lambda2_rho[idx] = 0.0f;
                } else {
                    grid.rdg[idx] = static_cast<float>(rdg);
                    grid.sign_lambda2_rho[idx] = static_cast<float>(sign_term);
                }
            }
        }
    }

    return grid;
}

}  // namespace sbox::analysis
