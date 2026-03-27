#include "core/gaussian_eval.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace sbox::basis {

namespace {

constexpr double kSqrt3 = 1.73205080756887729353;
constexpr double kSqrt15 = 3.87298334620741688518;
constexpr double kSqrt3Over8 = 0.61237243569579452455;
constexpr double kSqrt5Over8 = 0.79056941504209483299;

double pow_int(double value, int exponent) {
    if (exponent <= 0) {
        return 1.0;
    }

    double result = 1.0;
    for (int i = 0; i < exponent; ++i) {
        result *= value;
    }
    return result;
}

int shell_basis_count(int angular_momentum, bool spherical) {
    switch (angular_momentum) {
    case 0:
        return 1;
    case 1:
        return 3;
    case 2:
        return spherical ? 5 : 6;
    case 3:
        return spherical ? 7 : 10;
    default:
        throw std::runtime_error("Unsupported shell angular momentum");
    }
}

double radial_sum_for_shell(const BasisShell& shell, double r2) {
    double radial = 0.0;
    for (const GaussianPrimitive& primitive : shell.primitives) {
        radial += primitive.coefficient * std::exp(-primitive.exponent * r2);
    }
    return radial;
}

int evaluate_shell_impl(const BasisShell& shell,
                        const Eigen::Vector3d& center,
                        const Eigen::Vector3d& point,
                        int basis_offset,
                        std::vector<double>& out_values,
                        bool spherical) {
    const Eigen::Vector3d dr = point - center;
    const double dx = dr.x();
    const double dy = dr.y();
    const double dz = dr.z();
    const double r2 = dr.squaredNorm();
    const double radial = radial_sum_for_shell(shell, r2);

    const int count = shell_basis_count(shell.angular_momentum, spherical);
    if (basis_offset < 0 || basis_offset + count > static_cast<int>(out_values.size())) {
        throw std::runtime_error("Basis output buffer too small for shell");
    }

    switch (shell.angular_momentum) {
    case 0:
        out_values[static_cast<std::size_t>(basis_offset)] = radial;
        return 1;
    case 1:
        out_values[static_cast<std::size_t>(basis_offset + 0)] = dx * radial;
        out_values[static_cast<std::size_t>(basis_offset + 1)] = dy * radial;
        out_values[static_cast<std::size_t>(basis_offset + 2)] = dz * radial;
        return 3;
    case 2: {
        const double xx = dx * dx * radial;
        const double yy = dy * dy * radial;
        const double zz = dz * dz * radial;
        const double xy = dx * dy * radial;
        const double xz = dx * dz * radial;
        const double yz = dy * dz * radial;

        if (!spherical) {
            out_values[static_cast<std::size_t>(basis_offset + 0)] = xx;
            out_values[static_cast<std::size_t>(basis_offset + 1)] = yy;
            out_values[static_cast<std::size_t>(basis_offset + 2)] = zz;
            out_values[static_cast<std::size_t>(basis_offset + 3)] = xy;
            out_values[static_cast<std::size_t>(basis_offset + 4)] = xz;
            out_values[static_cast<std::size_t>(basis_offset + 5)] = yz;
            return 6;
        }

        out_values[static_cast<std::size_t>(basis_offset + 0)] = 0.5 * (2.0 * zz - xx - yy);
        out_values[static_cast<std::size_t>(basis_offset + 1)] = kSqrt3 * xz;
        out_values[static_cast<std::size_t>(basis_offset + 2)] = kSqrt3 * yz;
        out_values[static_cast<std::size_t>(basis_offset + 3)] = 0.5 * kSqrt3 * (xx - yy);
        out_values[static_cast<std::size_t>(basis_offset + 4)] = kSqrt3 * xy;
        return 5;
    }
    case 3: {
        const double xx = dx * dx;
        const double yy = dy * dy;
        const double zz = dz * dz;

        const double xxx = dx * xx * radial;
        const double yyy = dy * yy * radial;
        const double zzz = dz * zz * radial;
        const double xyy = dx * yy * radial;
        const double xxy = xx * dy * radial;
        const double xxz = xx * dz * radial;
        const double xzz = dx * zz * radial;
        const double yzz = dy * zz * radial;
        const double yyz = yy * dz * radial;
        const double xyz = dx * dy * dz * radial;

        if (!spherical) {
            out_values[static_cast<std::size_t>(basis_offset + 0)] = xxx;
            out_values[static_cast<std::size_t>(basis_offset + 1)] = yyy;
            out_values[static_cast<std::size_t>(basis_offset + 2)] = zzz;
            out_values[static_cast<std::size_t>(basis_offset + 3)] = xyy;
            out_values[static_cast<std::size_t>(basis_offset + 4)] = xxy;
            out_values[static_cast<std::size_t>(basis_offset + 5)] = xxz;
            out_values[static_cast<std::size_t>(basis_offset + 6)] = xzz;
            out_values[static_cast<std::size_t>(basis_offset + 7)] = yzz;
            out_values[static_cast<std::size_t>(basis_offset + 8)] = yyz;
            out_values[static_cast<std::size_t>(basis_offset + 9)] = xyz;
            return 10;
        }

        out_values[static_cast<std::size_t>(basis_offset + 0)] =
            0.5 * dz * (2.0 * zz - 3.0 * xx - 3.0 * yy) * radial;
        out_values[static_cast<std::size_t>(basis_offset + 1)] =
            kSqrt3Over8 * dx * (4.0 * zz - xx - yy) * radial;
        out_values[static_cast<std::size_t>(basis_offset + 2)] =
            kSqrt3Over8 * dy * (4.0 * zz - xx - yy) * radial;
        out_values[static_cast<std::size_t>(basis_offset + 3)] =
            0.5 * kSqrt15 * dz * (xx - yy) * radial;
        out_values[static_cast<std::size_t>(basis_offset + 4)] =
            kSqrt15 * dx * dy * dz * radial;
        out_values[static_cast<std::size_t>(basis_offset + 5)] =
            kSqrt5Over8 * dx * (xx - 3.0 * yy) * radial;
        out_values[static_cast<std::size_t>(basis_offset + 6)] =
            kSqrt5Over8 * dy * (3.0 * xx - yy) * radial;
        return 7;
    }
    default:
        throw std::runtime_error("Unsupported shell angular momentum");
    }
}

}  // namespace

double evaluate_primitive(double alpha,
                          double coeff,
                          int lx,
                          int ly,
                          int lz,
                          double dx,
                          double dy,
                          double dz) {
    const double r2 = dx * dx + dy * dy + dz * dz;
    return coeff * pow_int(dx, lx) * pow_int(dy, ly) * pow_int(dz, lz) * std::exp(-alpha * r2);
}

double evaluate_shell(const BasisShell& shell,
                      const Eigen::Vector3d& center,
                      const Eigen::Vector3d& point,
                      int basis_offset,
                      std::vector<double>& out_values) {
    bool spherical = false;
    if (shell.angular_momentum >= 2) {
        const int remaining = static_cast<int>(out_values.size()) - basis_offset;
        spherical = remaining == shell_basis_count(shell.angular_momentum, true);
    }

    return static_cast<double>(
        evaluate_shell_impl(shell, center, point, basis_offset, out_values, spherical));
}

void evaluate_basis_at_point(const MOData& mo_data,
                             const Eigen::Vector3d& point,
                             Eigen::VectorXd& basis_values) {
    const int num_basis = mo_data.basis.num_basis_functions();
    basis_values = Eigen::VectorXd::Zero(num_basis);
    std::vector<double> values(static_cast<std::size_t>(num_basis), 0.0);

    int basis_offset = 0;
    for (const BasisShell& shell : mo_data.basis.shells) {
        if (shell.atom_index < 0 || shell.atom_index >= static_cast<int>(mo_data.atom_positions.size())) {
            throw std::runtime_error("Basis shell atom index out of range");
        }

        const int count = evaluate_shell_impl(shell,
                                              mo_data.atom_positions[static_cast<std::size_t>(shell.atom_index)],
                                              point,
                                              basis_offset,
                                              values,
                                              mo_data.basis.spherical);
        basis_offset += count;
    }

    if (basis_offset != num_basis) {
        throw std::runtime_error("Basis evaluation produced a size mismatch");
    }

    for (int i = 0; i < num_basis; ++i) {
        basis_values(i) = values[static_cast<std::size_t>(i)];
    }
}

double evaluate_mo_at_point(const MOData& mo_data, int mo_index, const Eigen::Vector3d& point) {
    if (mo_index < 0 || mo_index >= mo_data.coefficients.cols()) {
        throw std::runtime_error("MO index out of range");
    }

    Eigen::VectorXd basis_values;
    evaluate_basis_at_point(mo_data, point, basis_values);
    return basis_values.dot(mo_data.coefficients.col(mo_index));
}

double evaluate_mo_density_at_point(const MOData& mo_data, int mo_index, const Eigen::Vector3d& point) {
    const double psi = evaluate_mo_at_point(mo_data, mo_index, point);
    return psi * psi;
}

Eigen::VectorXd evaluate_mo_on_grid(const MOData& mo_data,
                                    int mo_index,
                                    const Eigen::Vector3d& origin,
                                    const Eigen::Vector3d& step,
                                    int nx,
                                    int ny,
                                    int nz) {
    if (nx < 0 || ny < 0 || nz < 0) {
        throw std::runtime_error("Grid dimensions must be non-negative");
    }

    Eigen::VectorXd values(nx * ny * nz);
    int index = 0;
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                const Eigen::Vector3d point =
                    origin + step.cwiseProduct(Eigen::Vector3d(ix, iy, iz));
                values(index++) = evaluate_mo_at_point(mo_data, mo_index, point);
            }
        }
    }
    return values;
}

}  // namespace sbox::basis
