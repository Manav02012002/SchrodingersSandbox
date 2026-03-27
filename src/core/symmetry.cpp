#include "core/symmetry.h"

#include "core/elements.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace sbox::chem {
namespace {

constexpr double kAxisDuplicateTol = 1e-3;
constexpr double kNormTol = 1e-8;
constexpr double kRightAngleTol = 2e-2;
constexpr double kMomentTol = 1e-3;

struct CenteredMolecule {
    std::vector<int> Z;
    std::vector<Eigen::Vector3d> pos;
};

Eigen::Vector3d canonicalize_axis(Eigen::Vector3d axis) {
    const double norm = axis.norm();
    if (norm <= kNormTol) {
        return Eigen::Vector3d::Zero();
    }
    axis /= norm;
    for (int i = 0; i < 3; ++i) {
        if (std::abs(axis[i]) > kNormTol) {
            if (axis[i] < 0.0) {
                axis = -axis;
            }
            break;
        }
    }
    return axis;
}

bool axes_equivalent(const Eigen::Vector3d& a, const Eigen::Vector3d& b, double tol = kAxisDuplicateTol) {
    return (a - b).norm() < tol || (a + b).norm() < tol;
}

void add_unique_axis(std::vector<Eigen::Vector3d>& axes, const Eigen::Vector3d& axis) {
    const Eigen::Vector3d canon = canonicalize_axis(axis);
    if (canon.norm() <= kNormTol) {
        return;
    }
    for (const Eigen::Vector3d& existing : axes) {
        if (axes_equivalent(existing, canon)) {
            return;
        }
    }
    axes.push_back(canon);
}

CenteredMolecule center_molecule(const MolecularSystem& mol) {
    CenteredMolecule centered;
    centered.Z.reserve(static_cast<std::size_t>(mol.num_atoms()));
    centered.pos.reserve(static_cast<std::size_t>(mol.num_atoms()));

    const Eigen::Vector3d com = mol.center_of_mass();
    for (int i = 0; i < mol.num_atoms(); ++i) {
        centered.Z.push_back(mol.atom(i).Z);
        centered.pos.push_back(mol.atom(i).position - com);
    }
    return centered;
}

bool is_symmetry_operation(const CenteredMolecule& mol, const Eigen::Matrix3d& op, double tol) {
    std::vector<bool> used(mol.pos.size(), false);

    for (std::size_t i = 0; i < mol.pos.size(); ++i) {
        const Eigen::Vector3d transformed = op * mol.pos[i];
        int best_match = -1;
        double best_dist = std::numeric_limits<double>::max();

        for (std::size_t j = 0; j < mol.pos.size(); ++j) {
            if (used[j] || mol.Z[j] != mol.Z[i]) {
                continue;
            }
            const double dist = (transformed - mol.pos[j]).norm();
            if (dist < best_dist) {
                best_dist = dist;
                best_match = static_cast<int>(j);
            }
        }

        if (best_match < 0 || best_dist > tol) {
            return false;
        }
        used[static_cast<std::size_t>(best_match)] = true;
    }

    return true;
}

Eigen::Matrix3d rotation_matrix(const Eigen::Vector3d& axis, double angle) {
    const Eigen::Vector3d u = axis.normalized();
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const Eigen::Matrix3d ux = (Eigen::Matrix3d() <<
        0.0, -u.z(), u.y(),
        u.z(), 0.0, -u.x(),
        -u.y(), u.x(), 0.0).finished();
    return c * Eigen::Matrix3d::Identity() + (1.0 - c) * (u * u.transpose()) + s * ux;
}

Eigen::Matrix3d reflection_matrix(const Eigen::Vector3d& normal) {
    const Eigen::Vector3d n = normal.normalized();
    return Eigen::Matrix3d::Identity() - 2.0 * (n * n.transpose());
}

bool check_cn(const CenteredMolecule& mol, const Eigen::Vector3d& axis, int n, double tol) {
    return is_symmetry_operation(mol, rotation_matrix(axis, 2.0 * M_PI / static_cast<double>(n)), tol);
}

bool check_sn(const CenteredMolecule& mol, const Eigen::Vector3d& axis, int n, double tol) {
    const Eigen::Matrix3d op = reflection_matrix(axis) * rotation_matrix(axis, 2.0 * M_PI / static_cast<double>(n));
    return is_symmetry_operation(mol, op, tol);
}

bool all_collinear(const CenteredMolecule& mol, double tol) {
    Eigen::Vector3d axis = Eigen::Vector3d::Zero();
    for (const Eigen::Vector3d& pos : mol.pos) {
        if (pos.norm() > tol) {
            axis = pos.normalized();
            break;
        }
    }
    if (axis.norm() <= kNormTol) {
        return true;
    }

    for (const Eigen::Vector3d& pos : mol.pos) {
        if (pos.cross(axis).norm() > tol) {
            return false;
        }
    }
    return true;
}

Eigen::Matrix3d inertia_tensor(const CenteredMolecule& mol) {
    Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero();
    for (std::size_t i = 0; i < mol.pos.size(); ++i) {
        const double mass = sbox::elements::get_element(mol.Z[i]).atomic_mass;
        const Eigen::Vector3d& r = mol.pos[i];
        tensor += mass * ((r.squaredNorm() * Eigen::Matrix3d::Identity()) - (r * r.transpose()));
    }
    return tensor;
}

PointGroup cn_group(int n, bool sigma_v) {
    switch (n) {
        case 2: return sigma_v ? PointGroup::C2v : PointGroup::C2;
        case 3: return sigma_v ? PointGroup::C3v : PointGroup::C3;
        case 4: return sigma_v ? PointGroup::C4v : PointGroup::C4;
        case 5: return sigma_v ? PointGroup::C5v : PointGroup::C5;
        case 6: return sigma_v ? PointGroup::C6v : PointGroup::C6;
        default: return PointGroup::Unknown;
    }
}

PointGroup cnh_group(int n) {
    switch (n) {
        case 2: return PointGroup::C2h;
        case 3: return PointGroup::C3h;
        case 4: return PointGroup::C4h;
        case 5: return PointGroup::C5h;
        case 6: return PointGroup::C6h;
        default: return PointGroup::Unknown;
    }
}

PointGroup dn_group(int n) {
    switch (n) {
        case 2: return PointGroup::D2;
        case 3: return PointGroup::D3;
        case 4: return PointGroup::D4;
        case 5: return PointGroup::D5;
        case 6: return PointGroup::D6;
        default: return PointGroup::Unknown;
    }
}

PointGroup dnh_group(int n) {
    switch (n) {
        case 2: return PointGroup::D2h;
        case 3: return PointGroup::D3h;
        case 4: return PointGroup::D4h;
        case 5: return PointGroup::D5h;
        case 6: return PointGroup::D6h;
        default: return PointGroup::Unknown;
    }
}

PointGroup dnd_group(int n) {
    switch (n) {
        case 2: return PointGroup::D2d;
        case 3: return PointGroup::D3d;
        case 4: return PointGroup::D4d;
        case 5: return PointGroup::D5d;
        case 6: return PointGroup::D6d;
        default: return PointGroup::Unknown;
    }
}

std::vector<Eigen::Vector3d> build_candidate_axes(const CenteredMolecule& mol, const Eigen::Matrix3d& eigenvectors) {
    std::vector<Eigen::Vector3d> axes;

    for (int i = 0; i < 3; ++i) {
        add_unique_axis(axes, eigenvectors.col(i));
    }

    for (const Eigen::Vector3d& pos : mol.pos) {
        add_unique_axis(axes, pos);
    }

    for (std::size_t i = 0; i < mol.pos.size(); ++i) {
        for (std::size_t j = i + 1; j < mol.pos.size(); ++j) {
            add_unique_axis(axes, mol.pos[i] + mol.pos[j]);
            add_unique_axis(axes, mol.pos[i] - mol.pos[j]);
            add_unique_axis(axes, mol.pos[i].cross(mol.pos[j]));
        }
    }

    return axes;
}

std::vector<Eigen::Vector3d> build_vertical_plane_normals(const CenteredMolecule& mol,
                                                          const Eigen::Vector3d& principal_axis,
                                                          const std::vector<Eigen::Vector3d>& c2_perp_axes) {
    std::vector<Eigen::Vector3d> normals;

    for (const Eigen::Vector3d& axis : c2_perp_axes) {
        add_unique_axis(normals, principal_axis.cross(axis));
    }

    std::vector<Eigen::Vector3d> projected;
    for (const Eigen::Vector3d& pos : mol.pos) {
        Eigen::Vector3d proj = pos - pos.dot(principal_axis) * principal_axis;
        if (proj.norm() > kNormTol) {
            proj.normalize();
            projected.push_back(proj);
            add_unique_axis(normals, principal_axis.cross(proj));
        }
    }

    for (std::size_t i = 0; i < projected.size(); ++i) {
        for (std::size_t j = i + 1; j < projected.size(); ++j) {
            add_unique_axis(normals, principal_axis.cross(projected[i] + projected[j]));
            add_unique_axis(normals, principal_axis.cross(projected[i] - projected[j]));
        }
    }

    return normals;
}

bool has_vertical_plane(const CenteredMolecule& mol,
                        const Eigen::Vector3d& principal_axis,
                        const std::vector<Eigen::Vector3d>& normals) {
    for (const Eigen::Vector3d& normal : normals) {
        if (std::abs(normal.dot(principal_axis)) > 1e-2) {
            continue;
        }
        if (is_symmetry_operation(mol, reflection_matrix(normal), 0.3)) {
            return true;
        }
    }
    return false;
}

bool has_sigma_d_plane(const CenteredMolecule& mol,
                       const Eigen::Vector3d& principal_axis,
                       const std::vector<Eigen::Vector3d>& c2_perp_axes,
                       const std::vector<Eigen::Vector3d>& normals) {
    if (c2_perp_axes.empty()) {
        return false;
    }

    for (const Eigen::Vector3d& normal : normals) {
        if (std::abs(normal.dot(principal_axis)) > 1e-2) {
            continue;
        }
        if (!is_symmetry_operation(mol, reflection_matrix(normal), 0.3)) {
            continue;
        }

        const Eigen::Vector3d plane_dir = canonicalize_axis(principal_axis.cross(normal));
        bool contains_c2 = false;
        for (const Eigen::Vector3d& c2_axis : c2_perp_axes) {
            if (std::abs(std::abs(plane_dir.dot(c2_axis)) - 1.0) < 5e-2) {
                contains_c2 = true;
                break;
            }
        }
        if (!contains_c2) {
            return true;
        }
    }

    return false;
}

double effective_tol(const CenteredMolecule& mol, double tol) {
    double scale = 0.0;
    for (const Eigen::Vector3d& pos : mol.pos) {
        scale = std::max(scale, pos.norm());
    }
    return std::max(tol, 0.02 * std::max(1.0, scale));
}

}  // namespace

std::string point_group_name(PointGroup pg) {
    switch (pg) {
        case PointGroup::C1: return "C1";
        case PointGroup::Ci: return "Ci";
        case PointGroup::Cs: return "Cs";
        case PointGroup::C2: return "C2";
        case PointGroup::C3: return "C3";
        case PointGroup::C4: return "C4";
        case PointGroup::C5: return "C5";
        case PointGroup::C6: return "C6";
        case PointGroup::C2v: return "C2v";
        case PointGroup::C3v: return "C3v";
        case PointGroup::C4v: return "C4v";
        case PointGroup::C5v: return "C5v";
        case PointGroup::C6v: return "C6v";
        case PointGroup::C2h: return "C2h";
        case PointGroup::C3h: return "C3h";
        case PointGroup::C4h: return "C4h";
        case PointGroup::C5h: return "C5h";
        case PointGroup::C6h: return "C6h";
        case PointGroup::D2: return "D2";
        case PointGroup::D3: return "D3";
        case PointGroup::D4: return "D4";
        case PointGroup::D5: return "D5";
        case PointGroup::D6: return "D6";
        case PointGroup::D2h: return "D2h";
        case PointGroup::D3h: return "D3h";
        case PointGroup::D4h: return "D4h";
        case PointGroup::D5h: return "D5h";
        case PointGroup::D6h: return "D6h";
        case PointGroup::D2d: return "D2d";
        case PointGroup::D3d: return "D3d";
        case PointGroup::D4d: return "D4d";
        case PointGroup::D5d: return "D5d";
        case PointGroup::D6d: return "D6d";
        case PointGroup::S4: return "S4";
        case PointGroup::S6: return "S6";
        case PointGroup::S8: return "S8";
        case PointGroup::T: return "T";
        case PointGroup::Td: return "Td";
        case PointGroup::Th: return "Th";
        case PointGroup::O: return "O";
        case PointGroup::Oh: return "Oh";
        case PointGroup::Ih: return "Ih";
        case PointGroup::Cinfv: return "Cinfv";
        case PointGroup::Dinfh: return "Dinfh";
        case PointGroup::Unknown: return "Unknown";
    }
    return "Unknown";
}

PointGroup detect_point_group(const MolecularSystem& mol, double tolerance) {
    if (mol.num_atoms() <= 1) {
        return PointGroup::C1;
    }

    if (mol.num_atoms() == 2) {
        return mol.atom(0).Z == mol.atom(1).Z ? PointGroup::Dinfh : PointGroup::Cinfv;
    }

    const CenteredMolecule centered = center_molecule(mol);
    const double tol = effective_tol(centered, tolerance);

    if (all_collinear(centered, tol)) {
        const bool inversion = is_symmetry_operation(centered, -Eigen::Matrix3d::Identity(), tol);
        return inversion ? PointGroup::Dinfh : PointGroup::Cinfv;
    }

    const Eigen::Matrix3d inertia = inertia_tensor(centered);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(inertia);
    const Eigen::Vector3d moments = solver.eigenvalues();
    const Eigen::Matrix3d eigenvectors = solver.eigenvectors();

    const bool spherical_top =
        std::abs(moments[0] - moments[1]) < kMomentTol && std::abs(moments[1] - moments[2]) < kMomentTol;

    const bool inversion = is_symmetry_operation(centered, -Eigen::Matrix3d::Identity(), tol);
    std::vector<Eigen::Vector3d> candidate_axes = build_candidate_axes(centered, eigenvectors);

    int c3_count = 0;
    int c4_count = 0;
    int c5_count = 0;
    if (spherical_top) {
        for (const Eigen::Vector3d& axis : candidate_axes) {
            if (check_cn(centered, axis, 3, tol)) {
                ++c3_count;
            }
            if (check_cn(centered, axis, 4, tol)) {
                ++c4_count;
            }
            if (check_cn(centered, axis, 5, tol)) {
                ++c5_count;
            }
        }

        if (c5_count >= 6) {
            return PointGroup::Ih;
        }
        if (c3_count >= 4 && c4_count >= 3) {
            return inversion ? PointGroup::Oh : PointGroup::O;
        }
        if (c3_count >= 4) {
            return inversion ? PointGroup::Th : PointGroup::Td;
        }
    }

    int principal_n = 1;
    Eigen::Vector3d principal_axis = eigenvectors.col(2);
    for (int n = 6; n >= 2; --n) {
        for (const Eigen::Vector3d& axis : candidate_axes) {
            if (check_cn(centered, axis, n, tol)) {
                principal_n = n;
                principal_axis = axis;
                break;
            }
        }
        if (principal_n == n) {
            break;
        }
    }

    if (principal_n == 1) {
        std::vector<Eigen::Vector3d> plane_normals = build_candidate_axes(centered, eigenvectors);
        for (const Eigen::Vector3d& normal : plane_normals) {
            if (is_symmetry_operation(centered, reflection_matrix(normal), tol)) {
                return PointGroup::Cs;
            }
        }
        return inversion ? PointGroup::Ci : PointGroup::C1;
    }

    principal_axis = canonicalize_axis(principal_axis);

    std::vector<Eigen::Vector3d> c2_perp_axes;
    for (const Eigen::Vector3d& axis : candidate_axes) {
        if (std::abs(std::abs(axis.dot(principal_axis)) - 1.0) < 5e-2) {
            continue;
        }
        if (std::abs(axis.dot(principal_axis)) > kRightAngleTol) {
            continue;
        }
        if (check_cn(centered, axis, 2, tol)) {
            add_unique_axis(c2_perp_axes, axis);
        }
    }

    const bool has_n_c2_perp = static_cast<int>(c2_perp_axes.size()) >= principal_n;
    const bool sigma_h = is_symmetry_operation(centered, reflection_matrix(principal_axis), tol);
    const std::vector<Eigen::Vector3d> vertical_normals =
        build_vertical_plane_normals(centered, principal_axis, c2_perp_axes);
    const bool sigma_v = has_vertical_plane(centered, principal_axis, vertical_normals);
    const bool sigma_d = has_sigma_d_plane(centered, principal_axis, c2_perp_axes, vertical_normals);

    if (has_n_c2_perp) {
        if (sigma_h) {
            return dnh_group(principal_n);
        }
        if (sigma_d || inversion) {
            return dnd_group(principal_n);
        }
        return dn_group(principal_n);
    }

    if (sigma_h || inversion) {
        return cnh_group(principal_n);
    }
    if (sigma_v) {
        return cn_group(principal_n, true);
    }
    if (principal_n == 2 && check_sn(centered, principal_axis, 4, tol)) {
        return PointGroup::S4;
    }
    if (principal_n == 3 && check_sn(centered, principal_axis, 6, tol)) {
        return PointGroup::S6;
    }
    return cn_group(principal_n, false);
}

}  // namespace sbox::chem
