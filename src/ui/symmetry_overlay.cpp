#include "ui/symmetry_overlay.h"

#include "core/elements.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace sbox::ui {

namespace {

constexpr double kNormTol = 1.0e-8;
constexpr double kAxisDuplicateTol = 1.0e-3;
constexpr double kScreenEpsilon = 1.0e-5;
constexpr float kPiF = 3.14159265358979323846f;
constexpr double kPi = 3.14159265358979323846;

struct CenteredMolecule {
    std::vector<int> atomic_numbers;
    std::vector<Eigen::Vector3d> positions;
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
};

struct FrameAxes {
    Eigen::Vector3d principal = Eigen::Vector3d::UnitZ();
    Eigen::Vector3d secondary = Eigen::Vector3d::UnitX();
    Eigen::Vector3d tertiary = Eigen::Vector3d::UnitY();
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

CenteredMolecule center_molecule(const sbox::chem::MolecularSystem& mol) {
    CenteredMolecule centered;
    centered.center = mol.center_of_mass();
    centered.atomic_numbers.reserve(static_cast<std::size_t>(mol.num_atoms()));
    centered.positions.reserve(static_cast<std::size_t>(mol.num_atoms()));
    for (int i = 0; i < mol.num_atoms(); ++i) {
        centered.atomic_numbers.push_back(mol.atom(i).Z);
        centered.positions.push_back(mol.atom(i).position - centered.center);
    }
    return centered;
}

Eigen::Matrix3d inertia_tensor(const CenteredMolecule& centered) {
    Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero();
    for (std::size_t i = 0; i < centered.positions.size(); ++i) {
        const double mass = sbox::elements::get_element(centered.atomic_numbers[i]).atomic_mass;
        const Eigen::Vector3d& r = centered.positions[i];
        tensor += mass * ((r.squaredNorm() * Eigen::Matrix3d::Identity()) - (r * r.transpose()));
    }
    return tensor;
}

FrameAxes compute_frame_axes(const CenteredMolecule& centered) {
    FrameAxes frame;
    if (centered.positions.empty()) {
        return frame;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(inertia_tensor(centered));
    if (solver.info() != Eigen::Success) {
        return frame;
    }

    std::array<int, 3> order = {0, 1, 2};
    const Eigen::Vector3d eigenvalues = solver.eigenvalues();
    std::sort(order.begin(), order.end(), [&](int a, int b) { return eigenvalues[a] < eigenvalues[b]; });

    frame.principal = canonicalize_axis(solver.eigenvectors().col(order[0]));
    frame.secondary = canonicalize_axis(solver.eigenvectors().col(order[1]));
    frame.tertiary = canonicalize_axis(solver.eigenvectors().col(order[2]));

    if (frame.principal.norm() <= kNormTol) frame.principal = Eigen::Vector3d::UnitZ();
    if (frame.secondary.norm() <= kNormTol || std::abs(frame.secondary.dot(frame.principal)) > 0.99) {
        frame.secondary = canonicalize_axis(frame.principal.unitOrthogonal());
    }
    frame.tertiary = canonicalize_axis(frame.principal.cross(frame.secondary));
    if (frame.tertiary.norm() <= kNormTol) {
        frame.tertiary = canonicalize_axis(frame.principal.unitOrthogonal());
    }
    frame.secondary = canonicalize_axis(frame.tertiary.cross(frame.principal));
    return frame;
}

Eigen::Matrix3d rotation_matrix(const Eigen::Vector3d& axis, double angle) {
    const Eigen::Vector3d unit = axis.normalized();
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const Eigen::Matrix3d skew = (Eigen::Matrix3d() <<
        0.0, -unit.z(), unit.y(),
        unit.z(), 0.0, -unit.x(),
        -unit.y(), unit.x(), 0.0).finished();
    return c * Eigen::Matrix3d::Identity() + (1.0 - c) * (unit * unit.transpose()) + s * skew;
}

Eigen::Matrix3d reflection_matrix(const Eigen::Vector3d& normal) {
    const Eigen::Vector3d unit = normal.normalized();
    return Eigen::Matrix3d::Identity() - 2.0 * (unit * unit.transpose());
}

bool is_symmetry_operation(const CenteredMolecule& centered, const Eigen::Matrix3d& op, double tolerance) {
    std::vector<bool> used(centered.positions.size(), false);
    for (std::size_t i = 0; i < centered.positions.size(); ++i) {
        const Eigen::Vector3d transformed = op * centered.positions[i];
        int best_index = -1;
        double best_distance = std::numeric_limits<double>::max();
        for (std::size_t j = 0; j < centered.positions.size(); ++j) {
            if (used[j] || centered.atomic_numbers[j] != centered.atomic_numbers[i]) {
                continue;
            }
            const double distance = (transformed - centered.positions[j]).norm();
            if (distance < best_distance) {
                best_distance = distance;
                best_index = static_cast<int>(j);
            }
        }
        if (best_index < 0 || best_distance > tolerance) {
            return false;
        }
        used[static_cast<std::size_t>(best_index)] = true;
    }
    return true;
}

bool check_cn(const CenteredMolecule& centered, const Eigen::Vector3d& axis, int order, double tolerance) {
    return is_symmetry_operation(centered, rotation_matrix(axis, 2.0 * kPi / static_cast<double>(order)), tolerance);
}

bool check_sigma(const CenteredMolecule& centered, const Eigen::Vector3d& normal, double tolerance) {
    return is_symmetry_operation(centered, reflection_matrix(normal), tolerance);
}

bool check_sn(const CenteredMolecule& centered, const Eigen::Vector3d& axis, int order, double tolerance) {
    const Eigen::Matrix3d op = reflection_matrix(axis) * rotation_matrix(axis, 2.0 * kPi / static_cast<double>(order));
    return is_symmetry_operation(centered, op, tolerance);
}

bool check_inversion(const CenteredMolecule& centered, double tolerance) {
    return is_symmetry_operation(centered, -Eigen::Matrix3d::Identity(), tolerance);
}

std::vector<Eigen::Vector3d> build_candidate_axes(const CenteredMolecule& centered, const FrameAxes& frame) {
    std::vector<Eigen::Vector3d> axes;
    add_unique_axis(axes, frame.principal);
    add_unique_axis(axes, frame.secondary);
    add_unique_axis(axes, frame.tertiary);

    for (const Eigen::Vector3d& pos : centered.positions) {
        add_unique_axis(axes, pos);
    }

    for (std::size_t i = 0; i < centered.positions.size(); ++i) {
        for (std::size_t j = i + 1; j < centered.positions.size(); ++j) {
            const Eigen::Vector3d a = centered.positions[i];
            const Eigen::Vector3d b = centered.positions[j];
            add_unique_axis(axes, a + b);
            add_unique_axis(axes, a - b);
            add_unique_axis(axes, a.cross(b));
        }
    }

    for (std::size_t i = 0; i < centered.positions.size(); ++i) {
        for (std::size_t j = i + 1; j < centered.positions.size(); ++j) {
            for (std::size_t k = j + 1; k < centered.positions.size(); ++k) {
                add_unique_axis(axes, centered.positions[i] + centered.positions[j] + centered.positions[k]);
            }
        }
    }

    return axes;
}

std::vector<Eigen::Vector3d> build_candidate_plane_normals(const CenteredMolecule& centered, const FrameAxes& frame) {
    std::vector<Eigen::Vector3d> normals;
    add_unique_axis(normals, frame.principal);
    add_unique_axis(normals, frame.secondary);
    add_unique_axis(normals, frame.tertiary);

    for (const Eigen::Vector3d& pos : centered.positions) {
        add_unique_axis(normals, pos);
        add_unique_axis(normals, frame.principal.cross(pos));
        add_unique_axis(normals, frame.secondary.cross(pos));
        add_unique_axis(normals, frame.tertiary.cross(pos));
    }

    for (std::size_t i = 0; i < centered.positions.size(); ++i) {
        for (std::size_t j = i + 1; j < centered.positions.size(); ++j) {
            const Eigen::Vector3d a = centered.positions[i];
            const Eigen::Vector3d b = centered.positions[j];
            add_unique_axis(normals, a.cross(b));
            add_unique_axis(normals, a + b);
            add_unique_axis(normals, a - b);
        }
    }

    return normals;
}

template <typename Predicate>
std::vector<Eigen::Vector3d> filter_axes(const std::vector<Eigen::Vector3d>& axes, Predicate&& predicate) {
    std::vector<Eigen::Vector3d> out;
    out.reserve(axes.size());
    for (const Eigen::Vector3d& axis : axes) {
        if (predicate(axis)) {
            add_unique_axis(out, axis);
        }
    }
    return out;
}

std::vector<Eigen::Vector3d> find_cn_axes(const CenteredMolecule& centered,
                                          const std::vector<Eigen::Vector3d>& candidates,
                                          int order,
                                          double tolerance) {
    return filter_axes(candidates, [&](const Eigen::Vector3d& axis) {
        return check_cn(centered, axis, order, tolerance);
    });
}

std::vector<Eigen::Vector3d> find_sn_axes(const CenteredMolecule& centered,
                                          const std::vector<Eigen::Vector3d>& candidates,
                                          int order,
                                          double tolerance) {
    return filter_axes(candidates, [&](const Eigen::Vector3d& axis) {
        return check_sn(centered, axis, order, tolerance);
    });
}

std::vector<Eigen::Vector3d> find_sigma_normals(const CenteredMolecule& centered,
                                                const std::vector<Eigen::Vector3d>& candidates,
                                                double tolerance) {
    return filter_axes(candidates, [&](const Eigen::Vector3d& normal) {
        return check_sigma(centered, normal, tolerance);
    });
}

double axis_score(const Eigen::Vector3d& axis, const Eigen::Vector3d& reference) {
    return std::abs(canonicalize_axis(axis).dot(canonicalize_axis(reference)));
}

std::vector<Eigen::Vector3d> sort_by_reference(std::vector<Eigen::Vector3d> axes,
                                               const Eigen::Vector3d& reference,
                                               bool descending = true) {
    std::sort(axes.begin(), axes.end(), [&](const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
        const double sa = axis_score(a, reference);
        const double sb = axis_score(b, reference);
        return descending ? sa > sb : sa < sb;
    });
    return axes;
}

std::vector<Eigen::Vector3d> select_count(const std::vector<Eigen::Vector3d>& axes, std::size_t count) {
    if (axes.size() <= count) {
        return axes;
    }
    return std::vector<Eigen::Vector3d>(axes.begin(), axes.begin() + static_cast<std::ptrdiff_t>(count));
}

void append_axes(std::vector<SymmetryElement>& elements,
                 const std::vector<Eigen::Vector3d>& axes,
                 SymmetryElement::Type type,
                 int order,
                 const std::string& prefix,
                 const Eigen::Vector3d& center) {
    for (const Eigen::Vector3d& axis : axes) {
        SymmetryElement element;
        element.type = type;
        element.axis = canonicalize_axis(axis);
        element.position = center;
        element.order = order;
        element.label = prefix + std::to_string(order);
        elements.push_back(element);
    }
}

void append_planes(std::vector<SymmetryElement>& elements,
                   const std::vector<Eigen::Vector3d>& normals,
                   const std::string& label,
                   const Eigen::Vector3d& center) {
    for (const Eigen::Vector3d& normal : normals) {
        SymmetryElement element;
        element.type = SymmetryElement::Type::Sigma;
        element.axis = canonicalize_axis(normal);
        element.position = center;
        element.label = label;
        elements.push_back(element);
    }
}

ImU32 color_cn() {
    return IM_COL32(255, 215, 0, 255);
}

ImU32 color_sigma_fill() {
    return IM_COL32(135, 206, 235, 38);
}

ImU32 color_sigma_line() {
    return IM_COL32(135, 206, 235, 220);
}

ImU32 color_inversion() {
    return IM_COL32(255, 0, 255, 255);
}

ImU32 color_sn() {
    return IM_COL32(0, 255, 255, 255);
}

bool world_to_screen(const Eigen::Vector3f& world_pos,
                     const Eigen::Matrix4f& vp_matrix,
                     const ImVec2& viewport_pos,
                     const ImVec2& viewport_size,
                     ImVec2& screen_pos) {
    const Eigen::Vector4f clip = vp_matrix * Eigen::Vector4f(world_pos.x(), world_pos.y(), world_pos.z(), 1.0f);
    if (clip.w() <= kScreenEpsilon) {
        return false;
    }
    const Eigen::Vector3f ndc = clip.head<3>() / clip.w();
    if (ndc.z() < -1.5f || ndc.z() > 1.5f) {
        return false;
    }
    screen_pos.x = viewport_pos.x + (ndc.x() * 0.5f + 0.5f) * viewport_size.x;
    screen_pos.y = viewport_pos.y + (1.0f - (ndc.y() * 0.5f + 0.5f)) * viewport_size.y;
    return true;
}

void draw_arrowhead(ImDrawList* draw_list, const ImVec2& tip, const ImVec2& tail, ImU32 color) {
    const ImVec2 dir = ImVec2(tip.x - tail.x, tip.y - tail.y);
    const float len = std::sqrt(dir.x * dir.x + dir.y * dir.y);
    if (len < 1.0f) {
        return;
    }
    const ImVec2 unit = ImVec2(dir.x / len, dir.y / len);
    const ImVec2 perp = ImVec2(-unit.y, unit.x);
    const float head_length = 10.0f;
    const float head_width = 5.0f;
    const ImVec2 base = ImVec2(tip.x - unit.x * head_length, tip.y - unit.y * head_length);
    draw_list->AddTriangleFilled(
        tip,
        ImVec2(base.x + perp.x * head_width, base.y + perp.y * head_width),
        ImVec2(base.x - perp.x * head_width, base.y - perp.y * head_width),
        color);
}

void draw_dashed_line(ImDrawList* draw_list,
                      const ImVec2& start,
                      const ImVec2& end,
                      ImU32 color,
                      float thickness) {
    const ImVec2 delta = ImVec2(end.x - start.x, end.y - start.y);
    const float length = std::sqrt(delta.x * delta.x + delta.y * delta.y);
    if (length < 1.0f) {
        return;
    }
    const ImVec2 unit = ImVec2(delta.x / length, delta.y / length);
    const float dash = 10.0f;
    const float gap = 6.0f;
    for (float t = 0.0f; t < length; t += dash + gap) {
        const float next = std::min(length, t + dash);
        const ImVec2 a = ImVec2(start.x + unit.x * t, start.y + unit.y * t);
        const ImVec2 b = ImVec2(start.x + unit.x * next, start.y + unit.y * next);
        draw_list->AddLine(a, b, color, thickness);
    }
}

void draw_rotation_glyph(ImDrawList* draw_list,
                         const SymmetryElement& element,
                         const Eigen::Vector3f& center,
                         const Eigen::Matrix4f& vp_matrix,
                         const ImVec2& viewport_pos,
                         const ImVec2& viewport_size,
                         float molecule_radius,
                         ImU32 color) {
    if (element.order <= 1) {
        return;
    }

    Eigen::Vector3f axis = element.axis.cast<float>().normalized();
    Eigen::Vector3f u = axis.unitOrthogonal().normalized();
    Eigen::Vector3f v = axis.cross(u).normalized();
    const float arc_radius = std::max(0.25f * molecule_radius, 0.6f);
    const float angle = 2.0f * kPiF / static_cast<float>(element.order);
    const float start_angle = 0.25f * kPiF;
    const int segments = 20;

    std::vector<ImVec2> points;
    points.reserve(static_cast<std::size_t>(segments + 1));
    for (int i = 0; i <= segments; ++i) {
        const float t = static_cast<float>(i) / static_cast<float>(segments);
        const float theta = start_angle + angle * t;
        const Eigen::Vector3f world = center + arc_radius * (std::cos(theta) * u + std::sin(theta) * v);
        ImVec2 screen{};
        if (!world_to_screen(world, vp_matrix, viewport_pos, viewport_size, screen)) {
            return;
        }
        points.push_back(screen);
    }

    for (std::size_t i = 1; i < points.size(); ++i) {
        draw_list->AddLine(points[i - 1], points[i], color, 1.5f);
    }
    draw_arrowhead(draw_list, points.back(), points[points.size() - 2], color);
}

void draw_plane(ImDrawList* draw_list,
                const SymmetryElement& element,
                const Eigen::Matrix4f& vp_matrix,
                const ImVec2& viewport_pos,
                const ImVec2& viewport_size,
                float molecule_radius) {
    const Eigen::Vector3f normal = element.axis.cast<float>().normalized();
    Eigen::Vector3f basis_u = normal.unitOrthogonal().normalized();
    Eigen::Vector3f basis_v = normal.cross(basis_u).normalized();
    const Eigen::Vector3f center = element.position.cast<float>();
    const float half_extent = std::max(1.5f * molecule_radius, 1.0f);

    const std::array<Eigen::Vector3f, 4> corners = {{
        center + (-basis_u - basis_v) * half_extent,
        center + ( basis_u - basis_v) * half_extent,
        center + ( basis_u + basis_v) * half_extent,
        center + (-basis_u + basis_v) * half_extent,
    }};

    std::array<ImVec2, 4> screen{};
    for (std::size_t i = 0; i < corners.size(); ++i) {
        if (!world_to_screen(corners[i], vp_matrix, viewport_pos, viewport_size, screen[i])) {
            return;
        }
    }

    draw_list->AddQuadFilled(screen[0], screen[1], screen[2], screen[3], color_sigma_fill());
    draw_list->AddQuad(screen[0], screen[1], screen[2], screen[3], color_sigma_line(), 1.5f);
    draw_list->AddText(screen[2], color_sigma_line(), element.label.c_str());
}

}  // namespace

std::vector<SymmetryElement> extract_symmetry_elements(const sbox::chem::MolecularSystem& mol,
                                                       sbox::chem::PointGroup pg,
                                                       double tolerance) {
    std::vector<SymmetryElement> elements;
    if (mol.num_atoms() == 0) {
        return elements;
    }

    const CenteredMolecule centered = center_molecule(mol);
    const FrameAxes frame = compute_frame_axes(centered);
    const std::vector<Eigen::Vector3d> candidate_axes = build_candidate_axes(centered, frame);
    const std::vector<Eigen::Vector3d> candidate_normals = build_candidate_plane_normals(centered, frame);
    const std::vector<Eigen::Vector3d> c2_axes = find_cn_axes(centered, candidate_axes, 2, tolerance);
    const std::vector<Eigen::Vector3d> c3_axes = find_cn_axes(centered, candidate_axes, 3, tolerance);
    const std::vector<Eigen::Vector3d> c4_axes = find_cn_axes(centered, candidate_axes, 4, tolerance);
    const std::vector<Eigen::Vector3d> c6_axes = find_cn_axes(centered, candidate_axes, 6, tolerance);
    const std::vector<Eigen::Vector3d> s4_axes = find_sn_axes(centered, candidate_axes, 4, tolerance);
    const std::vector<Eigen::Vector3d> s6_axes = find_sn_axes(centered, candidate_axes, 6, tolerance);
    const std::vector<Eigen::Vector3d> sigma_normals = find_sigma_normals(centered, candidate_normals, tolerance);
    const bool has_inversion = check_inversion(centered, tolerance);

    auto vertical_planes = [&](const Eigen::Vector3d& axis) {
        auto planes = filter_axes(sigma_normals, [&](const Eigen::Vector3d& normal) {
            return std::abs(canonicalize_axis(normal).dot(canonicalize_axis(axis))) < 0.15;
        });
        return sort_by_reference(std::move(planes), frame.secondary, false);
    };

    auto horizontal_planes = [&](const Eigen::Vector3d& axis) {
        auto planes = filter_axes(sigma_normals, [&](const Eigen::Vector3d& normal) {
            return std::abs(canonicalize_axis(normal).dot(canonicalize_axis(axis))) > 0.9;
        });
        return sort_by_reference(std::move(planes), axis, true);
    };

    switch (pg) {
    case sbox::chem::PointGroup::C2v: {
        const std::vector<Eigen::Vector3d> primary = select_count(sort_by_reference(c2_axes, frame.principal, true), 1);
        if (!primary.empty()) {
            append_axes(elements, primary, SymmetryElement::Type::Cn, 2, "C", centered.center);
            append_planes(elements, select_count(vertical_planes(primary.front()), 2), "sigma_v", centered.center);
        }
        break;
    }
    case sbox::chem::PointGroup::C3v: {
        const std::vector<Eigen::Vector3d> primary = select_count(sort_by_reference(c3_axes, frame.principal, true), 1);
        if (!primary.empty()) {
            append_axes(elements, primary, SymmetryElement::Type::Cn, 3, "C", centered.center);
            append_planes(elements, select_count(vertical_planes(primary.front()), 3), "sigma_v", centered.center);
        }
        break;
    }
    case sbox::chem::PointGroup::D6h: {
        const std::vector<Eigen::Vector3d> primary = select_count(sort_by_reference(c6_axes, frame.principal, true), 1);
        if (!primary.empty()) {
            const Eigen::Vector3d axis = primary.front();
            append_axes(elements, primary, SymmetryElement::Type::Cn, 6, "C", centered.center);
            auto perp_c2 = filter_axes(c2_axes, [&](const Eigen::Vector3d& candidate) {
                return std::abs(canonicalize_axis(candidate).dot(canonicalize_axis(axis))) < 0.15;
            });
            perp_c2 = sort_by_reference(std::move(perp_c2), frame.secondary, false);
            append_axes(elements, select_count(perp_c2, 6), SymmetryElement::Type::Cn, 2, "C", centered.center);
            append_planes(elements, select_count(horizontal_planes(axis), 1), "sigma_h", centered.center);
            append_planes(elements, select_count(vertical_planes(axis), 6), "sigma_v", centered.center);
            append_axes(elements, select_count(sort_by_reference(s6_axes, axis, true), 1), SymmetryElement::Type::Sn, 6, "S", centered.center);
        }
        if (has_inversion) {
            SymmetryElement inversion;
            inversion.type = SymmetryElement::Type::Inversion;
            inversion.position = centered.center;
            inversion.label = "i";
            elements.push_back(inversion);
        }
        break;
    }
    case sbox::chem::PointGroup::Td: {
        append_axes(elements, select_count(sort_by_reference(c3_axes, frame.principal, false), 4), SymmetryElement::Type::Cn, 3, "C", centered.center);
        auto td_c2 = sort_by_reference(c2_axes, frame.secondary, false);
        append_axes(elements, select_count(td_c2, 3), SymmetryElement::Type::Cn, 2, "C", centered.center);
        append_planes(elements, select_count(sort_by_reference(sigma_normals, frame.principal, false), 6), "sigma_d", centered.center);
        auto td_s4 = sort_by_reference(s4_axes, frame.secondary, false);
        append_axes(elements, select_count(td_s4, 3), SymmetryElement::Type::Sn, 4, "S", centered.center);
        break;
    }
    case sbox::chem::PointGroup::Oh: {
        const auto oh_c4 = select_count(sort_by_reference(c4_axes, frame.secondary, false), 3);
        append_axes(elements, oh_c4, SymmetryElement::Type::Cn, 4, "C", centered.center);
        append_axes(elements, select_count(sort_by_reference(c3_axes, frame.principal, false), 4), SymmetryElement::Type::Cn, 3, "C", centered.center);
        auto filtered_c2 = filter_axes(c2_axes, [&](const Eigen::Vector3d& axis) {
            for (const Eigen::Vector3d& c4 : oh_c4) {
                if (axis_score(axis, c4) > 0.95) {
                    return false;
                }
            }
            return true;
        });
        append_axes(elements, select_count(sort_by_reference(std::move(filtered_c2), frame.tertiary, false), 6), SymmetryElement::Type::Cn, 2, "C", centered.center);
        append_planes(elements, select_count(sort_by_reference(sigma_normals, frame.secondary, false), 9), "sigma", centered.center);
        if (has_inversion) {
            SymmetryElement inversion;
            inversion.type = SymmetryElement::Type::Inversion;
            inversion.position = centered.center;
            inversion.label = "i";
            elements.push_back(inversion);
        }
        break;
    }
    default: {
        if (!c6_axes.empty()) {
            append_axes(elements, select_count(sort_by_reference(c6_axes, frame.principal, true), 1), SymmetryElement::Type::Cn, 6, "C", centered.center);
        } else if (!c4_axes.empty()) {
            append_axes(elements, select_count(sort_by_reference(c4_axes, frame.principal, true), 1), SymmetryElement::Type::Cn, 4, "C", centered.center);
        } else if (!c3_axes.empty()) {
            append_axes(elements, select_count(sort_by_reference(c3_axes, frame.principal, true), 1), SymmetryElement::Type::Cn, 3, "C", centered.center);
        } else if (!c2_axes.empty()) {
            append_axes(elements, select_count(sort_by_reference(c2_axes, frame.principal, true), 1), SymmetryElement::Type::Cn, 2, "C", centered.center);
        }
        if (!sigma_normals.empty()) {
            append_planes(elements, select_count(sort_by_reference(sigma_normals, frame.principal, false), 3), "sigma", centered.center);
        }
        if (has_inversion) {
            SymmetryElement inversion;
            inversion.type = SymmetryElement::Type::Inversion;
            inversion.position = centered.center;
            inversion.label = "i";
            elements.push_back(inversion);
        }
        break;
    }
    }

    return elements;
}

void draw_symmetry_overlays(const std::vector<SymmetryElement>& elements,
                            const sbox::chem::MolecularSystem& mol,
                            const Eigen::Matrix4f& view_matrix,
                            const Eigen::Matrix4f& proj_matrix,
                            const ImVec2& viewport_pos,
                            const ImVec2& viewport_size,
                            float molecule_radius) {
    if (elements.empty() || mol.num_atoms() == 0 || viewport_size.x <= 1.0f || viewport_size.y <= 1.0f) {
        return;
    }

    ImDrawList* draw_list = ImGui::GetForegroundDrawList();
    const Eigen::Matrix4f vp_matrix = proj_matrix * view_matrix;
    const Eigen::Vector3f center = mol.center_of_mass().cast<float>();
    const float line_extent = std::max(2.0f * molecule_radius, 2.0f);

    for (const SymmetryElement& element : elements) {
        if (element.type == SymmetryElement::Type::Sigma) {
            draw_plane(draw_list, element, vp_matrix, viewport_pos, viewport_size, molecule_radius);
            continue;
        }

        if (element.type == SymmetryElement::Type::Inversion) {
            ImVec2 screen{};
            if (!world_to_screen(element.position.cast<float>(), vp_matrix, viewport_pos, viewport_size, screen)) {
                continue;
            }
            draw_list->AddCircleFilled(screen, 6.0f, color_inversion());
            draw_list->AddText(ImVec2(screen.x + 8.0f, screen.y - 10.0f), color_inversion(), element.label.c_str());
            continue;
        }

        const Eigen::Vector3f axis = element.axis.cast<float>().normalized();
        const Eigen::Vector3f base = element.position.cast<float>();
        const Eigen::Vector3f start_world = base - axis * line_extent;
        const Eigen::Vector3f end_world = base + axis * line_extent;
        ImVec2 start_screen{};
        ImVec2 end_screen{};
        if (!world_to_screen(start_world, vp_matrix, viewport_pos, viewport_size, start_screen) ||
            !world_to_screen(end_world, vp_matrix, viewport_pos, viewport_size, end_screen)) {
            continue;
        }

        const ImU32 color = element.type == SymmetryElement::Type::Sn ? color_sn() : color_cn();
        if (element.type == SymmetryElement::Type::Sn) {
            draw_dashed_line(draw_list, start_screen, end_screen, color, 2.0f);
        } else {
            draw_list->AddLine(start_screen, end_screen, color, 2.5f);
        }
        draw_arrowhead(draw_list, end_screen, start_screen, color);
        draw_arrowhead(draw_list, start_screen, end_screen, color);
        draw_list->AddText(ImVec2(end_screen.x + 8.0f, end_screen.y + 4.0f), color, element.label.c_str());
        draw_rotation_glyph(draw_list, element, base, vp_matrix, viewport_pos, viewport_size, molecule_radius, color);
    }
}

}  // namespace sbox::ui
