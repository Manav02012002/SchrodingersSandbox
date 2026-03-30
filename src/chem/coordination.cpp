#include "chem/coordination.h"

#include "core/covalent_radii.h"

#include <Eigen/Geometry>

#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <tuple>

namespace sbox::chem {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kSqrt3 = 1.73205080756887729353;

Eigen::Vector3d unit(double x, double y, double z) {
    Eigen::Vector3d v(x, y, z);
    const double n = v.norm();
    if (n > 1.0e-12) {
        return v / n;
    }
    return Eigen::Vector3d::UnitZ();
}

std::vector<Eigen::Vector3d> pentagon_xy() {
    std::vector<Eigen::Vector3d> dirs;
    dirs.reserve(5);
    for (int i = 0; i < 5; ++i) {
        const double a = 2.0 * kPi * static_cast<double>(i) / 5.0;
        dirs.emplace_back(std::cos(a), std::sin(a), 0.0);
    }
    return dirs;
}

std::vector<Eigen::Vector3d> square_antiprismatic_dirs() {
    std::vector<Eigen::Vector3d> dirs;
    dirs.reserve(8);
    const double z = 0.5;
    const double r = std::sqrt(1.0 - z * z);
    for (int i = 0; i < 4; ++i) {
        const double a = 2.0 * kPi * static_cast<double>(i) / 4.0;
        dirs.push_back(unit(r * std::cos(a), r * std::sin(a), -z));
    }
    for (int i = 0; i < 4; ++i) {
        const double a = 2.0 * kPi * static_cast<double>(i) / 4.0 + kPi / 4.0;
        dirs.push_back(unit(r * std::cos(a), r * std::sin(a), z));
    }
    return dirs;
}

const std::vector<CoordinationTemplate>& templates() {
    static const std::vector<CoordinationTemplate> kTemplates = {
        {CoordinationGeometry::Linear, "Linear", 2,
         {Eigen::Vector3d(0.0, 0.0, 1.0), Eigen::Vector3d(0.0, 0.0, -1.0)}},
        {CoordinationGeometry::TrigonalPlanar, "Trigonal Planar", 3,
         {Eigen::Vector3d(1.0, 0.0, 0.0),
          Eigen::Vector3d(-0.5, kSqrt3 / 2.0, 0.0),
          Eigen::Vector3d(-0.5, -kSqrt3 / 2.0, 0.0)}},
        {CoordinationGeometry::TShaped, "T-Shaped", 3,
         {Eigen::Vector3d(1.0, 0.0, 0.0),
          Eigen::Vector3d(-1.0, 0.0, 0.0),
          Eigen::Vector3d(0.0, 0.0, 1.0)}},
        {CoordinationGeometry::Tetrahedral, "Tetrahedral", 4,
         {unit(1.0, 1.0, 1.0),
          unit(1.0, -1.0, -1.0),
          unit(-1.0, 1.0, -1.0),
          unit(-1.0, -1.0, 1.0)}},
        {CoordinationGeometry::SquarePlanar, "Square Planar", 4,
         {Eigen::Vector3d(1.0, 0.0, 0.0),
          Eigen::Vector3d(0.0, 1.0, 0.0),
          Eigen::Vector3d(-1.0, 0.0, 0.0),
          Eigen::Vector3d(0.0, -1.0, 0.0)}},
        {CoordinationGeometry::SeeSaw, "See-Saw", 4,
         {Eigen::Vector3d(0.0, 0.0, 1.0),
          Eigen::Vector3d(0.0, 0.0, -1.0),
          Eigen::Vector3d(1.0, 0.0, 0.0),
          Eigen::Vector3d(-0.5, kSqrt3 / 2.0, 0.0)}},
        {CoordinationGeometry::TrigonalBipyramidal, "Trigonal Bipyramidal", 5,
         {Eigen::Vector3d(0.0, 0.0, 1.0),
          Eigen::Vector3d(0.0, 0.0, -1.0),
          Eigen::Vector3d(1.0, 0.0, 0.0),
          Eigen::Vector3d(-0.5, kSqrt3 / 2.0, 0.0),
          Eigen::Vector3d(-0.5, -kSqrt3 / 2.0, 0.0)}},
        {CoordinationGeometry::SquarePyramidal, "Square Pyramidal", 5,
         {Eigen::Vector3d(0.0, 0.0, 1.0),
          Eigen::Vector3d(1.0, 0.0, 0.0),
          Eigen::Vector3d(0.0, 1.0, 0.0),
          Eigen::Vector3d(-1.0, 0.0, 0.0),
          Eigen::Vector3d(0.0, -1.0, 0.0)}},
        {CoordinationGeometry::Octahedral, "Octahedral", 6,
         {Eigen::Vector3d(1.0, 0.0, 0.0),
          Eigen::Vector3d(-1.0, 0.0, 0.0),
          Eigen::Vector3d(0.0, 1.0, 0.0),
          Eigen::Vector3d(0.0, -1.0, 0.0),
          Eigen::Vector3d(0.0, 0.0, 1.0),
          Eigen::Vector3d(0.0, 0.0, -1.0)}},
        {CoordinationGeometry::PentagonalBipyramidal, "Pentagonal Bipyramidal", 7,
         [] {
             std::vector<Eigen::Vector3d> dirs = {Eigen::Vector3d(0.0, 0.0, 1.0), Eigen::Vector3d(0.0, 0.0, -1.0)};
             const std::vector<Eigen::Vector3d> pent = pentagon_xy();
             dirs.insert(dirs.end(), pent.begin(), pent.end());
             return dirs;
         }()},
        {CoordinationGeometry::SquareAntiprismatic, "Square Antiprismatic", 8, square_antiprismatic_dirs()},
    };
    return kTemplates;
}

Eigen::Vector3d ligand_direction(const MolecularSystem& ligand) {
    if (ligand.num_atoms() <= 1) {
        return Eigen::Vector3d::UnitZ();
    }
    const Eigen::Vector3d donor = ligand.atom(0).position;
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    for (int i = 1; i < ligand.num_atoms(); ++i) {
        centroid += ligand.atom(i).position;
    }
    centroid /= static_cast<double>(ligand.num_atoms() - 1);
    Eigen::Vector3d dir = centroid - donor;
    if (dir.norm() < 1.0e-10) {
        dir = Eigen::Vector3d::UnitZ();
    }
    return dir.normalized();
}

double angle_deg(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    const double denom = a.norm() * b.norm();
    if (denom < 1.0e-12) {
        return 0.0;
    }
    const double c = std::clamp(a.dot(b) / denom, -1.0, 1.0);
    return std::acos(c) * 180.0 / kPi;
}

bool near(double value, double target, double tolerance) {
    return std::abs(value - target) <= tolerance;
}

}  // namespace

const CoordinationTemplate& get_template(CoordinationGeometry geom) {
    for (const CoordinationTemplate& tmpl : templates()) {
        if (tmpl.geometry == geom) {
            return tmpl;
        }
    }
    throw std::invalid_argument("Unknown coordination geometry");
}

double default_ml_bond_length(int metal_Z, int ligand_Z) {
    static const std::map<std::pair<int, int>, double> lengths = {
        {{26, 7}, 3.78}, {{26, 8}, 3.97}, {{26, 17}, 4.35}, {{26, 16}, 4.35},
        {{27, 7}, 3.70}, {{27, 8}, 3.89},
        {{28, 7}, 3.59}, {{28, 8}, 3.78},
        {{29, 7}, 3.78}, {{29, 8}, 3.70},
        {{30, 7}, 3.89}, {{30, 8}, 3.97},
        {{24, 7}, 3.89}, {{24, 8}, 3.97},
        {{25, 7}, 3.89}, {{25, 8}, 3.97},
    };
    const auto it = lengths.find({metal_Z, ligand_Z});
    if (it != lengths.end()) {
        return it->second;
    }
    return covalent_radius(metal_Z) + covalent_radius(ligand_Z);
}

Eigen::Matrix3d orient_ligand(const Eigen::Vector3d& current_direction,
                              const Eigen::Vector3d& target_direction) {
    Eigen::Vector3d from = current_direction.norm() > 1.0e-12 ? current_direction.normalized() : Eigen::Vector3d::UnitZ();
    Eigen::Vector3d to = target_direction.norm() > 1.0e-12 ? target_direction.normalized() : Eigen::Vector3d::UnitZ();

    const double dot = std::clamp(from.dot(to), -1.0, 1.0);
    if (dot > 1.0 - 1.0e-12) {
        return Eigen::Matrix3d::Identity();
    }
    if (dot < -1.0 + 1.0e-12) {
        Eigen::Vector3d axis = from.unitOrthogonal().normalized();
        return Eigen::AngleAxisd(kPi, axis).toRotationMatrix();
    }

    Eigen::Vector3d axis = from.cross(to);
    const double angle = std::acos(dot);
    return Eigen::AngleAxisd(angle, axis.normalized()).toRotationMatrix();
}

MolecularSystem assemble_complex(int metal_Z,
                                 int metal_charge,
                                 CoordinationGeometry geometry,
                                 const std::vector<LigandSpec>& ligands) {
    const CoordinationTemplate& tmpl = get_template(geometry);
    if (static_cast<int>(ligands.size()) != tmpl.coordination_number) {
        throw std::invalid_argument("Ligand count does not match coordination geometry");
    }

    MolecularSystem complex;
    complex.set_name("Coordination Complex");
    complex.set_charge(metal_charge);
    complex.add_atom({metal_Z, Eigen::Vector3d::Zero(), "", metal_charge});

    for (int i = 0; i < tmpl.coordination_number; ++i) {
        const LigandSpec& spec = ligands[static_cast<std::size_t>(i)];
        if (spec.ligand.num_atoms() == 0) {
            throw std::invalid_argument("Ligand must contain at least the donor atom");
        }

        const Eigen::Vector3d direction = tmpl.directions[static_cast<std::size_t>(i)].normalized();
        const double bond_length = spec.bond_length > 0.0 ? spec.bond_length : default_ml_bond_length(metal_Z, spec.donor_Z);
        const Eigen::Vector3d donor_target = direction * bond_length;
        const Eigen::Vector3d donor_origin = spec.ligand.atom(0).position;
        const Eigen::Matrix3d rotation = orient_ligand(ligand_direction(spec.ligand), direction);

        const int atom_offset = complex.num_atoms();
        for (int atom_idx = 0; atom_idx < spec.ligand.num_atoms(); ++atom_idx) {
            Atom atom = spec.ligand.atom(atom_idx);
            atom.position = donor_target + rotation * (atom.position - donor_origin);
            complex.add_atom(atom);
        }
        for (const Bond& bond : spec.ligand.bonds()) {
            complex.add_bond(atom_offset + bond.atom_i, atom_offset + bond.atom_j, bond.order);
        }
        complex.add_bond(0, atom_offset, BondOrder::Single);
    }

    return complex;
}

CoordinationGeometry detect_coordination_geometry(const MolecularSystem& mol,
                                                 int metal_index,
                                                 double tolerance) {
    const std::vector<int> neigh = mol.neighbors(metal_index);
    const int cn = static_cast<int>(neigh.size());
    if (cn <= 2) {
        if (cn == 2) {
            const double ang = mol.angle(neigh[0], metal_index, neigh[1]) * 180.0 / kPi;
            if (near(ang, 180.0, tolerance)) {
                return CoordinationGeometry::Linear;
            }
        }
        return CoordinationGeometry::Linear;
    }

    std::vector<Eigen::Vector3d> dirs;
    dirs.reserve(neigh.size());
    for (int idx : neigh) {
        Eigen::Vector3d v = mol.atom(idx).position - mol.atom(metal_index).position;
        if (v.norm() > 1.0e-12) {
            dirs.push_back(v.normalized());
        }
    }

    std::vector<double> angles;
    for (std::size_t i = 0; i < dirs.size(); ++i) {
        for (std::size_t j = i + 1; j < dirs.size(); ++j) {
            angles.push_back(angle_deg(dirs[i], dirs[j]));
        }
    }

    auto count_near = [&](double target) {
        return static_cast<int>(std::count_if(angles.begin(), angles.end(), [&](double a) { return near(a, target, tolerance); }));
    };

    if (cn == 3) {
        const int linear_pairs = count_near(180.0);
        return linear_pairs >= 1 ? CoordinationGeometry::TShaped : CoordinationGeometry::TrigonalPlanar;
    }

    if (cn == 4) {
        const int near_90 = count_near(90.0);
        const int near_180 = count_near(180.0);
        const int near_109 = count_near(109.47);

        Eigen::Vector3d normal = Eigen::Vector3d::Zero();
        if (dirs.size() >= 3) {
            normal = (dirs[1] - dirs[0]).cross(dirs[2] - dirs[0]);
        }
        bool coplanar = normal.norm() > 1.0e-8;
        if (coplanar) {
            normal.normalize();
            for (const Eigen::Vector3d& dir : dirs) {
                if (std::abs(dir.dot(normal)) > 0.3) {
                    coplanar = false;
                    break;
                }
            }
        }

        if (near_90 >= 4 && near_180 >= 2 && coplanar) {
            return CoordinationGeometry::SquarePlanar;
        }
        if (near_109 >= 4) {
            return CoordinationGeometry::Tetrahedral;
        }
        return CoordinationGeometry::Tetrahedral;
    }

    if (cn == 5) {
        const int near_120 = count_near(120.0);
        const int near_180 = count_near(180.0);
        if (near_120 >= 3 && near_180 >= 1) {
            return CoordinationGeometry::TrigonalBipyramidal;
        }
        return CoordinationGeometry::SquarePyramidal;
    }

    if (cn == 6) {
        const int near_90 = count_near(90.0);
        const int near_180 = count_near(180.0);
        if (near_90 >= 8 && near_180 >= 3) {
            return CoordinationGeometry::Octahedral;
        }
        return CoordinationGeometry::Octahedral;
    }

    if (cn == 7) {
        return CoordinationGeometry::PentagonalBipyramidal;
    }
    if (cn >= 8) {
        return CoordinationGeometry::SquareAntiprismatic;
    }
    return CoordinationGeometry::Octahedral;
}

}  // namespace sbox::chem
