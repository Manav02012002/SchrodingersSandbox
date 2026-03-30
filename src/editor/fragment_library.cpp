#include "editor/fragment_library.h"

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace sbox::editor {

namespace {

using sbox::chem::Atom;
using sbox::chem::BondOrder;
using sbox::chem::MolecularSystem;

constexpr double kPi = 3.14159265358979323846;
constexpr double kCH = 2.06;
constexpr double kNH = 1.91;
constexpr double kOH = 1.81;
constexpr double kSH = 2.53;
constexpr double kCN = 2.21;
constexpr double kCO = 2.26;

int add_atom(MolecularSystem& mol, int Z, const Eigen::Vector3d& p, int charge = 0) {
    return mol.add_atom(Atom{Z, p, "", charge});
}

void add_methyl_hydrogens(MolecularSystem& mol, int carbon_index, const Eigen::Vector3d& away_from_attachment) {
    Eigen::Vector3d axis = away_from_attachment.norm() > 1.0e-8 ? away_from_attachment.normalized() : Eigen::Vector3d::UnitZ();
    Eigen::Vector3d ref = std::abs(axis.z()) < 0.9 ? Eigen::Vector3d::UnitZ() : Eigen::Vector3d::UnitX();
    Eigen::Vector3d u = axis.cross(ref).normalized();
    Eigen::Vector3d v = axis.cross(u).normalized();
    const double theta = std::acos(-1.0 / 3.0);
    for (int i = 0; i < 3; ++i) {
        const double phi = 2.0 * kPi * static_cast<double>(i) / 3.0;
        const Eigen::Vector3d dir = std::cos(theta) * axis + std::sin(theta) * (std::cos(phi) * u + std::sin(phi) * v);
        const int h = add_atom(mol, 1, mol.atom(carbon_index).position + kCH * dir);
        mol.add_bond(carbon_index, h, BondOrder::Single);
    }
}

Fragment make_methyl() {
    Fragment f;
    f.name = "methyl";
    f.category = "Functional Groups";
    f.attachment_atom = 0;
    auto& mol = f.molecule;
    const int c = add_atom(mol, 6, Eigen::Vector3d::Zero());
    add_methyl_hydrogens(mol, c, Eigen::Vector3d(0.0, 0.0, -1.0));
    return f;
}

Fragment make_hydroxyl() {
    Fragment f;
    f.name = "hydroxyl";
    f.category = "Functional Groups";
    f.attachment_atom = 0;
    auto& mol = f.molecule;
    const int o = add_atom(mol, 8, Eigen::Vector3d::Zero());
    const int h = add_atom(mol, 1, Eigen::Vector3d(0.0, 0.0, kOH));
    mol.add_bond(o, h, BondOrder::Single);
    return f;
}

Fragment make_amino() {
    Fragment f;
    f.name = "amino";
    f.category = "Functional Groups";
    f.attachment_atom = 0;
    auto& mol = f.molecule;
    const int n = add_atom(mol, 7, Eigen::Vector3d::Zero());
    const double theta = 107.0 * kPi / 180.0;
    const Eigen::Vector3d d1(std::sin(theta / 2.0), 0.0, std::cos(theta / 2.0));
    const Eigen::Vector3d d2(-std::sin(theta / 2.0), 0.0, std::cos(theta / 2.0));
    const int h1 = add_atom(mol, 1, kNH * d1);
    const int h2 = add_atom(mol, 1, kNH * d2);
    mol.add_bond(n, h1, BondOrder::Single);
    mol.add_bond(n, h2, BondOrder::Single);
    return f;
}

Fragment make_carboxyl() {
    Fragment f;
    f.name = "carboxyl";
    f.category = "Functional Groups";
    f.attachment_atom = 0;
    auto& mol = f.molecule;
    const int c = add_atom(mol, 6, Eigen::Vector3d::Zero());
    const int od = add_atom(mol, 8, Eigen::Vector3d(0.0, 0.0, kCO));
    const int os = add_atom(mol, 8, Eigen::Vector3d(kCO, 0.0, 0.0));
    const int h = add_atom(mol, 1, Eigen::Vector3d(kCO + kOH, 0.0, 0.0));
    mol.add_bond(c, od, BondOrder::Double);
    mol.add_bond(c, os, BondOrder::Single);
    mol.add_bond(os, h, BondOrder::Single);
    return f;
}

Fragment make_aldehyde() {
    Fragment f;
    f.name = "aldehyde";
    f.category = "Functional Groups";
    f.attachment_atom = 0;
    auto& mol = f.molecule;
    const int c = add_atom(mol, 6, Eigen::Vector3d::Zero());
    const int o = add_atom(mol, 8, Eigen::Vector3d(0.0, 0.0, kCO));
    const int h = add_atom(mol, 1, Eigen::Vector3d(kCH, 0.0, 0.0));
    mol.add_bond(c, o, BondOrder::Double);
    mol.add_bond(c, h, BondOrder::Single);
    return f;
}

Fragment make_nitro() {
    Fragment f;
    f.name = "nitro";
    f.category = "Functional Groups";
    f.attachment_atom = 0;
    auto& mol = f.molecule;
    const int n = add_atom(mol, 7, Eigen::Vector3d::Zero(), 1);
    const int o1 = add_atom(mol, 8, Eigen::Vector3d(0.0, 0.0, kCO));
    const int o2 = add_atom(mol, 8, Eigen::Vector3d(kCO * std::cos(2.0 * kPi / 3.0), kCO * std::sin(2.0 * kPi / 3.0), 0.0), -1);
    mol.add_bond(n, o1, BondOrder::Double);
    mol.add_bond(n, o2, BondOrder::Single);
    return f;
}

Fragment make_cyano() {
    Fragment f;
    f.name = "cyano";
    f.category = "Functional Groups";
    f.attachment_atom = 0;
    auto& mol = f.molecule;
    const int c = add_atom(mol, 6, Eigen::Vector3d::Zero());
    const int n = add_atom(mol, 7, Eigen::Vector3d(0.0, 0.0, kCN));
    mol.add_bond(c, n, BondOrder::Triple);
    return f;
}

Fragment make_thiol() {
    Fragment f;
    f.name = "thiol";
    f.category = "Functional Groups";
    f.attachment_atom = 0;
    auto& mol = f.molecule;
    const int s = add_atom(mol, 16, Eigen::Vector3d::Zero());
    const int h = add_atom(mol, 1, Eigen::Vector3d(0.0, 0.0, kSH));
    mol.add_bond(s, h, BondOrder::Single);
    return f;
}

Fragment make_benzene() {
    Fragment f;
    f.name = "benzene";
    f.category = "Ring Systems";
    auto& mol = f.molecule;
    const double r = 2.64;
    const double rh = r + kCH;
    for (int i = 0; i < 6; ++i) {
        const double a = i * kPi / 3.0;
        add_atom(mol, 6, Eigen::Vector3d(r * std::cos(a), r * std::sin(a), 0.0));
    }
    for (int i = 0; i < 6; ++i) {
        const double a = i * kPi / 3.0;
        add_atom(mol, 1, Eigen::Vector3d(rh * std::cos(a), rh * std::sin(a), 0.0));
    }
    for (int i = 0; i < 6; ++i) {
        mol.add_bond(i, (i + 1) % 6, BondOrder::Aromatic);
        mol.add_bond(i, 6 + i, BondOrder::Single);
    }
    return f;
}

Fragment make_cyclohexane() {
    Fragment f;
    f.name = "cyclohexane";
    f.category = "Ring Systems";
    auto& mol = f.molecule;
    const std::vector<Eigen::Vector3d> carbons = {
        {2.6, 0.0, 1.0}, {1.3, 2.25, -1.0}, {-1.3, 2.25, 1.0},
        {-2.6, 0.0, -1.0}, {-1.3, -2.25, 1.0}, {1.3, -2.25, -1.0}
    };
    for (const auto& p : carbons) add_atom(mol, 6, p);
    for (int i = 0; i < 6; ++i) mol.add_bond(i, (i + 1) % 6, BondOrder::Single);
    for (int i = 0; i < 6; ++i) {
        const Eigen::Vector3d radial = carbons[i].normalized();
        const Eigen::Vector3d axial = (i % 2 == 0) ? Eigen::Vector3d(0.0, 0.0, 1.0) : Eigen::Vector3d(0.0, 0.0, -1.0);
        int h1 = add_atom(mol, 1, carbons[i] + kCH * radial);
        int h2 = add_atom(mol, 1, carbons[i] + kCH * axial);
        mol.add_bond(i, h1, BondOrder::Single);
        mol.add_bond(i, h2, BondOrder::Single);
    }
    return f;
}

Fragment make_cyclopentadienyl() {
    Fragment f;
    f.name = "cyclopentadienyl";
    f.category = "Ring Systems";
    auto& mol = f.molecule;
    const double r = 2.45;
    const double rh = r + kCH;
    for (int i = 0; i < 5; ++i) {
        const double a = 2.0 * kPi * static_cast<double>(i) / 5.0;
        add_atom(mol, 6, Eigen::Vector3d(r * std::cos(a), r * std::sin(a), 0.0));
    }
    for (int i = 0; i < 5; ++i) {
        const double a = 2.0 * kPi * static_cast<double>(i) / 5.0;
        add_atom(mol, 1, Eigen::Vector3d(rh * std::cos(a), rh * std::sin(a), 0.0));
    }
    for (int i = 0; i < 5; ++i) {
        mol.add_bond(i, (i + 1) % 5, BondOrder::Aromatic);
        mol.add_bond(i, 5 + i, BondOrder::Single);
    }
    return f;
}

Fragment make_pyridine() {
    Fragment f;
    f.name = "pyridine";
    f.category = "Ring Systems";
    auto& mol = f.molecule;
    const double r = 2.62;
    const double rh = r + kCH;
    for (int i = 0; i < 6; ++i) {
        const double a = i * kPi / 3.0;
        add_atom(mol, i == 0 ? 7 : 6, Eigen::Vector3d(r * std::cos(a), r * std::sin(a), 0.0));
    }
    for (int i = 1; i < 6; ++i) {
        const double a = i * kPi / 3.0;
        add_atom(mol, 1, Eigen::Vector3d(rh * std::cos(a), rh * std::sin(a), 0.0));
    }
    for (int i = 0; i < 6; ++i) mol.add_bond(i, (i + 1) % 6, BondOrder::Aromatic);
    for (int i = 1; i < 6; ++i) mol.add_bond(i, 5 + i - 1, BondOrder::Single);
    return f;
}

Fragment make_imidazole() {
    Fragment f;
    f.name = "imidazole";
    f.category = "Ring Systems";
    auto& mol = f.molecule;
    const std::vector<std::pair<int, Eigen::Vector3d>> ring = {
        {7, {2.2, 0.0, 0.0}}, {6, {0.7, 2.1, 0.0}}, {7, {-1.8, 1.3, 0.0}},
        {6, {-1.8, -1.3, 0.0}}, {6, {0.7, -2.1, 0.0}}
    };
    for (const auto& [z, p] : ring) add_atom(mol, z, p);
    for (int i = 0; i < 5; ++i) mol.add_bond(i, (i + 1) % 5, BondOrder::Aromatic);
    add_atom(mol, 1, Eigen::Vector3d(3.7, 0.0, 0.0));
    add_atom(mol, 1, Eigen::Vector3d(1.2, 3.7, 0.0));
    add_atom(mol, 1, Eigen::Vector3d(-3.3, -2.4, 0.0));
    add_atom(mol, 1, Eigen::Vector3d(1.2, -3.7, 0.0));
    mol.add_bond(0, 5, BondOrder::Single);
    mol.add_bond(1, 6, BondOrder::Single);
    mol.add_bond(3, 7, BondOrder::Single);
    mol.add_bond(4, 8, BondOrder::Single);
    return f;
}

Fragment make_water() {
    Fragment f;
    f.name = "water";
    f.category = "Ligands";
    f.attachment_atom = 0;
    auto& mol = f.molecule;
    const int o = add_atom(mol, 8, Eigen::Vector3d::Zero());
    const double angle = 104.5 * kPi / 180.0;
    const Eigen::Vector3d h1(0.0, std::sin(angle / 2.0) * kOH, -std::cos(angle / 2.0) * kOH);
    const Eigen::Vector3d h2(0.0, -std::sin(angle / 2.0) * kOH, -std::cos(angle / 2.0) * kOH);
    const int a = add_atom(mol, 1, h1);
    const int b = add_atom(mol, 1, h2);
    mol.add_bond(o, a, BondOrder::Single);
    mol.add_bond(o, b, BondOrder::Single);
    return f;
}

Fragment make_ammonia() {
    Fragment f;
    f.name = "ammonia";
    f.category = "Ligands";
    f.attachment_atom = 0;
    auto& mol = f.molecule;
    const int n = add_atom(mol, 7, Eigen::Vector3d::Zero());
    add_methyl_hydrogens(mol, n, Eigen::Vector3d(0.0, 0.0, -1.0));
    return f;
}

Fragment make_carbon_monoxide() {
    Fragment f;
    f.name = "carbon monoxide";
    f.category = "Ligands";
    f.attachment_atom = 0;
    auto& mol = f.molecule;
    const int c = add_atom(mol, 6, Eigen::Vector3d::Zero());
    const int o = add_atom(mol, 8, Eigen::Vector3d(0.0, 0.0, kCO));
    mol.add_bond(c, o, BondOrder::Triple);
    return f;
}

Fragment make_cyanide() {
    Fragment f;
    f.name = "cyanide";
    f.category = "Ligands";
    f.attachment_atom = 0;
    auto& mol = f.molecule;
    mol.set_charge(-1);
    const int c = add_atom(mol, 6, Eigen::Vector3d::Zero(), -1);
    const int n = add_atom(mol, 7, Eigen::Vector3d(0.0, 0.0, kCN));
    mol.add_bond(c, n, BondOrder::Triple);
    return f;
}

Fragment make_en() {
    Fragment f;
    f.name = "ethylenediamine";
    f.category = "Ligands";
    f.attachment_atom = 0;
    auto& mol = f.molecule;
    const int n1 = add_atom(mol, 7, Eigen::Vector3d(-3.0, 0.0, 0.0));
    const int c1 = add_atom(mol, 6, Eigen::Vector3d(-1.2, 0.0, 0.0));
    const int c2 = add_atom(mol, 6, Eigen::Vector3d(1.2, 0.0, 0.0));
    const int n2 = add_atom(mol, 7, Eigen::Vector3d(3.0, 0.0, 0.0));
    mol.add_bond(n1, c1, BondOrder::Single);
    mol.add_bond(c1, c2, BondOrder::Single);
    mol.add_bond(c2, n2, BondOrder::Single);
    add_atom(mol, 1, Eigen::Vector3d(-3.8, 1.4, 0.0));
    add_atom(mol, 1, Eigen::Vector3d(-3.8, -1.4, 0.0));
    add_atom(mol, 1, Eigen::Vector3d(-1.2, 1.8, 0.9));
    add_atom(mol, 1, Eigen::Vector3d(-1.2, -1.8, -0.9));
    add_atom(mol, 1, Eigen::Vector3d(1.2, 1.8, -0.9));
    add_atom(mol, 1, Eigen::Vector3d(1.2, -1.8, 0.9));
    add_atom(mol, 1, Eigen::Vector3d(3.8, 1.4, 0.0));
    add_atom(mol, 1, Eigen::Vector3d(3.8, -1.4, 0.0));
    mol.add_bond(n1, 4, BondOrder::Single);
    mol.add_bond(n1, 5, BondOrder::Single);
    mol.add_bond(c1, 6, BondOrder::Single);
    mol.add_bond(c1, 7, BondOrder::Single);
    mol.add_bond(c2, 8, BondOrder::Single);
    mol.add_bond(c2, 9, BondOrder::Single);
    mol.add_bond(n2, 10, BondOrder::Single);
    mol.add_bond(n2, 11, BondOrder::Single);
    return f;
}

Fragment make_glycine() {
    Fragment f;
    f.name = "glycine";
    f.category = "Amino Acids";
    auto& mol = f.molecule;
    const int n = add_atom(mol, 7, Eigen::Vector3d(-4.0, 0.0, 0.0));
    const int ca = add_atom(mol, 6, Eigen::Vector3d(-1.8, 0.0, 0.0));
    const int c = add_atom(mol, 6, Eigen::Vector3d(0.8, 0.0, 0.0));
    const int od = add_atom(mol, 8, Eigen::Vector3d(2.8, 1.0, 0.0));
    const int os = add_atom(mol, 8, Eigen::Vector3d(2.8, -1.0, 0.0));
    const int ho = add_atom(mol, 1, Eigen::Vector3d(4.2, -1.0, 0.0));
    mol.add_bond(n, ca, BondOrder::Single);
    mol.add_bond(ca, c, BondOrder::Single);
    mol.add_bond(c, od, BondOrder::Double);
    mol.add_bond(c, os, BondOrder::Single);
    mol.add_bond(os, ho, BondOrder::Single);
    add_atom(mol, 1, Eigen::Vector3d(-4.8, 1.4, 0.0));
    add_atom(mol, 1, Eigen::Vector3d(-4.8, -1.4, 0.0));
    add_atom(mol, 1, Eigen::Vector3d(-1.8, 1.5, 1.0));
    add_atom(mol, 1, Eigen::Vector3d(-1.8, -1.5, -1.0));
    mol.add_bond(n, 6, BondOrder::Single);
    mol.add_bond(n, 7, BondOrder::Single);
    mol.add_bond(ca, 8, BondOrder::Single);
    mol.add_bond(ca, 9, BondOrder::Single);
    return f;
}

Fragment make_alanine() {
    Fragment f;
    f.name = "alanine";
    f.category = "Amino Acids";
    auto& mol = f.molecule;
    const int n = add_atom(mol, 7, Eigen::Vector3d(-4.2, 0.0, 0.0));
    const int ca = add_atom(mol, 6, Eigen::Vector3d(-2.0, 0.0, 0.0));
    const int c = add_atom(mol, 6, Eigen::Vector3d(0.6, 0.0, 0.0));
    const int cb = add_atom(mol, 6, Eigen::Vector3d(-2.0, 2.4, 0.0));
    const int od = add_atom(mol, 8, Eigen::Vector3d(2.6, 1.0, 0.0));
    const int os = add_atom(mol, 8, Eigen::Vector3d(2.6, -1.0, 0.0));
    const int ho = add_atom(mol, 1, Eigen::Vector3d(4.0, -1.0, 0.0));
    mol.add_bond(n, ca, BondOrder::Single);
    mol.add_bond(ca, c, BondOrder::Single);
    mol.add_bond(ca, cb, BondOrder::Single);
    mol.add_bond(c, od, BondOrder::Double);
    mol.add_bond(c, os, BondOrder::Single);
    mol.add_bond(os, ho, BondOrder::Single);
    add_atom(mol, 1, Eigen::Vector3d(-5.0, 1.4, 0.0));
    add_atom(mol, 1, Eigen::Vector3d(-5.0, -1.4, 0.0));
    add_atom(mol, 1, Eigen::Vector3d(-2.0, -1.6, 1.0));
    add_atom(mol, 1, Eigen::Vector3d(-0.6, 3.2, 0.8));
    add_atom(mol, 1, Eigen::Vector3d(-2.0, 2.4, -1.8));
    add_atom(mol, 1, Eigen::Vector3d(-3.4, 3.2, 0.8));
    mol.add_bond(n, 7, BondOrder::Single);
    mol.add_bond(n, 8, BondOrder::Single);
    mol.add_bond(ca, 9, BondOrder::Single);
    mol.add_bond(cb, 10, BondOrder::Single);
    mol.add_bond(cb, 11, BondOrder::Single);
    mol.add_bond(cb, 12, BondOrder::Single);
    return f;
}

Fragment make_methanol() {
    Fragment f;
    f.name = "methanol";
    f.category = "Solvents";
    auto& mol = f.molecule;
    const int c = add_atom(mol, 6, Eigen::Vector3d::Zero());
    const int o = add_atom(mol, 8, Eigen::Vector3d(0.0, 0.0, 2.7));
    mol.add_bond(c, o, BondOrder::Single);
    add_methyl_hydrogens(mol, c, Eigen::Vector3d(0.0, 0.0, -1.0));
    const int h = add_atom(mol, 1, Eigen::Vector3d(0.0, 0.0, 2.7 + kOH));
    mol.add_bond(o, h, BondOrder::Single);
    return f;
}

Fragment make_dmso() {
    Fragment f;
    f.name = "dmso";
    f.category = "Solvents";
    auto& mol = f.molecule;
    const int s = add_atom(mol, 16, Eigen::Vector3d::Zero());
    const int o = add_atom(mol, 8, Eigen::Vector3d(0.0, 0.0, 2.8));
    const int c1 = add_atom(mol, 6, Eigen::Vector3d(2.9, 0.0, -0.8));
    const int c2 = add_atom(mol, 6, Eigen::Vector3d(-2.9, 0.0, -0.8));
    mol.add_bond(s, o, BondOrder::Double);
    mol.add_bond(s, c1, BondOrder::Single);
    mol.add_bond(s, c2, BondOrder::Single);
    add_methyl_hydrogens(mol, c1, Eigen::Vector3d(1.0, 0.0, -0.2));
    add_methyl_hydrogens(mol, c2, Eigen::Vector3d(-1.0, 0.0, -0.2));
    return f;
}

}  // namespace

FragmentLibrary::FragmentLibrary() {
    build_library();
}

const std::vector<Fragment>& FragmentLibrary::all() const {
    return fragments_;
}

std::vector<const Fragment*> FragmentLibrary::by_category(const std::string& category) const {
    std::vector<const Fragment*> out;
    for (const auto& fragment : fragments_) {
        if (fragment.category == category) {
            out.push_back(&fragment);
        }
    }
    return out;
}

const Fragment* FragmentLibrary::find(const std::string& name) const {
    for (const auto& fragment : fragments_) {
        if (fragment.name == name) {
            return &fragment;
        }
    }
    return nullptr;
}

std::vector<std::string> FragmentLibrary::categories() const {
    std::vector<std::string> out;
    for (const auto& fragment : fragments_) {
        if (std::find(out.begin(), out.end(), fragment.category) == out.end()) {
            out.push_back(fragment.category);
        }
    }
    return out;
}

sbox::chem::MolecularSystem FragmentLibrary::place(const Fragment& fragment,
                                                   const Eigen::Vector3d& position) const {
    MolecularSystem placed = fragment.molecule;
    Eigen::Vector3d anchor = Eigen::Vector3d::Zero();
    if (fragment.attachment_atom >= 0 && fragment.attachment_atom < placed.num_atoms()) {
        anchor = placed.atom(fragment.attachment_atom).position;
    }
    const Eigen::Vector3d shift = position - anchor;
    for (int i = 0; i < placed.num_atoms(); ++i) {
        placed.atom(i).position += shift;
    }
    return placed;
}

void FragmentLibrary::build_library() {
    fragments_.clear();
    fragments_.push_back(make_methyl());
    fragments_.push_back(make_hydroxyl());
    fragments_.push_back(make_amino());
    fragments_.push_back(make_carboxyl());
    fragments_.push_back(make_aldehyde());
    fragments_.push_back(make_nitro());
    fragments_.push_back(make_cyano());
    fragments_.push_back(make_thiol());
    fragments_.push_back(make_benzene());
    fragments_.push_back(make_cyclohexane());
    fragments_.push_back(make_cyclopentadienyl());
    fragments_.push_back(make_pyridine());
    fragments_.push_back(make_imidazole());
    fragments_.push_back(make_glycine());
    fragments_.push_back(make_alanine());
    fragments_.push_back(make_water());
    fragments_.push_back(make_ammonia());
    fragments_.push_back(make_carbon_monoxide());
    fragments_.push_back(make_cyanide());
    fragments_.push_back(make_en());
    fragments_.push_back(make_methanol());
    fragments_.push_back(make_dmso());
}

}  // namespace sbox::editor
