#include "chem/ligand_library.h"

#include "core/elements.h"

#include <Eigen/Geometry>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <stdexcept>

namespace sbox::chem {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kCH = 2.06;
constexpr double kNH = 1.91;
constexpr double kOH = 1.81;
constexpr double kPH = 2.67;
constexpr double kCN = 2.21;
constexpr double kCO = 2.13;
constexpr double kCC = 2.89;
constexpr double kNO = 2.30;
constexpr double kNN = 2.35;
constexpr double kCS = 3.05;
constexpr double kCO_SINGLE = 2.70;
constexpr double kCO_DOUBLE = 2.36;

int add_atom(MolecularSystem& mol, int Z, const Eigen::Vector3d& pos, int charge = 0) {
    return mol.add_atom(Atom{Z, pos, "", charge});
}

std::string lower_copy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

Eigen::Vector3d tetra_dir(int i) {
    static const Eigen::Vector3d dirs[] = {
        Eigen::Vector3d(1.0, 1.0, 1.0).normalized(),
        Eigen::Vector3d(1.0, -1.0, -1.0).normalized(),
        Eigen::Vector3d(-1.0, 1.0, -1.0).normalized(),
        Eigen::Vector3d(-1.0, -1.0, 1.0).normalized(),
    };
    return dirs[i % 4];
}

void add_trigonal_hydrogens(MolecularSystem& mol, int center, double bond_length, const Eigen::Vector3d& away) {
    const Eigen::Vector3d axis = away.norm() > 1.0e-8 ? away.normalized() : Eigen::Vector3d::UnitZ();
    const Eigen::Vector3d ref = std::abs(axis.z()) < 0.9 ? Eigen::Vector3d::UnitZ() : Eigen::Vector3d::UnitX();
    const Eigen::Vector3d u = axis.cross(ref).normalized();
    const Eigen::Vector3d v = axis.cross(u).normalized();
    const double theta = std::acos(-1.0 / 3.0);
    for (int i = 0; i < 3; ++i) {
        const double phi = 2.0 * kPi * static_cast<double>(i) / 3.0;
        const Eigen::Vector3d dir = std::cos(theta) * axis + std::sin(theta) * (std::cos(phi) * u + std::sin(phi) * v);
        const int h = add_atom(mol, 1, mol.atom(center).position + bond_length * dir);
        mol.add_bond(center, h, BondOrder::Single);
    }
}

LigandEntry make_water() {
    LigandEntry e;
    e.name = "Water";
    e.abbreviation = "h2o";
    e.formula = "H₂O";
    e.category = "Neutral";
    e.denticity = LigandDenticity::Monodentate;
    e.charge = 0;
    e.field_strength = "weak-moderate";
    e.typical_delta_oct_eV = 1.2;
    auto& mol = e.molecule;
    const int o = add_atom(mol, 8, Eigen::Vector3d::Zero());
    const double half = 104.5 * kPi / 360.0;
    add_atom(mol, 1, Eigen::Vector3d(kOH * std::sin(half), 0.0, kOH * std::cos(half)));
    add_atom(mol, 1, Eigen::Vector3d(-kOH * std::sin(half), 0.0, kOH * std::cos(half)));
    mol.add_bond(o, 1, BondOrder::Single);
    mol.add_bond(o, 2, BondOrder::Single);
    e.donor_atoms = {0};
    return e;
}

LigandEntry make_ammonia() {
    LigandEntry e;
    e.name = "Ammonia";
    e.abbreviation = "nh3";
    e.formula = "NH₃";
    e.category = "Neutral";
    e.denticity = LigandDenticity::Monodentate;
    e.field_strength = "moderate-strong";
    e.typical_delta_oct_eV = 1.5;
    auto& mol = e.molecule;
    const int n = add_atom(mol, 7, Eigen::Vector3d::Zero());
    add_trigonal_hydrogens(mol, n, kNH, Eigen::Vector3d(0.0, 0.0, -1.0));
    e.donor_atoms = {0};
    return e;
}

LigandEntry make_carbon_monoxide() {
    LigandEntry e;
    e.name = "Carbon Monoxide";
    e.abbreviation = "co";
    e.formula = "CO";
    e.category = "Neutral";
    e.denticity = LigandDenticity::Monodentate;
    e.field_strength = "very strong";
    e.typical_delta_oct_eV = 2.1;
    auto& mol = e.molecule;
    add_atom(mol, 6, Eigen::Vector3d::Zero());
    add_atom(mol, 8, Eigen::Vector3d(0.0, 0.0, kCO));
    mol.add_bond(0, 1, BondOrder::Triple);
    e.donor_atoms = {0};
    return e;
}

LigandEntry make_phosphine() {
    LigandEntry e;
    e.name = "Phosphine";
    e.abbreviation = "ph3";
    e.formula = "PH₃";
    e.category = "Neutral";
    e.denticity = LigandDenticity::Monodentate;
    e.field_strength = "strong";
    e.typical_delta_oct_eV = 1.7;
    auto& mol = e.molecule;
    const int p = add_atom(mol, 15, Eigen::Vector3d::Zero());
    add_trigonal_hydrogens(mol, p, kPH, Eigen::Vector3d(0.0, 0.0, -1.0));
    e.donor_atoms = {0};
    return e;
}

LigandEntry make_pyridine() {
    LigandEntry e;
    e.name = "Pyridine";
    e.abbreviation = "py";
    e.formula = "C₅H₅N";
    e.category = "Neutral";
    e.denticity = LigandDenticity::Monodentate;
    e.field_strength = "moderate";
    e.typical_delta_oct_eV = 1.4;
    auto& mol = e.molecule;
    const double r = 2.62;
    const double rh = r + kCH;
    for (int i = 0; i < 6; ++i) {
        const double a = static_cast<double>(i) * kPi / 3.0;
        add_atom(mol, i == 0 ? 7 : 6, Eigen::Vector3d(r * std::cos(a), r * std::sin(a), 0.0));
    }
    for (int i = 1; i < 6; ++i) {
        const double a = static_cast<double>(i) * kPi / 3.0;
        add_atom(mol, 1, Eigen::Vector3d(rh * std::cos(a), rh * std::sin(a), 0.0));
    }
    for (int i = 0; i < 6; ++i) mol.add_bond(i, (i + 1) % 6, BondOrder::Aromatic);
    for (int i = 1; i < 6; ++i) mol.add_bond(i, 5 + i - 1, BondOrder::Single);
    e.donor_atoms = {0};
    return e;
}

LigandEntry make_acetonitrile() {
    LigandEntry e;
    e.name = "Acetonitrile";
    e.abbreviation = "mecn";
    e.formula = "CH₃CN";
    e.category = "Neutral";
    e.denticity = LigandDenticity::Monodentate;
    e.field_strength = "moderate-strong";
    e.typical_delta_oct_eV = 1.45;
    auto& mol = e.molecule;
    const int n = add_atom(mol, 7, Eigen::Vector3d::Zero());
    const int c_sp = add_atom(mol, 6, Eigen::Vector3d(0.0, 0.0, kCN));
    const int c_me = add_atom(mol, 6, Eigen::Vector3d(0.0, 0.0, kCN + kCC));
    mol.add_bond(0, 1, BondOrder::Triple);
    mol.add_bond(1, 2, BondOrder::Single);
    for (int i = 0; i < 3; ++i) {
        const Eigen::Vector3d dir = tetra_dir(i);
        const int h = add_atom(mol, 1, mol.atom(c_me).position + kCH * dir);
        mol.add_bond(c_me, h, BondOrder::Single);
    }
    e.donor_atoms = {0};
    return e;
}

LigandEntry make_halide(const std::string& name, const std::string& abbr, const std::string& formula, int Z,
                        const std::string& strength, double delta) {
    LigandEntry e;
    e.name = name;
    e.abbreviation = abbr;
    e.formula = formula;
    e.category = "Anionic";
    e.denticity = LigandDenticity::Monodentate;
    e.charge = -1;
    e.field_strength = strength;
    e.typical_delta_oct_eV = delta;
    add_atom(e.molecule, Z, Eigen::Vector3d::Zero(), -1);
    e.donor_atoms = {0};
    return e;
}

LigandEntry make_hydroxide() {
    LigandEntry e;
    e.name = "Hydroxide";
    e.abbreviation = "oh";
    e.formula = "OH⁻";
    e.category = "Anionic";
    e.denticity = LigandDenticity::Monodentate;
    e.charge = -1;
    e.field_strength = "weak-moderate";
    e.typical_delta_oct_eV = 1.3;
    add_atom(e.molecule, 8, Eigen::Vector3d::Zero(), -1);
    add_atom(e.molecule, 1, Eigen::Vector3d(0.0, 0.0, kOH));
    e.molecule.add_bond(0, 1, BondOrder::Single);
    e.donor_atoms = {0};
    return e;
}

LigandEntry make_cyanide() {
    LigandEntry e;
    e.name = "Cyanide";
    e.abbreviation = "cn";
    e.formula = "CN⁻";
    e.category = "Anionic";
    e.denticity = LigandDenticity::Monodentate;
    e.charge = -1;
    e.field_strength = "very strong";
    e.typical_delta_oct_eV = 2.0;
    add_atom(e.molecule, 6, Eigen::Vector3d::Zero(), -1);
    add_atom(e.molecule, 7, Eigen::Vector3d(0.0, 0.0, kCN));
    e.molecule.add_bond(0, 1, BondOrder::Triple);
    e.donor_atoms = {0};
    return e;
}

LigandEntry make_thiocyanate() {
    LigandEntry e;
    e.name = "Thiocyanate";
    e.abbreviation = "scn";
    e.formula = "SCN⁻";
    e.category = "Anionic";
    e.denticity = LigandDenticity::Monodentate;
    e.charge = -1;
    e.field_strength = "weak";
    e.typical_delta_oct_eV = 1.3;
    add_atom(e.molecule, 16, Eigen::Vector3d::Zero(), -1);
    add_atom(e.molecule, 6, Eigen::Vector3d(0.0, 0.0, kCS));
    add_atom(e.molecule, 7, Eigen::Vector3d(0.0, 0.0, kCS + kCN));
    e.molecule.add_bond(0, 1, BondOrder::Single);
    e.molecule.add_bond(1, 2, BondOrder::Triple);
    e.donor_atoms = {0, 2};
    return e;
}

LigandEntry make_nitrite() {
    LigandEntry e;
    e.name = "Nitrite";
    e.abbreviation = "no2";
    e.formula = "NO₂⁻";
    e.category = "Anionic";
    e.denticity = LigandDenticity::Monodentate;
    e.charge = -1;
    e.field_strength = "strong";
    e.typical_delta_oct_eV = 1.7;
    const int n = add_atom(e.molecule, 7, Eigen::Vector3d::Zero());
    const double angle = 115.0 * kPi / 180.0;
    const Eigen::Vector3d o1(kNO * std::sin(angle / 2.0), 0.0, kNO * std::cos(angle / 2.0));
    const Eigen::Vector3d o2(-kNO * std::sin(angle / 2.0), 0.0, kNO * std::cos(angle / 2.0));
    add_atom(e.molecule, 8, o1, 0);
    add_atom(e.molecule, 8, o2, -1);
    e.molecule.add_bond(n, 1, BondOrder::Double);
    e.molecule.add_bond(n, 2, BondOrder::Single);
    e.donor_atoms = {0};
    return e;
}

LigandEntry make_azide() {
    LigandEntry e;
    e.name = "Azide";
    e.abbreviation = "n3";
    e.formula = "N₃⁻";
    e.category = "Anionic";
    e.denticity = LigandDenticity::Monodentate;
    e.charge = -1;
    e.field_strength = "moderate";
    e.typical_delta_oct_eV = 1.35;
    add_atom(e.molecule, 7, Eigen::Vector3d::Zero(), -1);
    add_atom(e.molecule, 7, Eigen::Vector3d(0.0, 0.0, kNN));
    add_atom(e.molecule, 7, Eigen::Vector3d(0.0, 0.0, 2.0 * kNN));
    e.molecule.add_bond(0, 1, BondOrder::Double);
    e.molecule.add_bond(1, 2, BondOrder::Double);
    e.donor_atoms = {0};
    return e;
}

LigandEntry make_en() {
    LigandEntry e;
    e.name = "Ethylenediamine";
    e.abbreviation = "en";
    e.formula = "H₂NCH₂CH₂NH₂";
    e.category = "Polydentate";
    e.denticity = LigandDenticity::Bidentate;
    e.charge = 0;
    e.field_strength = "strong";
    e.typical_delta_oct_eV = 1.55;
    auto& mol = e.molecule;
    const int n1 = add_atom(mol, 7, Eigen::Vector3d(-3.0, 0.0, 0.0));
    const int c1 = add_atom(mol, 6, Eigen::Vector3d(-1.2, 0.0, 0.6));
    const int c2 = add_atom(mol, 6, Eigen::Vector3d(1.2, 0.0, -0.6));
    const int n2 = add_atom(mol, 7, Eigen::Vector3d(3.0, 0.0, 0.0));
    mol.add_bond(n1, c1, BondOrder::Single);
    mol.add_bond(c1, c2, BondOrder::Single);
    mol.add_bond(c2, n2, BondOrder::Single);
    add_atom(mol, 1, Eigen::Vector3d(-3.8, 1.4, 0.0));
    add_atom(mol, 1, Eigen::Vector3d(-3.8, -1.4, 0.0));
    add_atom(mol, 1, Eigen::Vector3d(-1.3, 1.7, 1.3));
    add_atom(mol, 1, Eigen::Vector3d(-1.0, -1.7, 1.2));
    add_atom(mol, 1, Eigen::Vector3d(1.3, 1.7, -1.3));
    add_atom(mol, 1, Eigen::Vector3d(1.0, -1.7, -1.2));
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
    e.donor_atoms = {0, 3};
    return e;
}

LigandEntry make_bipy() {
    LigandEntry e;
    e.name = "Bipyridine";
    e.abbreviation = "bipy";
    e.formula = "C₁₀H₈N₂";
    e.category = "Polydentate";
    e.denticity = LigandDenticity::Bidentate;
    e.charge = 0;
    e.field_strength = "strong";
    e.typical_delta_oct_eV = 1.6;
    auto& mol = e.molecule;
    const double r = 2.62;
    for (int ring = 0; ring < 2; ++ring) {
        const double xoff = ring == 0 ? -3.2 : 3.2;
        const double rot = ring == 0 ? 0.0 : kPi;
        for (int i = 0; i < 6; ++i) {
            const double a = rot + static_cast<double>(i) * kPi / 3.0;
            add_atom(mol, i == 0 ? 7 : 6, Eigen::Vector3d(xoff + r * std::cos(a), r * std::sin(a), 0.0));
        }
    }
    for (int i = 0; i < 6; ++i) mol.add_bond(i, (i + 1) % 6, BondOrder::Aromatic);
    for (int i = 6; i < 12; ++i) mol.add_bond(i, 6 + ((i - 6 + 1) % 6), BondOrder::Aromatic);
    mol.add_bond(3, 9, BondOrder::Single);
    e.donor_atoms = {0, 6};
    return e;
}

LigandEntry make_acac() {
    LigandEntry e;
    e.name = "Acetylacetonate";
    e.abbreviation = "acac";
    e.formula = "acac⁻";
    e.category = "Polydentate";
    e.denticity = LigandDenticity::Bidentate;
    e.charge = -1;
    e.field_strength = "moderate";
    e.typical_delta_oct_eV = 1.45;
    auto& mol = e.molecule;
    const int o1 = add_atom(mol, 8, Eigen::Vector3d(-5.0, 0.8, 0.0), -1);
    const int c1 = add_atom(mol, 6, Eigen::Vector3d(-2.6, 0.0, 0.0));
    const int c2 = add_atom(mol, 6, Eigen::Vector3d(0.0, 0.0, 0.0));
    const int c3 = add_atom(mol, 6, Eigen::Vector3d(2.6, 0.0, 0.0));
    const int o2 = add_atom(mol, 8, Eigen::Vector3d(5.0, 0.8, 0.0));
    const int m1 = add_atom(mol, 6, Eigen::Vector3d(-4.0, -2.1, 0.0));
    const int m2 = add_atom(mol, 6, Eigen::Vector3d(4.0, -2.1, 0.0));
    mol.add_bond(o1, c1, BondOrder::Single);
    mol.add_bond(c1, c2, BondOrder::Single);
    mol.add_bond(c2, c3, BondOrder::Double);
    mol.add_bond(c3, o2, BondOrder::Single);
    mol.add_bond(c1, m1, BondOrder::Single);
    mol.add_bond(c3, m2, BondOrder::Single);
    e.donor_atoms = {0, 4};
    return e;
}

LigandEntry make_oxalate() {
    LigandEntry e;
    e.name = "Oxalate";
    e.abbreviation = "ox";
    e.formula = "C₂O₄²⁻";
    e.category = "Polydentate";
    e.denticity = LigandDenticity::Bidentate;
    e.charge = -2;
    e.field_strength = "weak";
    e.typical_delta_oct_eV = 1.1;
    auto& mol = e.molecule;
    const int o1 = add_atom(mol, 8, Eigen::Vector3d(-4.4, 1.4, 0.0), -1);
    const int c1 = add_atom(mol, 6, Eigen::Vector3d(-2.3, 0.5, 0.0));
    const int o2 = add_atom(mol, 8, Eigen::Vector3d(-2.2, -1.9, 0.0));
    const int c2 = add_atom(mol, 6, Eigen::Vector3d(2.3, 0.5, 0.0));
    const int o3 = add_atom(mol, 8, Eigen::Vector3d(2.2, -1.9, 0.0));
    const int o4 = add_atom(mol, 8, Eigen::Vector3d(4.4, 1.4, 0.0), -1);
    mol.add_bond(o1, c1, BondOrder::Single);
    mol.add_bond(c1, o2, BondOrder::Double);
    mol.add_bond(c1, c2, BondOrder::Single);
    mol.add_bond(c2, o3, BondOrder::Double);
    mol.add_bond(c2, o4, BondOrder::Single);
    e.donor_atoms = {0, 5};
    return e;
}

LigandEntry make_phen() {
    LigandEntry e;
    e.name = "Phenanthroline";
    e.abbreviation = "phen";
    e.formula = "C₁₂H₈N₂";
    e.category = "Polydentate";
    e.denticity = LigandDenticity::Bidentate;
    e.charge = 0;
    e.field_strength = "strong";
    e.typical_delta_oct_eV = 1.65;
    auto& mol = e.molecule;
    const std::vector<std::pair<int, Eigen::Vector3d>> atoms = {
        {7, {-5.2, 1.0, 0.0}}, {6, {-3.0, 2.1, 0.0}}, {6, {-0.6, 1.4, 0.0}}, {6, {1.8, 2.1, 0.0}},
        {7, {4.0, 1.0, 0.0}}, {6, {-3.0, -1.1, 0.0}}, {6, {-0.6, -1.8, 0.0}}, {6, {1.8, -1.1, 0.0}},
        {6, {-0.6, 3.9, 0.0}}, {6, {1.8, 4.6, 0.0}}, {6, {-0.6, -4.3, 0.0}}, {6, {1.8, -3.6, 0.0}}
    };
    for (const auto& [z, p] : atoms) add_atom(mol, z, p);
    const std::vector<std::pair<int, int>> bonds = {
        {0,1},{1,2},{2,3},{3,4},{0,5},{5,6},{6,7},{7,4},{2,6},{2,8},{8,9},{9,3},{6,10},{10,11},{11,7}
    };
    for (auto [a,b] : bonds) mol.add_bond(a,b,BondOrder::Aromatic);
    e.donor_atoms = {0, 4};
    return e;
}

LigandEntry make_dien() {
    LigandEntry e;
    e.name = "Diethylenetriamine";
    e.abbreviation = "dien";
    e.formula = "dien";
    e.category = "Polydentate";
    e.denticity = LigandDenticity::Tridentate;
    e.charge = 0;
    e.field_strength = "strong";
    e.typical_delta_oct_eV = 1.58;
    auto& mol = e.molecule;
    const std::vector<std::pair<int, Eigen::Vector3d>> atoms = {
        {7, {-6.0, 0.0, 0.0}}, {6, {-3.5, 0.8, 0.8}}, {6, {-1.0, 0.0, -0.8}},
        {7, {1.5, 0.8, 0.0}}, {6, {4.0, 0.0, 0.8}}, {6, {6.5, -0.8, -0.8}}, {7, {9.0, 0.0, 0.0}}
    };
    for (const auto& [z,p] : atoms) add_atom(mol,z,p);
    for (int i = 0; i < 6; ++i) mol.add_bond(i, i+1, BondOrder::Single);
    e.donor_atoms = {0, 3, 6};
    return e;
}

LigandEntry make_edta() {
    LigandEntry e;
    e.name = "EDTA";
    e.abbreviation = "edta";
    e.formula = "EDTA⁴⁻";
    e.category = "Polydentate";
    e.denticity = LigandDenticity::Hexadentate;
    e.charge = -4;
    e.field_strength = "strong";
    e.typical_delta_oct_eV = 1.8;
    auto& mol = e.molecule;
    const std::vector<std::pair<int, Eigen::Vector3d>> donor_atoms = {
        {7, {-5.0, 1.5, 0.0}}, {7, {5.0, 1.5, 0.0}},
        {8, {-3.8, -3.2, 0.0}}, {8, {-1.2, -4.0, 0.0}},
        {8, {1.2, -4.0, 0.0}}, {8, {3.8, -3.2, 0.0}}
    };
    for (const auto& [z,p] : donor_atoms) add_atom(mol,z,p, z == 8 ? -1 : 0);
    add_atom(mol, 6, Eigen::Vector3d(-2.5, 0.0, 0.0));
    add_atom(mol, 6, Eigen::Vector3d(0.0, 1.2, 0.0));
    add_atom(mol, 6, Eigen::Vector3d(2.5, 0.0, 0.0));
    for (int i = 0; i < 6; ++i) {
        const int carbon = 6 + std::min(i / 2, 2);
        mol.add_bond(i, carbon, BondOrder::Single);
    }
    mol.add_bond(6, 7, BondOrder::Single);
    mol.add_bond(7, 8, BondOrder::Single);
    e.donor_atoms = {0, 1, 2, 3, 4, 5};
    return e;
}

}  // namespace

LigandLibrary::LigandLibrary() {
    build_library();
}

const std::vector<LigandEntry>& LigandLibrary::all() const {
    return ligands_;
}

std::vector<const LigandEntry*> LigandLibrary::by_category(const std::string& category) const {
    std::vector<const LigandEntry*> out;
    for (const auto& ligand : ligands_) {
        if (ligand.category == category) {
            out.push_back(&ligand);
        }
    }
    return out;
}

std::vector<const LigandEntry*> LigandLibrary::by_denticity(LigandDenticity d) const {
    std::vector<const LigandEntry*> out;
    for (const auto& ligand : ligands_) {
        if (ligand.denticity == d) {
            out.push_back(&ligand);
        }
    }
    return out;
}

const LigandEntry* LigandLibrary::find(const std::string& name_or_abbrev) const {
    const std::string key = lower_copy(name_or_abbrev);
    for (const auto& ligand : ligands_) {
        if (lower_copy(ligand.name) == key || lower_copy(ligand.abbreviation) == key) {
            return &ligand;
        }
    }
    return nullptr;
}

std::vector<std::string> LigandLibrary::categories() const {
    std::vector<std::string> out;
    for (const auto& ligand : ligands_) {
        if (std::find(out.begin(), out.end(), ligand.category) == out.end()) {
            out.push_back(ligand.category);
        }
    }
    return out;
}

LigandSpec LigandLibrary::to_ligand_spec(const LigandEntry& entry, int donor_index) const {
    if (donor_index < 0 || donor_index >= static_cast<int>(entry.donor_atoms.size())) {
        throw std::out_of_range("Invalid donor index");
    }
    const int donor_atom = entry.donor_atoms[static_cast<std::size_t>(donor_index)];
    if (donor_atom < 0 || donor_atom >= entry.molecule.num_atoms()) {
        throw std::out_of_range("Invalid donor atom");
    }

    LigandSpec spec;
    spec.name = entry.abbreviation;
    spec.donor_Z = entry.molecule.atom(donor_atom).Z;

    MolecularSystem reordered;
    reordered.set_name(entry.name);
    reordered.set_charge(entry.charge);
    reordered.set_multiplicity(entry.molecule.multiplicity());

    std::vector<int> old_to_new(entry.molecule.num_atoms(), -1);
    old_to_new[donor_atom] = reordered.add_atom(entry.molecule.atom(donor_atom));
    for (int i = 0; i < entry.molecule.num_atoms(); ++i) {
        if (i == donor_atom) continue;
        old_to_new[i] = reordered.add_atom(entry.molecule.atom(i));
    }
    for (const Bond& bond : entry.molecule.bonds()) {
        reordered.add_bond(old_to_new[bond.atom_i], old_to_new[bond.atom_j], bond.order);
    }
    spec.ligand = reordered;
    return spec;
}

void LigandLibrary::build_library() {
    ligands_.clear();
    ligands_.push_back(make_water());
    ligands_.push_back(make_ammonia());
    ligands_.push_back(make_carbon_monoxide());
    ligands_.push_back(make_phosphine());
    ligands_.push_back(make_pyridine());
    ligands_.push_back(make_acetonitrile());
    ligands_.push_back(make_halide("Chloride", "cl", "Cl⁻", 17, "weak", 1.1));
    ligands_.push_back(make_halide("Bromide", "br", "Br⁻", 35, "very weak", 1.0));
    ligands_.push_back(make_halide("Iodide", "i", "I⁻", 53, "very weak", 0.9));
    ligands_.push_back(make_halide("Fluoride", "f", "F⁻", 9, "weak", 1.2));
    ligands_.push_back(make_hydroxide());
    ligands_.push_back(make_cyanide());
    ligands_.push_back(make_thiocyanate());
    ligands_.push_back(make_nitrite());
    ligands_.push_back(make_azide());
    ligands_.push_back(make_en());
    ligands_.push_back(make_bipy());
    ligands_.push_back(make_acac());
    ligands_.push_back(make_oxalate());
    ligands_.push_back(make_phen());
    ligands_.push_back(make_dien());
    ligands_.push_back(make_edta());
}

}  // namespace sbox::chem
