#pragma once

#include <Eigen/Core>

#include <cstdint>
#include <string>
#include <vector>

namespace sbox::chem {

struct Atom {
    int Z;
    Eigen::Vector3d position;
    std::string label;
    int formal_charge = 0;
};

enum class BondOrder : uint8_t {
    Unknown = 0,
    Single = 1,
    Double = 2,
    Triple = 3,
    Aromatic = 4,
};

struct Bond {
    int atom_i;
    int atom_j;
    BondOrder order;
};

class MolecularSystem {
public:
    int num_atoms() const;
    int num_bonds() const;

    const Atom& atom(int i) const;
    Atom& atom(int i);
    const Bond& bond(int i) const;

    int add_atom(Atom a);
    void remove_atom(int i);
    int add_bond(int i, int j, BondOrder order);
    void remove_bond(int i);
    bool has_bond(int i, int j) const;
    void perceive_bonds(double tolerance = 1.15);

    std::vector<int> neighbors(int atom_i) const;
    int coordination_number(int atom_i) const;

    double distance(int i, int j) const;
    double angle(int i, int j, int k) const;
    double dihedral(int i, int j, int k, int l) const;

    Eigen::Vector3d center_of_mass() const;
    void center();
    void clear();

    const std::string& name() const;
    void set_name(const std::string& name);
    int charge() const;
    void set_charge(int charge);
    int multiplicity() const;
    void set_multiplicity(int multiplicity);
    const std::vector<Atom>& atoms() const;
    const std::vector<Bond>& bonds() const;

private:
    std::vector<Atom> atoms_;
    std::vector<Bond> bonds_;
    std::string name_;
    int charge_ = 0;
    int multiplicity_ = 1;
};

}  // namespace sbox::chem
