#pragma once

#include "core/molecular_system.h"

#include <Eigen/Core>

#include <memory>
#include <string>
#include <vector>

namespace sbox::editor {

class Command {
public:
    virtual ~Command() = default;
    virtual void execute(sbox::chem::MolecularSystem& mol) = 0;
    virtual void undo(sbox::chem::MolecularSystem& mol) = 0;
    virtual std::string description() const = 0;
};

class CommandStack {
public:
    void execute(std::unique_ptr<Command> cmd, sbox::chem::MolecularSystem& mol);

    void undo(sbox::chem::MolecularSystem& mol);
    void redo(sbox::chem::MolecularSystem& mol);

    bool can_undo() const;
    bool can_redo() const;

    std::string undo_description() const;
    std::string redo_description() const;

    int size() const;
    void clear();

private:
    std::vector<std::unique_ptr<Command>> history_;
    int current_ = -1;
};

class AddAtomCommand : public Command {
public:
    explicit AddAtomCommand(sbox::chem::Atom atom);
    void execute(sbox::chem::MolecularSystem& mol) override;
    void undo(sbox::chem::MolecularSystem& mol) override;
    std::string description() const override;
private:
    sbox::chem::Atom atom_;
    int added_index_ = -1;
};

class RemoveAtomCommand : public Command {
public:
    explicit RemoveAtomCommand(int atom_index);
    void execute(sbox::chem::MolecularSystem& mol) override;
    void undo(sbox::chem::MolecularSystem& mol) override;
    std::string description() const override;
private:
    int index_;
    std::vector<sbox::chem::Atom> saved_atoms_;
    std::vector<sbox::chem::Bond> saved_bonds_;
    std::string saved_name_;
    int saved_charge_ = 0;
    int saved_multiplicity_ = 1;
};

class MoveAtomCommand : public Command {
public:
    MoveAtomCommand(int atom_index, Eigen::Vector3d new_position);
    void execute(sbox::chem::MolecularSystem& mol) override;
    void undo(sbox::chem::MolecularSystem& mol) override;
    std::string description() const override;
private:
    int index_;
    Eigen::Vector3d new_pos_;
    Eigen::Vector3d old_pos_ = Eigen::Vector3d::Zero();
};

class MoveAtomsCommand : public Command {
public:
    MoveAtomsCommand(std::vector<int> indices, std::vector<Eigen::Vector3d> new_positions);
    void execute(sbox::chem::MolecularSystem& mol) override;
    void undo(sbox::chem::MolecularSystem& mol) override;
    std::string description() const override;
private:
    std::vector<int> indices_;
    std::vector<Eigen::Vector3d> new_positions_;
    std::vector<Eigen::Vector3d> old_positions_;
};

class AddBondCommand : public Command {
public:
    AddBondCommand(int atom_i, int atom_j, sbox::chem::BondOrder order);
    void execute(sbox::chem::MolecularSystem& mol) override;
    void undo(sbox::chem::MolecularSystem& mol) override;
    std::string description() const override;
private:
    int atom_i_;
    int atom_j_;
    sbox::chem::BondOrder order_;
    int added_bond_index_ = -1;
};

class RemoveBondCommand : public Command {
public:
    explicit RemoveBondCommand(int bond_index);
    void execute(sbox::chem::MolecularSystem& mol) override;
    void undo(sbox::chem::MolecularSystem& mol) override;
    std::string description() const override;
private:
    int bond_index_;
    sbox::chem::Bond saved_bond_{0, 0, sbox::chem::BondOrder::Unknown};
};

class ChangeBondOrderCommand : public Command {
public:
    ChangeBondOrderCommand(int bond_index, sbox::chem::BondOrder new_order);
    void execute(sbox::chem::MolecularSystem& mol) override;
    void undo(sbox::chem::MolecularSystem& mol) override;
    std::string description() const override;
private:
    int bond_index_;
    sbox::chem::BondOrder new_order_;
    sbox::chem::BondOrder old_order_ = sbox::chem::BondOrder::Unknown;
    int atom_i_ = -1;
    int atom_j_ = -1;
};

class SetChargeCommand : public Command {
public:
    SetChargeCommand(int atom_index, int new_formal_charge);
    void execute(sbox::chem::MolecularSystem& mol) override;
    void undo(sbox::chem::MolecularSystem& mol) override;
    std::string description() const override;
private:
    int index_;
    int new_charge_;
    int old_charge_ = 0;
};

class SetElementCommand : public Command {
public:
    SetElementCommand(int atom_index, int new_Z);
    void execute(sbox::chem::MolecularSystem& mol) override;
    void undo(sbox::chem::MolecularSystem& mol) override;
    std::string description() const override;
private:
    int index_;
    int new_Z_;
    int old_Z_ = 0;
};

class AddFragmentCommand : public Command {
public:
    AddFragmentCommand(const sbox::chem::MolecularSystem& fragment, const std::string& fragment_name);
    void execute(sbox::chem::MolecularSystem& mol) override;
    void undo(sbox::chem::MolecularSystem& mol) override;
    std::string description() const override;
private:
    sbox::chem::MolecularSystem fragment_;
    std::string name_;
    int first_atom_index_ = -1;
    int num_atoms_added_ = 0;
    int num_bonds_added_ = 0;
    std::vector<sbox::chem::Atom> saved_atoms_;
    std::vector<sbox::chem::Bond> saved_bonds_;
    std::string saved_name_;
    int saved_charge_ = 0;
    int saved_multiplicity_ = 1;
};

class AddHydrogensCommand : public Command {
public:
    explicit AddHydrogensCommand(int atom_index = -1);
    void execute(sbox::chem::MolecularSystem& mol) override;
    void undo(sbox::chem::MolecularSystem& mol) override;
    std::string description() const override;
private:
    int atom_index_;
    std::vector<sbox::chem::Atom> saved_atoms_;
    std::vector<sbox::chem::Bond> saved_bonds_;
    std::string saved_name_;
    int saved_charge_ = 0;
    int saved_multiplicity_ = 1;
};

class RemoveHydrogensCommand : public Command {
public:
    void execute(sbox::chem::MolecularSystem& mol) override;
    void undo(sbox::chem::MolecularSystem& mol) override;
    std::string description() const override;
private:
    std::vector<sbox::chem::Atom> saved_atoms_;
    std::vector<sbox::chem::Bond> saved_bonds_;
    std::string saved_name_;
    int saved_charge_ = 0;
    int saved_multiplicity_ = 1;
};

}  // namespace sbox::editor
