#include "editor/command.h"

#include "core/elements.h"
#include "core/valence.h"

#include <stdexcept>
#include <utility>

namespace sbox::editor {

namespace {

void restore_snapshot(sbox::chem::MolecularSystem& mol,
                      const std::vector<sbox::chem::Atom>& atoms,
                      const std::vector<sbox::chem::Bond>& bonds,
                      const std::string& name,
                      int charge,
                      int multiplicity) {
    mol.clear();
    mol.set_name(name);
    mol.set_charge(charge);
    mol.set_multiplicity(multiplicity);
    for (const auto& atom : atoms) {
        mol.add_atom(atom);
    }
    for (const auto& bond : bonds) {
        mol.add_bond(bond.atom_i, bond.atom_j, bond.order);
    }
}

const char* bond_order_name(sbox::chem::BondOrder order) {
    switch (order) {
    case sbox::chem::BondOrder::Unknown: return "Unknown";
    case sbox::chem::BondOrder::Single: return "Single";
    case sbox::chem::BondOrder::Double: return "Double";
    case sbox::chem::BondOrder::Triple: return "Triple";
    case sbox::chem::BondOrder::Aromatic: return "Aromatic";
    }
    return "Unknown";
}

}  // namespace

void CommandStack::execute(std::unique_ptr<Command> cmd, sbox::chem::MolecularSystem& mol) {
    if (!cmd) {
        throw std::invalid_argument("Cannot execute null command");
    }
    if (current_ + 1 < static_cast<int>(history_.size())) {
        history_.erase(history_.begin() + current_ + 1, history_.end());
    }
    cmd->execute(mol);
    history_.push_back(std::move(cmd));
    current_ = static_cast<int>(history_.size()) - 1;
}

void CommandStack::undo(sbox::chem::MolecularSystem& mol) {
    if (!can_undo()) {
        return;
    }
    history_[static_cast<std::size_t>(current_)]->undo(mol);
    --current_;
}

void CommandStack::redo(sbox::chem::MolecularSystem& mol) {
    if (!can_redo()) {
        return;
    }
    ++current_;
    history_[static_cast<std::size_t>(current_)]->execute(mol);
}

bool CommandStack::can_undo() const {
    return current_ >= 0 && current_ < static_cast<int>(history_.size());
}

bool CommandStack::can_redo() const {
    return current_ + 1 >= 0 && current_ + 1 < static_cast<int>(history_.size());
}

std::string CommandStack::undo_description() const {
    return can_undo() ? history_[static_cast<std::size_t>(current_)]->description() : "";
}

std::string CommandStack::redo_description() const {
    return can_redo() ? history_[static_cast<std::size_t>(current_ + 1)]->description() : "";
}

int CommandStack::size() const {
    return static_cast<int>(history_.size());
}

void CommandStack::clear() {
    history_.clear();
    current_ = -1;
}

AddAtomCommand::AddAtomCommand(sbox::chem::Atom atom)
    : atom_(std::move(atom)) {}

void AddAtomCommand::execute(sbox::chem::MolecularSystem& mol) {
    added_index_ = mol.add_atom(atom_);
}

void AddAtomCommand::undo(sbox::chem::MolecularSystem& mol) {
    if (added_index_ >= 0) {
        mol.remove_atom(added_index_);
    }
}

std::string AddAtomCommand::description() const {
    return std::string("Add ") + sbox::elements::get_element(atom_.Z).symbol;
}

RemoveAtomCommand::RemoveAtomCommand(int atom_index)
    : index_(atom_index) {}

void RemoveAtomCommand::execute(sbox::chem::MolecularSystem& mol) {
    saved_atoms_ = mol.atoms();
    saved_bonds_ = mol.bonds();
    saved_name_ = mol.name();
    saved_charge_ = mol.charge();
    saved_multiplicity_ = mol.multiplicity();
    mol.remove_atom(index_);
}

void RemoveAtomCommand::undo(sbox::chem::MolecularSystem& mol) {
    restore_snapshot(mol, saved_atoms_, saved_bonds_, saved_name_, saved_charge_, saved_multiplicity_);
}

std::string RemoveAtomCommand::description() const {
    return "Remove atom " + std::to_string(index_ + 1);
}

MoveAtomCommand::MoveAtomCommand(int atom_index, Eigen::Vector3d new_position)
    : index_(atom_index), new_pos_(std::move(new_position)) {}

void MoveAtomCommand::execute(sbox::chem::MolecularSystem& mol) {
    old_pos_ = mol.atom(index_).position;
    mol.atom(index_).position = new_pos_;
}

void MoveAtomCommand::undo(sbox::chem::MolecularSystem& mol) {
    mol.atom(index_).position = old_pos_;
}

std::string MoveAtomCommand::description() const {
    return "Move atom " + std::to_string(index_ + 1);
}

MoveAtomsCommand::MoveAtomsCommand(std::vector<int> indices, std::vector<Eigen::Vector3d> new_positions)
    : indices_(std::move(indices)), new_positions_(std::move(new_positions)) {
    if (indices_.size() != new_positions_.size()) {
        throw std::invalid_argument("MoveAtomsCommand requires matching index and position counts");
    }
}

void MoveAtomsCommand::execute(sbox::chem::MolecularSystem& mol) {
    old_positions_.clear();
    old_positions_.reserve(indices_.size());
    for (std::size_t i = 0; i < indices_.size(); ++i) {
        old_positions_.push_back(mol.atom(indices_[i]).position);
        mol.atom(indices_[i]).position = new_positions_[i];
    }
}

void MoveAtomsCommand::undo(sbox::chem::MolecularSystem& mol) {
    for (std::size_t i = 0; i < indices_.size(); ++i) {
        mol.atom(indices_[i]).position = old_positions_[i];
    }
}

std::string MoveAtomsCommand::description() const {
    return "Move " + std::to_string(indices_.size()) + " atoms";
}

AddBondCommand::AddBondCommand(int atom_i, int atom_j, sbox::chem::BondOrder order)
    : atom_i_(atom_i), atom_j_(atom_j), order_(order) {}

void AddBondCommand::execute(sbox::chem::MolecularSystem& mol) {
    added_bond_index_ = mol.add_bond(atom_i_, atom_j_, order_);
}

void AddBondCommand::undo(sbox::chem::MolecularSystem& mol) {
    if (added_bond_index_ >= 0) {
        mol.remove_bond(added_bond_index_);
    }
}

std::string AddBondCommand::description() const {
    return "Add bond " + std::to_string(atom_i_ + 1) + "-" + std::to_string(atom_j_ + 1);
}

RemoveBondCommand::RemoveBondCommand(int bond_index)
    : bond_index_(bond_index) {}

void RemoveBondCommand::execute(sbox::chem::MolecularSystem& mol) {
    saved_bond_ = mol.bond(bond_index_);
    mol.remove_bond(bond_index_);
}

void RemoveBondCommand::undo(sbox::chem::MolecularSystem& mol) {
    mol.add_bond(saved_bond_.atom_i, saved_bond_.atom_j, saved_bond_.order);
}

std::string RemoveBondCommand::description() const {
    return "Remove bond " + std::to_string(bond_index_ + 1);
}

ChangeBondOrderCommand::ChangeBondOrderCommand(int bond_index, sbox::chem::BondOrder new_order)
    : bond_index_(bond_index), new_order_(new_order) {}

void ChangeBondOrderCommand::execute(sbox::chem::MolecularSystem& mol) {
    const sbox::chem::Bond bond = mol.bond(bond_index_);
    atom_i_ = bond.atom_i;
    atom_j_ = bond.atom_j;
    old_order_ = bond.order;
    mol.remove_bond(bond_index_);
    bond_index_ = mol.add_bond(atom_i_, atom_j_, new_order_);
}

void ChangeBondOrderCommand::undo(sbox::chem::MolecularSystem& mol) {
    mol.remove_bond(bond_index_);
    bond_index_ = mol.add_bond(atom_i_, atom_j_, old_order_);
}

std::string ChangeBondOrderCommand::description() const {
    return "Change bond " + std::to_string(atom_i_ + 1) + "-" + std::to_string(atom_j_ + 1) +
           " to " + bond_order_name(new_order_);
}

SetChargeCommand::SetChargeCommand(int atom_index, int new_formal_charge)
    : index_(atom_index), new_charge_(new_formal_charge) {}

void SetChargeCommand::execute(sbox::chem::MolecularSystem& mol) {
    old_charge_ = mol.atom(index_).formal_charge;
    mol.atom(index_).formal_charge = new_charge_;
}

void SetChargeCommand::undo(sbox::chem::MolecularSystem& mol) {
    mol.atom(index_).formal_charge = old_charge_;
}

std::string SetChargeCommand::description() const {
    return "Set charge on atom " + std::to_string(index_ + 1);
}

SetElementCommand::SetElementCommand(int atom_index, int new_Z)
    : index_(atom_index), new_Z_(new_Z) {}

void SetElementCommand::execute(sbox::chem::MolecularSystem& mol) {
    old_Z_ = mol.atom(index_).Z;
    mol.atom(index_).Z = new_Z_;
}

void SetElementCommand::undo(sbox::chem::MolecularSystem& mol) {
    mol.atom(index_).Z = old_Z_;
}

std::string SetElementCommand::description() const {
    return "Change atom " + std::to_string(index_ + 1) + " to " + sbox::elements::get_element(new_Z_).symbol;
}

AddFragmentCommand::AddFragmentCommand(const sbox::chem::MolecularSystem& fragment, const std::string& fragment_name)
    : fragment_(fragment), name_(fragment_name) {}

void AddFragmentCommand::execute(sbox::chem::MolecularSystem& mol) {
    saved_atoms_ = mol.atoms();
    saved_bonds_ = mol.bonds();
    saved_name_ = mol.name();
    saved_charge_ = mol.charge();
    saved_multiplicity_ = mol.multiplicity();

    first_atom_index_ = mol.num_atoms();
    num_atoms_added_ = fragment_.num_atoms();
    num_bonds_added_ = fragment_.num_bonds();

    for (const auto& atom : fragment_.atoms()) {
        mol.add_atom(atom);
    }
    for (const auto& bond : fragment_.bonds()) {
        mol.add_bond(first_atom_index_ + bond.atom_i, first_atom_index_ + bond.atom_j, bond.order);
    }
}

void AddFragmentCommand::undo(sbox::chem::MolecularSystem& mol) {
    restore_snapshot(mol, saved_atoms_, saved_bonds_, saved_name_, saved_charge_, saved_multiplicity_);
}

std::string AddFragmentCommand::description() const {
    return "Add " + name_;
}

AddHydrogensCommand::AddHydrogensCommand(int atom_index)
    : atom_index_(atom_index) {}

void AddHydrogensCommand::execute(sbox::chem::MolecularSystem& mol) {
    saved_atoms_ = mol.atoms();
    saved_bonds_ = mol.bonds();
    saved_name_ = mol.name();
    saved_charge_ = mol.charge();
    saved_multiplicity_ = mol.multiplicity();
    if (atom_index_ >= 0) {
        sbox::chem::add_hydrogens(mol, atom_index_);
    } else {
        sbox::chem::add_hydrogens(mol);
    }
}

void AddHydrogensCommand::undo(sbox::chem::MolecularSystem& mol) {
    restore_snapshot(mol, saved_atoms_, saved_bonds_, saved_name_, saved_charge_, saved_multiplicity_);
}

std::string AddHydrogensCommand::description() const {
    return atom_index_ >= 0 ? "Add hydrogens to atom " + std::to_string(atom_index_ + 1) : "Add hydrogens";
}

void RemoveHydrogensCommand::execute(sbox::chem::MolecularSystem& mol) {
    saved_atoms_ = mol.atoms();
    saved_bonds_ = mol.bonds();
    saved_name_ = mol.name();
    saved_charge_ = mol.charge();
    saved_multiplicity_ = mol.multiplicity();
    sbox::chem::remove_hydrogens(mol);
}

void RemoveHydrogensCommand::undo(sbox::chem::MolecularSystem& mol) {
    restore_snapshot(mol, saved_atoms_, saved_bonds_, saved_name_, saved_charge_, saved_multiplicity_);
}

std::string RemoveHydrogensCommand::description() const {
    return "Remove hydrogens";
}

}  // namespace sbox::editor
