#include "editor/command.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace {

sbox::chem::MolecularSystem make_chain3() {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "", 0});
    mol.add_atom({6, Eigen::Vector3d(1.0, 0.0, 0.0), "", 0});
    mol.add_atom({6, Eigen::Vector3d(2.0, 0.0, 0.0), "", 0});
    mol.add_bond(0, 1, sbox::chem::BondOrder::Single);
    mol.add_bond(1, 2, sbox::chem::BondOrder::Single);
    return mol;
}

}  // namespace

TEST(CommandStackTest, AddAtomCommandUndoRedo) {
    sbox::chem::MolecularSystem mol;
    sbox::editor::CommandStack stack;

    stack.execute(std::make_unique<sbox::editor::AddAtomCommand>(
                      sbox::chem::Atom{1, Eigen::Vector3d::Zero(), "", 0}),
                  mol);
    EXPECT_EQ(mol.num_atoms(), 1);

    stack.undo(mol);
    EXPECT_EQ(mol.num_atoms(), 0);

    stack.redo(mol);
    EXPECT_EQ(mol.num_atoms(), 1);
}

TEST(CommandStackTest, RemoveAtomCommandRestoresSnapshot) {
    sbox::chem::MolecularSystem mol = make_chain3();
    sbox::editor::CommandStack stack;

    stack.execute(std::make_unique<sbox::editor::RemoveAtomCommand>(1), mol);
    EXPECT_EQ(mol.num_atoms(), 2);
    EXPECT_EQ(mol.num_bonds(), 0);

    stack.undo(mol);
    EXPECT_EQ(mol.num_atoms(), 3);
    EXPECT_EQ(mol.num_bonds(), 2);
}

TEST(CommandStackTest, MoveAtomCommandUndoRestoresPosition) {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({6, Eigen::Vector3d::Zero(), "", 0});
    sbox::editor::CommandStack stack;

    stack.execute(std::make_unique<sbox::editor::MoveAtomCommand>(0, Eigen::Vector3d(1.0, 2.0, 3.0)), mol);
    EXPECT_EQ(mol.atom(0).position, Eigen::Vector3d(1.0, 2.0, 3.0));

    stack.undo(mol);
    EXPECT_EQ(mol.atom(0).position, Eigen::Vector3d::Zero());
}

TEST(CommandStackTest, AddBondAndRemoveBondCommands) {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({6, Eigen::Vector3d::Zero(), "", 0});
    mol.add_atom({6, Eigen::Vector3d(1.0, 0.0, 0.0), "", 0});

    sbox::editor::CommandStack stack;
    stack.execute(std::make_unique<sbox::editor::AddBondCommand>(0, 1, sbox::chem::BondOrder::Single), mol);
    EXPECT_TRUE(mol.has_bond(0, 1));

    stack.execute(std::make_unique<sbox::editor::RemoveBondCommand>(0), mol);
    EXPECT_FALSE(mol.has_bond(0, 1));

    stack.undo(mol);
    EXPECT_TRUE(mol.has_bond(0, 1));
}

TEST(CommandStackTest, ChangeBondOrderCommandUndo) {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({6, Eigen::Vector3d::Zero(), "", 0});
    mol.add_atom({6, Eigen::Vector3d(1.0, 0.0, 0.0), "", 0});
    mol.add_bond(0, 1, sbox::chem::BondOrder::Single);

    sbox::editor::CommandStack stack;
    stack.execute(std::make_unique<sbox::editor::ChangeBondOrderCommand>(0, sbox::chem::BondOrder::Double), mol);
    EXPECT_EQ(mol.bond(0).order, sbox::chem::BondOrder::Double);

    stack.undo(mol);
    EXPECT_EQ(mol.bond(0).order, sbox::chem::BondOrder::Single);
}

TEST(CommandStackTest, StackUndoRedoAndClearsRedoHistory) {
    sbox::chem::MolecularSystem mol;
    sbox::editor::CommandStack stack;

    stack.execute(std::make_unique<sbox::editor::AddAtomCommand>(
                      sbox::chem::Atom{1, Eigen::Vector3d::Zero(), "", 0}),
                  mol);
    stack.execute(std::make_unique<sbox::editor::AddAtomCommand>(
                      sbox::chem::Atom{6, Eigen::Vector3d(1.0, 0.0, 0.0), "", 0}),
                  mol);
    stack.execute(std::make_unique<sbox::editor::MoveAtomCommand>(1, Eigen::Vector3d(2.0, 0.0, 0.0)), mol);

    EXPECT_TRUE(stack.can_undo());
    EXPECT_FALSE(stack.can_redo());
    EXPECT_EQ(stack.size(), 3);

    stack.undo(mol);
    stack.undo(mol);
    EXPECT_TRUE(stack.can_redo());
    EXPECT_EQ(mol.num_atoms(), 1);

    stack.redo(mol);
    EXPECT_EQ(mol.num_atoms(), 2);
    EXPECT_TRUE(stack.can_undo());
    EXPECT_TRUE(stack.can_redo());

    stack.execute(std::make_unique<sbox::editor::SetChargeCommand>(0, 1), mol);
    EXPECT_FALSE(stack.can_redo());
}

TEST(CommandStackTest, AddFragmentCommandUndoRedo) {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "", 0});
    mol.add_atom({6, Eigen::Vector3d(1.0, 0.0, 0.0), "", 0});
    mol.add_bond(0, 1, sbox::chem::BondOrder::Single);

    sbox::chem::MolecularSystem fragment;
    fragment.add_atom({1, Eigen::Vector3d(0.0, 1.0, 0.0), "", 0});
    fragment.add_atom({1, Eigen::Vector3d(1.0, 1.0, 0.0), "", 0});
    fragment.add_atom({8, Eigen::Vector3d(0.5, 2.0, 0.0), "", 0});
    fragment.add_bond(0, 2, sbox::chem::BondOrder::Single);
    fragment.add_bond(1, 2, sbox::chem::BondOrder::Single);

    sbox::editor::CommandStack stack;
    stack.execute(std::make_unique<sbox::editor::AddFragmentCommand>(fragment, "fragment"), mol);
    EXPECT_EQ(mol.num_atoms(), 5);

    stack.undo(mol);
    EXPECT_EQ(mol.num_atoms(), 2);

    stack.redo(mol);
    EXPECT_EQ(mol.num_atoms(), 5);
}
