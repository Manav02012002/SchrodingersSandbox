#include "editor/picking.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include <cmath>

namespace {

Eigen::Matrix4f perspective(float fovy_rad, float aspect, float z_near, float z_far) {
    const float f = 1.0f / std::tan(fovy_rad * 0.5f);
    Eigen::Matrix4f m = Eigen::Matrix4f::Zero();
    m(0, 0) = f / aspect;
    m(1, 1) = f;
    m(2, 2) = (z_far + z_near) / (z_near - z_far);
    m(2, 3) = (2.0f * z_far * z_near) / (z_near - z_far);
    m(3, 2) = -1.0f;
    return m;
}

}  // namespace

TEST(PickingTest, ScreenToRayCenterPointsAlongMinusZ) {
    const Eigen::Matrix4f proj = perspective(static_cast<float>(std::acos(-1.0) / 3.0), 1.0f, 0.1f, 100.0f);
    const Eigen::Matrix4f inv_vp = proj.inverse();

    const sbox::editor::Ray ray = sbox::editor::screen_to_ray(400.0f, 400.0f, ImVec2(0.0f, 0.0f), ImVec2(800.0f, 800.0f), inv_vp);
    EXPECT_NEAR(ray.direction.x(), 0.0f, 1e-4f);
    EXPECT_NEAR(ray.direction.y(), 0.0f, 1e-4f);
    EXPECT_LT(ray.direction.z(), -0.9f);
}

TEST(PickingTest, PickAtomHitsClosestSphere) {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "", 0});
    mol.add_atom({6, Eigen::Vector3d(3.0, 0.0, 0.0), "", 0});
    mol.add_atom({6, Eigen::Vector3d(0.0, 3.0, 0.0), "", 0});

    sbox::editor::Ray ray0{Eigen::Vector3f(0.0f, 0.0f, 10.0f), Eigen::Vector3f(0.0f, 0.0f, -1.0f)};
    sbox::editor::PickResult hit0 = sbox::editor::pick_atom(ray0, mol);
    EXPECT_EQ(hit0.type, sbox::editor::PickResult::Type::Atom);
    EXPECT_EQ(hit0.index, 0);

    const Eigen::Vector3f dir1 = (Eigen::Vector3f(3.0f, 0.0f, 0.0f) - Eigen::Vector3f(0.0f, 0.0f, 10.0f)).normalized();
    sbox::editor::Ray ray1{Eigen::Vector3f(0.0f, 0.0f, 10.0f), dir1};
    sbox::editor::PickResult hit1 = sbox::editor::pick_atom(ray1, mol);
    EXPECT_EQ(hit1.type, sbox::editor::PickResult::Type::Atom);
    EXPECT_EQ(hit1.index, 1);
}

TEST(PickingTest, PickAtomMissReturnsNone) {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "", 0});
    sbox::editor::Ray ray{Eigen::Vector3f(0.0f, 0.0f, 10.0f), Eigen::Vector3f(0.0f, 1.0f, 0.0f)};
    const sbox::editor::PickResult hit = sbox::editor::pick_atom(ray, mol);
    EXPECT_EQ(hit.type, sbox::editor::PickResult::Type::None);
}

TEST(PickingTest, PickBondHitsMidpointCylinder) {
    sbox::chem::MolecularSystem mol;
    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "", 0});
    mol.add_atom({6, Eigen::Vector3d(5.0, 0.0, 0.0), "", 0});
    mol.add_bond(0, 1, sbox::chem::BondOrder::Single);

    const Eigen::Vector3f origin(2.5f, 2.0f, 5.0f);
    const Eigen::Vector3f direction = (Eigen::Vector3f(2.5f, 0.0f, 0.0f) - origin).normalized();
    const sbox::editor::Ray ray{origin, direction};
    const sbox::editor::PickResult hit = sbox::editor::pick_bond(ray, mol);
    EXPECT_EQ(hit.type, sbox::editor::PickResult::Type::Bond);
    EXPECT_EQ(hit.index, 0);
}

TEST(PickingTest, SelectionToggleAndClear) {
    sbox::editor::Selection selection;
    EXPECT_TRUE(selection.empty());

    selection.toggle_atom(3);
    EXPECT_TRUE(selection.has_atom(3));
    EXPECT_EQ(selection.num_atoms(), 1);

    selection.toggle_atom(3);
    EXPECT_FALSE(selection.has_atom(3));
    EXPECT_EQ(selection.num_atoms(), 0);

    selection.toggle_bond(1);
    EXPECT_TRUE(selection.has_bond(1));
    EXPECT_EQ(selection.num_bonds(), 1);

    selection.clear();
    EXPECT_TRUE(selection.empty());
}
