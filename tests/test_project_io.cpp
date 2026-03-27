#include "core/molecular_system.h"
#include "io/project_io.h"

#include <json.hpp>

#include <gtest/gtest.h>

#include <Eigen/Core>

#include <filesystem>
#include <fstream>

namespace {

std::filesystem::path temp_project_path(const char* stem) {
    return std::filesystem::temp_directory_path()
        / (std::string("schrodingerssandbox_") + stem + ".json");
}

sbox::chem::MolecularSystem make_water() {
    sbox::chem::MolecularSystem mol;
    mol.set_name("water");
    mol.set_charge(0);
    mol.set_multiplicity(1);
    mol.add_atom({8, Eigen::Vector3d(0.0, 0.0, 0.0), "O"});
    mol.add_atom({1, Eigen::Vector3d(1.81, 0.0, 0.0), "H1"});
    mol.add_atom({1, Eigen::Vector3d(-0.45, 0.0, 1.75), "H2"});
    mol.add_bond(0, 1, sbox::chem::BondOrder::Single);
    mol.add_bond(0, 2, sbox::chem::BondOrder::Single);
    return mol;
}

}  // namespace

TEST(ProjectIoTest, RoundTripPreservesWater) {
    const sbox::chem::MolecularSystem original = make_water();
    const std::filesystem::path path = temp_project_path("project_roundtrip");

    sbox::io::save_project(path.string(), original);
    const sbox::chem::MolecularSystem loaded = sbox::io::load_project(path.string());
    std::filesystem::remove(path);

    EXPECT_EQ(loaded.name(), original.name());
    EXPECT_EQ(loaded.charge(), original.charge());
    EXPECT_EQ(loaded.multiplicity(), original.multiplicity());
    ASSERT_EQ(loaded.num_atoms(), original.num_atoms());
    ASSERT_EQ(loaded.num_bonds(), original.num_bonds());
    for (int i = 0; i < original.num_atoms(); ++i) {
        EXPECT_EQ(loaded.atom(i).Z, original.atom(i).Z);
        EXPECT_NEAR(loaded.atom(i).position.x(), original.atom(i).position.x(), 1e-12);
        EXPECT_NEAR(loaded.atom(i).position.y(), original.atom(i).position.y(), 1e-12);
        EXPECT_NEAR(loaded.atom(i).position.z(), original.atom(i).position.z(), 1e-12);
    }
}

TEST(ProjectIoTest, PreservesExtraState) {
    const sbox::chem::MolecularSystem original = make_water();
    const std::filesystem::path path = temp_project_path("project_state");
    const nlohmann::json extra_state = {
        {"render_mode", 2},
        {"iso_value", 0.01},
    };

    sbox::io::save_project(path.string(), original, extra_state);
    nlohmann::json loaded_state;
    (void)sbox::io::load_project(path.string(), &loaded_state);
    std::filesystem::remove(path);

    EXPECT_EQ(loaded_state.at("render_mode").get<int>(), 2);
    EXPECT_DOUBLE_EQ(loaded_state.at("iso_value").get<double>(), 0.01);
}

TEST(ProjectIoTest, MalformedJsonThrows) {
    const std::filesystem::path path = temp_project_path("project_malformed");
    {
        std::ofstream out(path);
        out << "{ invalid json";
    }

    EXPECT_THROW((void)sbox::io::load_project(path.string()), std::runtime_error);
    std::filesystem::remove(path);
}

TEST(ProjectIoTest, PreservesBondOrders) {
    sbox::chem::MolecularSystem mol;
    mol.set_name("bond orders");
    mol.add_atom({6, Eigen::Vector3d(0.0, 0.0, 0.0), "A"});
    mol.add_atom({6, Eigen::Vector3d(1.0, 0.0, 0.0), "B"});
    mol.add_atom({6, Eigen::Vector3d(2.0, 0.0, 0.0), "C"});
    mol.add_atom({6, Eigen::Vector3d(3.0, 0.0, 0.0), "D"});
    mol.add_bond(0, 1, sbox::chem::BondOrder::Single);
    mol.add_bond(1, 2, sbox::chem::BondOrder::Double);
    mol.add_bond(2, 3, sbox::chem::BondOrder::Triple);

    const std::filesystem::path path = temp_project_path("project_bond_orders");
    sbox::io::save_project(path.string(), mol);
    const sbox::chem::MolecularSystem loaded = sbox::io::load_project(path.string());
    std::filesystem::remove(path);

    ASSERT_EQ(loaded.num_bonds(), 3);
    EXPECT_EQ(loaded.bond(0).order, sbox::chem::BondOrder::Single);
    EXPECT_EQ(loaded.bond(1).order, sbox::chem::BondOrder::Double);
    EXPECT_EQ(loaded.bond(2).order, sbox::chem::BondOrder::Triple);
}
