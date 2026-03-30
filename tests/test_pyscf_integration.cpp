#include "backend/backend_manager.h"
#include "backend/python_env.h"
#include "core/molecular_system.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

#include <chrono>
#include <cmath>
#include <thread>
#include <vector>

namespace {

constexpr double kAngstromToBohr = 1.8897259886;

sbox::chem::MolecularSystem make_h2() {
    sbox::chem::MolecularSystem mol;
    mol.set_name("H2");
    mol.add_atom({1, Eigen::Vector3d(0.0, 0.0, -0.7), "", 0});
    mol.add_atom({1, Eigen::Vector3d(0.0, 0.0, 0.7), "", 0});
    mol.perceive_bonds();
    return mol;
}

sbox::chem::MolecularSystem make_h2o() {
    sbox::chem::MolecularSystem mol;
    mol.set_name("H2O");
    mol.add_atom({8, Eigen::Vector3d(0.0, 0.0, 0.117) * kAngstromToBohr, "", 0});
    mol.add_atom({1, Eigen::Vector3d(0.0, 0.757, -0.469) * kAngstromToBohr, "", 0});
    mol.add_atom({1, Eigen::Vector3d(0.0, -0.757, -0.469) * kAngstromToBohr, "", 0});
    mol.perceive_bonds();
    return mol;
}

sbox::chem::MolecularSystem make_ethanol() {
    sbox::chem::MolecularSystem mol;
    mol.set_name("Ethanol");
    const std::vector<std::pair<int, Eigen::Vector3d>> atoms = {
        {6, {0.000, 0.000, 0.000}},
        {6, {1.520, 0.000, 0.000}},
        {8, {2.040, 1.367, 0.000}},
        {1, {-0.360, -0.510, 0.890}},
        {1, {-0.360, -0.510, -0.890}},
        {1, {-0.360, 1.020, 0.000}},
        {1, {1.880, -0.510, 0.890}},
        {1, {1.880, -0.510, -0.890}},
        {1, {1.680, 1.850, 0.780}},
    };
    for (const auto& [z, pos_ang] : atoms) {
        mol.add_atom({z, pos_ang * kAngstromToBohr, "", 0});
    }
    mol.perceive_bonds();
    return mol;
}

sbox::chem::MolecularSystem make_benzene() {
    sbox::chem::MolecularSystem mol;
    mol.set_name("Benzene");
    const double pi = std::acos(-1.0);
    const double c_r = 1.397 * kAngstromToBohr;
    const double h_r = 2.479 * kAngstromToBohr;
    for (int i = 0; i < 6; ++i) {
        const double angle = i * pi / 3.0;
        mol.add_atom({6, Eigen::Vector3d(c_r * std::cos(angle), c_r * std::sin(angle), 0.0), "", 0});
    }
    for (int i = 0; i < 6; ++i) {
        const double angle = i * pi / 3.0;
        mol.add_atom({1, Eigen::Vector3d(h_r * std::cos(angle), h_r * std::sin(angle), 0.0), "", 0});
    }
    mol.perceive_bonds();
    return mol;
}

const sbox::backend::JobResult& wait_for_job(sbox::backend::BackendManager& backend,
                                             int job_id,
                                             std::chrono::seconds timeout) {
    const auto deadline = std::chrono::steady_clock::now() + timeout;
    while (std::chrono::steady_clock::now() < deadline) {
        const std::vector<int> completed = backend.poll_completed();
        for (int completed_id : completed) {
            if (completed_id == job_id) {
                const sbox::backend::JobResult* result = backend.result(job_id);
                EXPECT_NE(result, nullptr);
                return *result;
            }
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    FAIL() << "Timed out waiting for job " << job_id;
    static sbox::backend::JobResult unreachable;
    return unreachable;
}

class PySCFIntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
        env_.detect();
        env_.check_packages();
        if (!env_.has_pyscf()) {
            GTEST_SKIP() << "PySCF not available, skipping integration tests";
        }
        backend_.init(env_);
    }

    sbox::backend::PythonEnvironment env_;
    sbox::backend::BackendManager backend_;
};

}  // namespace

TEST_F(PySCFIntegrationTest, H2HfSto3g) {
    sbox::backend::JobSpec spec;
    spec.geometry = make_h2();
    spec.method = sbox::backend::Method::HF;
    spec.basis = sbox::backend::BasisSetType::STO_3G;
    spec.charge = 0;
    spec.multiplicity = 1;
    spec.properties = {
        sbox::backend::PropertyRequest::MullikenCharges,
        sbox::backend::PropertyRequest::DipoleMoment,
        sbox::backend::PropertyRequest::MoldenFile,
    };

    const int job_id = backend_.submit(spec);
    const auto& result = wait_for_job(backend_, job_id, std::chrono::seconds(60));

    EXPECT_TRUE(result.converged());
    EXPECT_GT(result.total_energy, -1.2);
    EXPECT_LT(result.total_energy, -1.0);
    ASSERT_EQ(result.mulliken_charges.size(), 2u);
    EXPECT_NEAR(result.mulliken_charges[0], 0.0, 1e-2);
    EXPECT_NEAR(result.mulliken_charges[1], 0.0, 1e-2);
    EXPECT_LT(result.dipole_moment.norm(), 0.01);
    EXPECT_TRUE(result.has_mo_data);
    ASSERT_EQ(result.mo_data.coefficients.rows(), 2);
    ASSERT_EQ(result.mo_data.coefficients.cols(), 2);
    EXPECT_EQ(result.homo_index(), 0);
    EXPECT_EQ(result.lumo_index(), 1);
}

TEST_F(PySCFIntegrationTest, H2OHfSto3g) {
    sbox::backend::JobSpec spec;
    spec.geometry = make_h2o();
    spec.method = sbox::backend::Method::HF;
    spec.basis = sbox::backend::BasisSetType::STO_3G;
    spec.charge = 0;
    spec.multiplicity = 1;
    spec.properties = {
        sbox::backend::PropertyRequest::MullikenCharges,
        sbox::backend::PropertyRequest::DipoleMoment,
        sbox::backend::PropertyRequest::MoldenFile,
    };

    const int job_id = backend_.submit(spec);
    const auto& result = wait_for_job(backend_, job_id, std::chrono::seconds(60));

    EXPECT_TRUE(result.converged());
    EXPECT_GT(result.total_energy, -76.0);
    EXPECT_LT(result.total_energy, -74.0);
    ASSERT_EQ(result.mulliken_charges.size(), 3u);
    EXPECT_LT(result.mulliken_charges[0], 0.0);
    EXPECT_GT(result.mulliken_charges[1], 0.0);
    EXPECT_GT(result.mulliken_charges[2], 0.0);
    EXPECT_GT(result.dipole_moment.norm(), 1.0);
    EXPECT_TRUE(result.has_mo_data);
    EXPECT_EQ(result.mo_data.coefficients.rows(), 7);
    EXPECT_GT(result.homo_lumo_gap_eV(), 10.0);
}

TEST_F(PySCFIntegrationTest, EthanolB3LYP631Gd) {
    sbox::backend::JobSpec spec;
    spec.geometry = make_ethanol();
    spec.method = sbox::backend::Method::DFT_B3LYP;
    spec.basis = sbox::backend::BasisSetType::B6_31Gd;
    spec.charge = 0;
    spec.multiplicity = 1;
    spec.properties = {
        sbox::backend::PropertyRequest::MullikenCharges,
        sbox::backend::PropertyRequest::DipoleMoment,
        sbox::backend::PropertyRequest::MoldenFile,
    };

    const int job_id = backend_.submit(spec);
    const auto& result = wait_for_job(backend_, job_id, std::chrono::seconds(120));

    EXPECT_TRUE(result.converged());
    EXPECT_LT(result.total_energy, -150.0);
    EXPECT_TRUE(result.has_mo_data);
    EXPECT_GT(result.mo_data.coefficients.cols(), 20);
}

TEST_F(PySCFIntegrationTest, JobCancellation) {
    sbox::backend::JobSpec spec;
    spec.geometry = make_benzene();
    spec.method = sbox::backend::Method::HF;
    spec.basis = sbox::backend::BasisSetType::cc_pVTZ;
    spec.charge = 0;
    spec.multiplicity = 1;
    spec.properties = {
        sbox::backend::PropertyRequest::MullikenCharges,
        sbox::backend::PropertyRequest::DipoleMoment,
        sbox::backend::PropertyRequest::MoldenFile,
    };

    const int job_id = backend_.submit(spec);
    std::this_thread::sleep_for(std::chrono::seconds(2));
    backend_.cancel(job_id);
    const auto& result = wait_for_job(backend_, job_id, std::chrono::seconds(60));

    EXPECT_EQ(result.status, sbox::backend::JobStatus::Cancelled);
}

TEST_F(PySCFIntegrationTest, InvalidMethodConfigurationFails) {
    sbox::backend::JobSpec spec;
    spec.geometry.clear();
    spec.geometry.set_name("H");
    spec.geometry.add_atom({1, Eigen::Vector3d::Zero(), "", 0});
    spec.method = sbox::backend::Method::HF;
    spec.basis = sbox::backend::BasisSetType::STO_3G;
    spec.charge = 0;
    spec.multiplicity = 4;
    spec.properties = {
        sbox::backend::PropertyRequest::MullikenCharges,
        sbox::backend::PropertyRequest::DipoleMoment,
        sbox::backend::PropertyRequest::MoldenFile,
    };

    const int job_id = backend_.submit(spec);
    const auto& result = wait_for_job(backend_, job_id, std::chrono::seconds(60));

    EXPECT_EQ(result.status, sbox::backend::JobStatus::Failed);
    EXPECT_FALSE(result.error_message.empty());
}
