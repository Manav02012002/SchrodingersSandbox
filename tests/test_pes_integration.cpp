#include "backend/backend_manager.h"
#include "backend/python_env.h"
#include "core/molecular_system.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <thread>
#include <vector>

namespace {

constexpr double kAngstromToBohr = 1.8897259886;
constexpr double kPi = 3.14159265358979323846;
constexpr double kDegToRad = kPi / 180.0;
constexpr double kRadToDeg = 180.0 / kPi;

sbox::chem::MolecularSystem make_h2() {
    sbox::chem::MolecularSystem mol;
    mol.set_name("H2");
    mol.add_atom({1, Eigen::Vector3d(0.0, 0.0, -0.7), "", 0});
    mol.add_atom({1, Eigen::Vector3d(0.0, 0.0, 0.7), "", 0});
    mol.perceive_bonds();
    return mol;
}

sbox::chem::MolecularSystem make_distorted_h2o() {
    sbox::chem::MolecularSystem mol;
    mol.set_name("H2O distorted");
    mol.add_atom({8, Eigen::Vector3d(0.0, 0.0, 0.0), "", 0});
    mol.add_atom({1, Eigen::Vector3d(0.0, 1.43, -0.55), "", 0});
    mol.add_atom({1, Eigen::Vector3d(0.0, -1.43, -0.55), "", 0});
    mol.perceive_bonds();
    return mol;
}

sbox::chem::MolecularSystem make_h3_reactant() {
    sbox::chem::MolecularSystem mol;
    mol.set_name("H + H2 reactant");
    mol.add_atom({1, Eigen::Vector3d(-3.0, 0.0, 0.0) * kAngstromToBohr, "", 0});
    mol.add_atom({1, Eigen::Vector3d(0.0, 0.0, 0.0) * kAngstromToBohr, "", 0});
    mol.add_atom({1, Eigen::Vector3d(0.7, 0.0, 0.0) * kAngstromToBohr, "", 0});
    mol.perceive_bonds();
    return mol;
}

sbox::chem::MolecularSystem make_h3_product() {
    sbox::chem::MolecularSystem mol;
    mol.set_name("H2 + H product");
    mol.add_atom({1, Eigen::Vector3d(-0.7, 0.0, 0.0) * kAngstromToBohr, "", 0});
    mol.add_atom({1, Eigen::Vector3d(0.0, 0.0, 0.0) * kAngstromToBohr, "", 0});
    mol.add_atom({1, Eigen::Vector3d(3.0, 0.0, 0.0) * kAngstromToBohr, "", 0});
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

class PESIntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
        env_.detect();
        env_.check_packages();
        if (!env_.has_pyscf()) {
            GTEST_SKIP() << "PySCF not available, skipping PES/NEB integration tests";
        }
        backend_.init(env_);
    }

    void RequireGeometric() {
        if (!env_.info().has_geometric) {
            GTEST_SKIP() << "geomeTRIC not available, skipping optimization-based integration test";
        }
    }

    sbox::backend::PythonEnvironment env_;
    sbox::backend::BackendManager backend_;
};

}  // namespace

TEST_F(PESIntegrationTest, OneDimensionalPESScanOnH2) {
    sbox::backend::JobSpec spec;
    spec.geometry = make_h2();
    spec.method = sbox::backend::Method::HF;
    spec.basis = sbox::backend::BasisSetType::STO_3G;
    spec.charge = 0;
    spec.multiplicity = 1;
    spec.run_pes_scan = true;
    spec.scan.is_2d = false;
    spec.scan.coord1.type = sbox::backend::JobSpec::ScanSpec::CoordType::Distance;
    spec.scan.coord1.atom_indices = {0, 1};
    spec.scan.coord1.start = 0.5;
    spec.scan.coord1.end = 3.0;
    spec.scan.coord1.steps = 10;

    const int job_id = backend_.submit(spec);
    const auto& result = wait_for_job(backend_, job_id, std::chrono::seconds(120));

    ASSERT_TRUE(result.has_scan);
    ASSERT_EQ(result.scan_result.energies.size(), 10u);
    ASSERT_EQ(result.scan_result.coord1_values.size(), 10u);
    for (double energy : result.scan_result.energies) {
        EXPECT_TRUE(std::isfinite(energy));
    }

    auto min_it = std::min_element(result.scan_result.energies.begin(), result.scan_result.energies.end());
    ASSERT_NE(min_it, result.scan_result.energies.end());
    const std::size_t min_idx = static_cast<std::size_t>(std::distance(result.scan_result.energies.begin(), min_it));
    const double min_distance_ang = result.scan_result.coord1_values[min_idx];
    EXPECT_GT(min_distance_ang, 0.6);
    EXPECT_LT(min_distance_ang, 0.9);
    EXPECT_GT(result.scan_result.energies.back(), *min_it);
}

TEST_F(PESIntegrationTest, GeometryOptimizationWithTrajectory) {
    RequireGeometric();

    sbox::backend::JobSpec spec;
    spec.geometry = make_distorted_h2o();
    spec.method = sbox::backend::Method::HF;
    spec.basis = sbox::backend::BasisSetType::STO_3G;
    spec.charge = 0;
    spec.multiplicity = 1;
    spec.optimize_geometry = true;
    spec.properties = {sbox::backend::PropertyRequest::MoldenFile};

    const int job_id = backend_.submit(spec);
    const auto& result = wait_for_job(backend_, job_id, std::chrono::seconds(180));

    EXPECT_TRUE(result.converged());
    EXPECT_TRUE(result.optimization_converged);
    ASSERT_GT(result.opt_history.size(), 1u);
    EXPECT_TRUE(result.has_trajectory);
    EXPECT_EQ(result.trajectory_frames.size(), result.opt_history.size());

    const sbox::chem::MolecularSystem& final_geom =
        result.has_optimized_geometry ? result.optimized_geometry : result.trajectory_frames.back();
    const double final_angle = final_geom.angle(1, 0, 2) * kRadToDeg;
    EXPECT_GT(final_angle, 100.0);
    EXPECT_LT(final_angle, 110.0);

    ASSERT_FALSE(result.opt_history.empty());
    EXPECT_LT(result.opt_history.back().energy, result.opt_history.front().energy);
}

TEST_F(PESIntegrationTest, ConstrainedOptimizationAppliesAngleConstraint) {
    RequireGeometric();

    sbox::backend::JobSpec spec;
    spec.geometry = make_distorted_h2o();
    spec.method = sbox::backend::Method::HF;
    spec.basis = sbox::backend::BasisSetType::STO_3G;
    spec.charge = 0;
    spec.multiplicity = 1;
    spec.optimize_geometry = true;
    spec.properties = {sbox::backend::PropertyRequest::MoldenFile};
    spec.constraints.fixed_angles.emplace_back(1, 0, 2, 109.47 * kDegToRad);

    const int job_id = backend_.submit(spec);
    const auto& result = wait_for_job(backend_, job_id, std::chrono::seconds(180));

    EXPECT_TRUE(result.converged());
    const sbox::chem::MolecularSystem& final_geom =
        result.has_optimized_geometry ? result.optimized_geometry : result.trajectory_frames.back();
    const double final_angle = final_geom.angle(1, 0, 2) * kRadToDeg;
    EXPECT_NEAR(final_angle, 109.47, 1.0);
}

TEST_F(PESIntegrationTest, SlowNEBProducesPathAndBarrier) {
    const char* skip_slow = std::getenv("SBOX_SKIP_SLOW_TESTS");
    if (skip_slow != nullptr && std::string(skip_slow) == "1") {
        GTEST_SKIP() << "Skipping slow NEB test due to SBOX_SKIP_SLOW_TESTS=1";
    }

    sbox::backend::JobSpec spec;
    spec.method = sbox::backend::Method::HF;
    spec.basis = sbox::backend::BasisSetType::STO_3G;
    spec.charge = 0;
    spec.multiplicity = 1;
    spec.run_neb = true;
    spec.neb.reactant = make_h3_reactant();
    spec.neb.product = make_h3_product();
    spec.neb.num_images = 5;
    spec.neb.max_neb_steps = 10;

    const int job_id = backend_.submit(spec);
    const auto& result = wait_for_job(backend_, job_id, std::chrono::seconds(300));

    ASSERT_TRUE(result.has_neb);
    EXPECT_EQ(result.neb_result.path_energies.size(), 5u);
    EXPECT_EQ(result.neb_result.path_geometries.size(), 5u);
    EXPECT_GE(result.neb_result.ts_index, 1);
    EXPECT_LT(result.neb_result.ts_index, 4);
    EXPECT_GT(result.neb_result.forward_barrier, 0.0);
}
