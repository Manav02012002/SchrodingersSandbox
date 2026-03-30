#pragma once

#include <string>
#include <vector>

#include <Eigen/Core>

#include "core/basis_set.h"
#include "core/molecular_system.h"
#include "io/cube_io.h"

namespace sbox::backend {

enum class Method {
    HF,
    UHF,
    DFT_B3LYP,
    DFT_PBE,
    DFT_PBE0,
    DFT_TPSS,
    DFT_M06_2X,
    MP2,
    CCSD,
    GFN2_XTB,
    GFN1_XTB,
    GFN_FF,
};

enum class BasisSetType {
    STO_3G,
    B3_21G,
    B6_31G,
    B6_31Gd,
    B6_31Gdp,
    B6_311Gdp,
    cc_pVDZ,
    cc_pVTZ,
    cc_pVQZ,
    def2_SVP,
    def2_TZVP,
    def2_TZVPP,
    aug_cc_pVDZ,
    aug_cc_pVTZ,
};

enum class JobStatus {
    Pending,
    Running,
    Converged,
    Failed,
    Cancelled,
    Timeout,
};

enum class PropertyRequest {
    MullikenCharges,
    LowdinCharges,
    DipoleMoment,
    MayerBondOrders,
    MoldenFile,
    CubeHOMO,
    CubeLUMO,
    CubeDensity,
    CubeESP,
    Frequencies,
    Optimization,
};

inline const char* method_to_string(Method m) {
    switch (m) {
    case Method::HF: return "hf";
    case Method::UHF: return "uhf";
    case Method::DFT_B3LYP: return "b3lyp";
    case Method::DFT_PBE: return "pbe";
    case Method::DFT_PBE0: return "pbe0";
    case Method::DFT_TPSS: return "tpss";
    case Method::DFT_M06_2X: return "m06-2x";
    case Method::MP2: return "mp2";
    case Method::CCSD: return "ccsd";
    case Method::GFN2_XTB: return "gfn2-xtb";
    case Method::GFN1_XTB: return "gfn1-xtb";
    case Method::GFN_FF: return "gfn-ff";
    }
    return "hf";
}

inline const char* basis_to_string(BasisSetType b) {
    switch (b) {
    case BasisSetType::STO_3G: return "sto-3g";
    case BasisSetType::B3_21G: return "3-21g";
    case BasisSetType::B6_31G: return "6-31g";
    case BasisSetType::B6_31Gd: return "6-31g*";
    case BasisSetType::B6_31Gdp: return "6-31g**";
    case BasisSetType::B6_311Gdp: return "6-311g**";
    case BasisSetType::cc_pVDZ: return "cc-pvdz";
    case BasisSetType::cc_pVTZ: return "cc-pvtz";
    case BasisSetType::cc_pVQZ: return "cc-pvqz";
    case BasisSetType::def2_SVP: return "def2-svp";
    case BasisSetType::def2_TZVP: return "def2-tzvp";
    case BasisSetType::def2_TZVPP: return "def2-tzvpp";
    case BasisSetType::aug_cc_pVDZ: return "aug-cc-pvdz";
    case BasisSetType::aug_cc_pVTZ: return "aug-cc-pvtz";
    }
    return "sto-3g";
}

inline const char* method_display_name(Method m) {
    switch (m) {
    case Method::HF: return "Hartree-Fock";
    case Method::UHF: return "Unrestricted Hartree-Fock";
    case Method::DFT_B3LYP: return "DFT (B3LYP)";
    case Method::DFT_PBE: return "DFT (PBE)";
    case Method::DFT_PBE0: return "DFT (PBE0)";
    case Method::DFT_TPSS: return "DFT (TPSS)";
    case Method::DFT_M06_2X: return "DFT (M06-2X)";
    case Method::MP2: return "MP2";
    case Method::CCSD: return "CCSD";
    case Method::GFN2_XTB: return "GFN2-xTB";
    case Method::GFN1_XTB: return "GFN1-xTB";
    case Method::GFN_FF: return "GFN-FF";
    }
    return "Hartree-Fock";
}

inline const char* basis_display_name(BasisSetType b) {
    switch (b) {
    case BasisSetType::STO_3G: return "STO-3G";
    case BasisSetType::B3_21G: return "3-21G";
    case BasisSetType::B6_31G: return "6-31G";
    case BasisSetType::B6_31Gd: return "6-31G*";
    case BasisSetType::B6_31Gdp: return "6-31G**";
    case BasisSetType::B6_311Gdp: return "6-311G**";
    case BasisSetType::cc_pVDZ: return "cc-pVDZ";
    case BasisSetType::cc_pVTZ: return "cc-pVTZ";
    case BasisSetType::cc_pVQZ: return "cc-pVQZ";
    case BasisSetType::def2_SVP: return "def2-SVP";
    case BasisSetType::def2_TZVP: return "def2-TZVP";
    case BasisSetType::def2_TZVPP: return "def2-TZVPP";
    case BasisSetType::aug_cc_pVDZ: return "aug-cc-pVDZ";
    case BasisSetType::aug_cc_pVTZ: return "aug-cc-pVTZ";
    }
    return "STO-3G";
}

inline bool method_needs_basis(Method m) {
    return m != Method::GFN2_XTB && m != Method::GFN1_XTB && m != Method::GFN_FF;
}

inline bool method_is_xtb(Method m) {
    return m == Method::GFN2_XTB || m == Method::GFN1_XTB || m == Method::GFN_FF;
}

struct JobSpec {
    sbox::chem::MolecularSystem geometry;

    Method method = Method::HF;
    BasisSetType basis = BasisSetType::STO_3G;

    int charge = 0;
    int multiplicity = 1;

    int max_scf_cycles = 200;
    double scf_convergence = 1e-8;

    std::vector<PropertyRequest> properties = {
        PropertyRequest::MullikenCharges,
        PropertyRequest::DipoleMoment,
        PropertyRequest::MoldenFile,
    };

    bool optimize_geometry = false;
    int max_opt_steps = 100;
    double opt_convergence = 1e-5;

    std::string solvent;

    std::string work_dir;

    int job_id = 0;
};

struct SCFIteration {
    int iteration = 0;
    double energy = 0.0;
    double delta_energy = 0.0;
    double gradient_norm = 0.0;
    double wall_time_ms = 0.0;
};

struct OptimizationStep {
    int step = 0;
    double energy = 0.0;
    double gradient_norm = 0.0;
    sbox::chem::MolecularSystem geometry;
};

struct JobResult {
    int job_id = 0;
    JobStatus status = JobStatus::Pending;
    std::string error_message;
    std::string work_dir;

    double total_energy = 0.0;

    sbox::basis::MOData mo_data;
    bool has_mo_data = false;

    sbox::chem::MolecularSystem optimized_geometry;
    bool has_optimized_geometry = false;

    Eigen::Vector3d dipole_moment = Eigen::Vector3d::Zero();
    std::vector<double> mulliken_charges;
    std::vector<double> lowdin_charges;
    Eigen::MatrixXd mayer_bond_orders;

    sbox::io::CubeData homo_cube;
    sbox::io::CubeData lumo_cube;
    sbox::io::CubeData density_cube;
    sbox::io::CubeData esp_cube;
    bool has_homo_cube = false;
    bool has_lumo_cube = false;
    bool has_density_cube = false;
    bool has_esp_cube = false;

    std::vector<double> frequencies_cm1;
    std::vector<double> ir_intensities;
    std::vector<Eigen::VectorXd> normal_modes;
    bool has_frequencies = false;

    std::vector<SCFIteration> scf_history;
    std::vector<OptimizationStep> opt_history;

    double wall_time_seconds = 0.0;

    bool converged() const { return status == JobStatus::Converged; }

    int homo_index() const {
        if (!has_mo_data || mo_data.occupations.size() == 0) {
            return -1;
        }
        int homo = -1;
        for (int i = 0; i < mo_data.occupations.size(); ++i) {
            if (mo_data.occupations(i) > 0.5) {
                homo = i;
            }
        }
        return homo;
    }

    int lumo_index() const {
        if (!has_mo_data || mo_data.energies.size() == 0) {
            return -1;
        }
        const int homo = homo_index();
        if (homo < 0) {
            return mo_data.energies.size() > 0 ? 0 : -1;
        }
        const int lumo = homo + 1;
        return lumo < mo_data.energies.size() ? lumo : -1;
    }

    double homo_lumo_gap_eV() const {
        const int homo = homo_index();
        const int lumo = lumo_index();
        if (homo < 0 || lumo < 0) {
            return 0.0;
        }
        return (mo_data.energies(lumo) - mo_data.energies(homo)) * 27.2114;
    }
};

}  // namespace sbox::backend
