#include "backend/backend_manager.h"

#include "core/molden_parser.h"

#include <json.hpp>

#include <Eigen/Core>

#include <chrono>
#include <csignal>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unistd.h>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/wait.h>

namespace sbox::backend {

namespace {

using json = nlohmann::json;

std::string property_to_string(PropertyRequest property) {
    switch (property) {
    case PropertyRequest::MullikenCharges: return "mulliken";
    case PropertyRequest::LowdinCharges: return "lowdin";
    case PropertyRequest::DipoleMoment: return "dipole";
    case PropertyRequest::MayerBondOrders: return "mayer";
    case PropertyRequest::MoldenFile: return "molden";
    case PropertyRequest::CubeHOMO: return "cube_homo";
    case PropertyRequest::CubeLUMO: return "cube_lumo";
    case PropertyRequest::CubeDensity: return "cube_density";
    case PropertyRequest::CubeESP: return "cube_esp";
    case PropertyRequest::Frequencies: return "frequencies";
    case PropertyRequest::Optimization: return "optimization";
    }
    return {};
}

std::filesystem::path detect_scripts_dir() {
    const std::filesystem::path cwd = std::filesystem::current_path();
    const std::vector<std::filesystem::path> candidates = {
        cwd / "data/scripts",
        cwd / "../data/scripts",
        cwd / "../../data/scripts",
    };
    for (const auto& candidate : candidates) {
        if (std::filesystem::exists(candidate / "pyscf_driver.py")) {
            return std::filesystem::weakly_canonical(candidate);
        }
    }
    return cwd / "data/scripts";
}

double vec3_component_or_zero(const json& j, std::size_t index) {
    if (!j.is_array() || j.size() <= index) {
        return 0.0;
    }
    return j[index].get<double>();
}

json load_json_file(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open JSON file: " + path.string());
    }
    json parsed;
    in >> parsed;
    return parsed;
}

sbox::chem::MolecularSystem geometry_from_json(const json& atoms_json) {
    sbox::chem::MolecularSystem mol;
    for (const auto& atom_entry : atoms_json) {
        if (!atom_entry.is_array() || atom_entry.size() != 2 || !atom_entry[1].is_array() || atom_entry[1].size() != 3) {
            throw std::runtime_error("Invalid optimized geometry entry in result.json");
        }
        sbox::chem::Atom atom;
        atom.Z = atom_entry[0].get<int>();
        atom.position = Eigen::Vector3d(
            atom_entry[1][0].get<double>(),
            atom_entry[1][1].get<double>(),
            atom_entry[1][2].get<double>());
        mol.add_atom(atom);
    }
    mol.perceive_bonds();
    return mol;
}

}  // namespace

BackendManager::BackendManager()
    : scripts_dir_(detect_scripts_dir().string()) {}

BackendManager::~BackendManager() {
    std::vector<RunningJob*> jobs;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        for (auto& [job_id, job] : running_jobs_) {
            (void)job_id;
            job->cancelled.store(true);
            if (job->pid > 0) {
                ::kill(job->pid, SIGTERM);
            }
            jobs.push_back(job.get());
        }
    }

    for (RunningJob* job : jobs) {
        if (job->future.valid()) {
            job->future.wait();
        }
    }
}

void BackendManager::init(const PythonEnvironment& env) {
    python_env_ = env;
}

int BackendManager::submit(const JobSpec& spec) {
    JobSpec job_spec = spec;

    std::unique_ptr<RunningJob> job = std::make_unique<RunningJob>();
    RunningJob* job_ptr = nullptr;
    int job_id = 0;

    {
        std::lock_guard<std::mutex> lock(mutex_);
        job->job_id = next_job_id_++;
        job_id = job->job_id;
    }

    job_spec.job_id = job_id;
    job->spec = job_spec;
    job->work_dir = create_work_dir(job_id);
    job->spec.work_dir = job->work_dir;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        auto [it, inserted] = running_jobs_.emplace(job_id, std::move(job));
        (void)inserted;
        job_ptr = it->second.get();
    }

    job_ptr->future = std::async(std::launch::async, [this, job_ptr]() {
        return run_job(job_ptr->spec, job_ptr->work_dir, job_ptr->cancelled);
    });

    return job_id;
}

bool BackendManager::is_running(int job_id) const {
    std::lock_guard<std::mutex> lock(mutex_);
    return running_jobs_.find(job_id) != running_jobs_.end();
}

JobStatus BackendManager::status(int job_id) const {
    std::lock_guard<std::mutex> lock(mutex_);
    const auto running_it = running_jobs_.find(job_id);
    if (running_it != running_jobs_.end()) {
        return JobStatus::Running;
    }
    const auto completed_it = completed_jobs_.find(job_id);
    if (completed_it != completed_jobs_.end()) {
        return completed_it->second.status;
    }
    return JobStatus::Pending;
}

const JobResult* BackendManager::result(int job_id) const {
    std::lock_guard<std::mutex> lock(mutex_);
    const auto it = completed_jobs_.find(job_id);
    return it != completed_jobs_.end() ? &it->second : nullptr;
}

std::vector<int> BackendManager::poll_completed() {
    std::lock_guard<std::mutex> lock(mutex_);

    for (auto it = running_jobs_.begin(); it != running_jobs_.end();) {
        if (it->second->future.valid()
            && it->second->future.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
            JobResult completed = it->second->future.get();
            completed_jobs_[it->first] = completed;
            newly_completed_.push_back(it->first);
            it = running_jobs_.erase(it);
        } else {
            ++it;
        }
    }

    std::vector<int> completed;
    completed.swap(newly_completed_);
    return completed;
}

BackendManager::Progress BackendManager::get_progress(int job_id) const {
    std::string work_dir;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        const auto it = running_jobs_.find(job_id);
        if (it == running_jobs_.end()) {
            return {};
        }
        work_dir = it->second->work_dir;
    }
    return parse_progress(work_dir);
}

void BackendManager::cancel(int job_id) {
    std::lock_guard<std::mutex> lock(mutex_);
    const auto it = running_jobs_.find(job_id);
    if (it == running_jobs_.end()) {
        return;
    }
    it->second->cancelled.store(true);
    if (it->second->pid > 0) {
        ::kill(it->second->pid, SIGTERM);
    }
}

void BackendManager::clear_job(int job_id) {
    std::lock_guard<std::mutex> lock(mutex_);
    completed_jobs_.erase(job_id);
}

bool BackendManager::can_run_pyscf() const {
    return python_env_.is_valid() && python_env_.has_pyscf()
        && std::filesystem::exists(std::filesystem::path(scripts_dir_) / "pyscf_driver.py");
}

bool BackendManager::can_run_xtb() const {
    const bool have_runtime = python_env_.is_valid() && (python_env_.has_tblite() || python_env_.info().has_xtb);
    return have_runtime
        && std::filesystem::exists(std::filesystem::path(scripts_dir_) / "xtb_driver.py");
}

std::string BackendManager::scripts_dir() const {
    return scripts_dir_;
}

JobResult BackendManager::run_job(const JobSpec& spec, const std::string& work_dir, std::atomic<bool>& cancelled) {
    JobResult result;
    result.job_id = spec.job_id;
    result.status = JobStatus::Running;
    result.work_dir = work_dir;

    try {
        if (!python_env_.is_valid()) {
            throw std::runtime_error("No valid Python environment configured");
        }

        write_job_json(spec, work_dir);
        const std::filesystem::path script_path = driver_script(spec.method);
        if (!std::filesystem::exists(script_path)) {
            throw std::runtime_error("Driver script not found: " + script_path.string());
        }

        const std::filesystem::path job_json_path = std::filesystem::path(work_dir) / "job.json";
        const std::filesystem::path log_path = std::filesystem::path(work_dir) / "subprocess.log";

        pid_t pid = ::fork();
        if (pid < 0) {
            throw std::runtime_error("Failed to fork backend driver process");
        }

        if (pid == 0) {
            const int log_fd = ::open(log_path.c_str(), O_CREAT | O_WRONLY | O_TRUNC, 0644);
            if (log_fd >= 0) {
                ::dup2(log_fd, STDOUT_FILENO);
                ::dup2(log_fd, STDERR_FILENO);
                ::close(log_fd);
            }

            ::chdir(work_dir.c_str());
            ::execl(
                python_env_.info().python_path.c_str(),
                python_env_.info().python_path.c_str(),
                script_path.c_str(),
                job_json_path.c_str(),
                static_cast<char*>(nullptr));
            _exit(127);
        }

        {
            std::lock_guard<std::mutex> lock(mutex_);
            const auto it = running_jobs_.find(spec.job_id);
            if (it != running_jobs_.end()) {
                it->second->pid = pid;
            }
        }

        int wait_status = 0;
        while (true) {
            const pid_t wait_rc = ::waitpid(pid, &wait_status, WNOHANG);
            if (wait_rc == pid) {
                break;
            }
            if (wait_rc < 0) {
                throw std::runtime_error("Failed while waiting for backend driver process");
            }

            if (cancelled.load()) {
                ::kill(pid, SIGTERM);
                ::waitpid(pid, &wait_status, 0);
                result.status = JobStatus::Cancelled;
                result.error_message = "Job cancelled";
                break;
            }

            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }

        {
            std::lock_guard<std::mutex> lock(mutex_);
            const auto it = running_jobs_.find(spec.job_id);
            if (it != running_jobs_.end()) {
                it->second->pid = -1;
            }
        }

        if (result.status != JobStatus::Cancelled) {
            if (WIFEXITED(wait_status) && WEXITSTATUS(wait_status) == 0) {
                result = parse_result(spec, work_dir);
            } else if (cancelled.load()) {
                result.status = JobStatus::Cancelled;
                result.error_message = "Job cancelled";
            } else {
                result = parse_result(spec, work_dir);
                if (result.status == JobStatus::Pending || result.status == JobStatus::Running) {
                    result.status = JobStatus::Failed;
                }
                if (result.error_message.empty()) {
                    result.error_message = "Backend driver exited abnormally";
                }
            }
        }
    } catch (const std::exception& e) {
        result.status = cancelled.load() ? JobStatus::Cancelled : JobStatus::Failed;
        result.error_message = e.what();
    }

    result.job_id = spec.job_id;

    return result;
}

std::string BackendManager::create_work_dir(int job_id) {
    const std::filesystem::path base = std::filesystem::temp_directory_path() / "sbox";
    std::filesystem::create_directories(base);

    std::filesystem::path work_dir = base / ("job_" + std::to_string(job_id));
    int suffix = 1;
    while (std::filesystem::exists(work_dir)) {
        work_dir = base / ("job_" + std::to_string(job_id) + "_" + std::to_string(suffix++));
    }
    std::filesystem::create_directories(work_dir);
    return work_dir.string();
}

void BackendManager::write_job_json(const JobSpec& spec, const std::string& work_dir) {
    json j;
    j["geometry"] = json::array();
    for (const auto& atom : spec.geometry.atoms()) {
        j["geometry"].push_back({
            atom.Z,
            {atom.position.x(), atom.position.y(), atom.position.z()},
        });
    }

    j["method"] = method_to_string(spec.method);
    j["basis"] = basis_to_string(spec.basis);
    j["charge"] = spec.charge;
    j["multiplicity"] = spec.multiplicity;
    j["max_scf_cycles"] = spec.max_scf_cycles;
    j["scf_convergence"] = spec.scf_convergence;
    j["properties"] = json::array();
    bool request_frequencies = false;
    for (PropertyRequest property : spec.properties) {
        j["properties"].push_back(property_to_string(property));
        request_frequencies = request_frequencies || property == PropertyRequest::Frequencies;
    }
    j["frequencies"] = request_frequencies;
    j["optimize"] = spec.optimize_geometry;
    j["max_opt_steps"] = spec.max_opt_steps;
    j["opt_convergence"] = spec.opt_convergence;
    j["solvent"] = spec.solvent;
    j["output_dir"] = work_dir;
    j["cube_resolution"] = 80;

    if (request_frequencies && !spec.optimize_geometry) {
        std::cerr << "Warning: frequency calculation requested without geometry optimization; "
                     "results will use the current geometry.\n";
    }

    std::ofstream out(std::filesystem::path(work_dir) / "job.json");
    if (!out) {
        throw std::runtime_error("Failed to write job.json to " + work_dir);
    }
    out << j.dump(2);
}

JobResult BackendManager::parse_result(const JobSpec& spec, const std::string& work_dir) {
    JobResult result;
    result.job_id = spec.job_id;
    result.work_dir = work_dir;

    const std::filesystem::path result_path = std::filesystem::path(work_dir) / "result.json";
    const json j = load_json_file(result_path);

    const bool success = j.value("success", false);
    result.status = success ? JobStatus::Converged : JobStatus::Failed;
    result.error_message = j.value("error", std::string{});
    result.total_energy = j.value("total_energy", 0.0);
    result.wall_time_seconds = j.value("wall_time", 0.0);

    result.dipole_moment = Eigen::Vector3d(
        vec3_component_or_zero(j.value("dipole_moment", json::array()), 0),
        vec3_component_or_zero(j.value("dipole_moment", json::array()), 1),
        vec3_component_or_zero(j.value("dipole_moment", json::array()), 2));

    if (j.contains("mulliken_charges")) {
        result.mulliken_charges = j["mulliken_charges"].get<std::vector<double>>();
    }
    if (j.contains("lowdin_charges")) {
        result.lowdin_charges = j["lowdin_charges"].get<std::vector<double>>();
    }

    if (j.contains("mayer_bond_orders") && j["mayer_bond_orders"].is_array() && !j["mayer_bond_orders"].empty()) {
        const std::size_t rows = j["mayer_bond_orders"].size();
        const std::size_t cols = j["mayer_bond_orders"][0].size();
        result.mayer_bond_orders = Eigen::MatrixXd::Zero(static_cast<int>(rows), static_cast<int>(cols));
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t k = 0; k < cols; ++k) {
                result.mayer_bond_orders(static_cast<int>(i), static_cast<int>(k)) =
                    j["mayer_bond_orders"][i][k].get<double>();
            }
        }
    }

    if (j.contains("orbital_energies") && j["orbital_energies"].is_array()) {
        const auto energies = j["orbital_energies"].get<std::vector<double>>();
        result.mo_data.energies.resize(static_cast<int>(energies.size()));
        for (int i = 0; i < result.mo_data.energies.size(); ++i) {
            result.mo_data.energies(i) = energies[static_cast<std::size_t>(i)];
        }
        result.has_mo_data = result.mo_data.energies.size() > 0;
    }
    if (j.contains("orbital_occupations") && j["orbital_occupations"].is_array()) {
        const auto occupations = j["orbital_occupations"].get<std::vector<double>>();
        result.mo_data.occupations.resize(static_cast<int>(occupations.size()));
        for (int i = 0; i < result.mo_data.occupations.size(); ++i) {
            result.mo_data.occupations(i) = occupations[static_cast<std::size_t>(i)];
        }
        result.has_mo_data = result.has_mo_data || result.mo_data.occupations.size() > 0;
    }

    if (j.contains("scf_history") && j["scf_history"].is_array()) {
        for (const auto& entry : j["scf_history"]) {
            SCFIteration iter;
            iter.iteration = entry.value("iteration", 0);
            iter.energy = entry.value("energy", 0.0);
            iter.delta_energy = entry.value("delta_energy", 0.0);
            iter.gradient_norm = entry.value("gradient_norm", 0.0);
            iter.wall_time_ms = entry.value("wall_time_ms", 0.0);
            result.scf_history.push_back(iter);
        }
    }

    if (j.contains("frequencies_cm1") && j["frequencies_cm1"].is_array()) {
        result.frequencies_cm1 = j["frequencies_cm1"].get<std::vector<double>>();
        result.has_frequencies = !result.frequencies_cm1.empty();
    } else if (j.contains("frequencies") && j["frequencies"].is_array()) {
        result.frequencies_cm1 = j["frequencies"].get<std::vector<double>>();
        result.has_frequencies = !result.frequencies_cm1.empty();
    }
    if (j.contains("ir_intensities") && j["ir_intensities"].is_array()) {
        result.ir_intensities = j["ir_intensities"].get<std::vector<double>>();
        result.has_frequencies = result.has_frequencies || !result.ir_intensities.empty();
    }
    if (j.contains("normal_modes") && j["normal_modes"].is_array()) {
        for (const auto& mode : j["normal_modes"]) {
            if (!mode.is_array()) {
                continue;
            }
            Eigen::VectorXd vec(static_cast<int>(mode.size()));
            for (int i = 0; i < vec.size(); ++i) {
                vec(i) = mode[static_cast<std::size_t>(i)].get<double>();
            }
            result.normal_modes.push_back(vec);
        }
    }

    const std::filesystem::path molden_path = std::filesystem::path(work_dir) / "result.molden";
    if (std::filesystem::exists(molden_path)) {
        try {
            result.mo_data = sbox::molden::parse_molden_file(molden_path.string());
            result.has_mo_data = true;
            result.mo_data.total_energy = result.total_energy;
        } catch (const std::exception& e) {
            if (result.error_message.empty()) {
                result.error_message = e.what();
            } else {
                result.error_message += "\n";
                result.error_message += e.what();
            }
        }
    }

    const std::filesystem::path density_path = std::filesystem::path(work_dir) / "density.cube";
    if (std::filesystem::exists(density_path)) {
        result.density_cube = sbox::io::read_cube(density_path.string());
        result.has_density_cube = true;
    }

    const std::filesystem::path homo_path = std::filesystem::path(work_dir) / "homo.cube";
    if (std::filesystem::exists(homo_path)) {
        result.homo_cube = sbox::io::read_cube(homo_path.string());
        result.has_homo_cube = true;
    }

    const std::filesystem::path lumo_path = std::filesystem::path(work_dir) / "lumo.cube";
    if (std::filesystem::exists(lumo_path)) {
        result.lumo_cube = sbox::io::read_cube(lumo_path.string());
        result.has_lumo_cube = true;
    }

    const std::filesystem::path esp_path = std::filesystem::path(work_dir) / "esp.cube";
    if (std::filesystem::exists(esp_path)) {
        result.esp_cube = sbox::io::read_cube(esp_path.string());
        result.has_esp_cube = true;
    }

    if (j.contains("optimized_geometry") && j["optimized_geometry"].is_array()) {
        result.optimized_geometry = geometry_from_json(j["optimized_geometry"]);
        result.optimized_geometry.set_charge(spec.charge);
        result.optimized_geometry.set_multiplicity(spec.multiplicity);
        result.has_optimized_geometry = true;
    }

    return result;
}

BackendManager::Progress BackendManager::parse_progress(const std::string& work_dir) const {
    const std::filesystem::path progress_path = std::filesystem::path(work_dir) / "progress.json";
    if (!std::filesystem::exists(progress_path)) {
        return {};
    }

    try {
        const json j = load_json_file(progress_path);
        Progress progress;
        progress.stage = j.value("stage", std::string{});
        progress.iteration = j.value("iteration", 0);
        progress.energy = j.value("energy", 0.0);
        progress.message = j.value("message", std::string{});
        return progress;
    } catch (...) {
        return {};
    }
}

std::string BackendManager::driver_script(Method method) const {
    const char* script_name = method_is_xtb(method) ? "xtb_driver.py" : "pyscf_driver.py";
    return (std::filesystem::path(scripts_dir_) / script_name).string();
}

}  // namespace sbox::backend
