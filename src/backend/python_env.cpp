#include "backend/python_env.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace sbox::backend {

namespace {

std::string trim(const std::string& input) {
    const std::size_t start = input.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
        return "";
    }
    const std::size_t end = input.find_last_not_of(" \t\r\n");
    return input.substr(start, end - start + 1);
}

std::string shell_quote(const std::string& value) {
    std::string quoted = "'";
    for (char ch : value) {
        if (ch == '\'') {
            quoted += "'\\''";
        } else {
            quoted.push_back(ch);
        }
    }
    quoted.push_back('\'');
    return quoted;
}

std::filesystem::path home_directory() {
    if (const char* home = std::getenv("HOME")) {
        return std::filesystem::path(home);
    }
    return {};
}

std::string python_version_command(const std::string& path) {
    return shell_quote(path)
        + " -c \"import sys; print(f'{sys.version_info[0]}.{sys.version_info[1]}.{sys.version_info[2]}')\"";
}

std::string resolve_python_path(const std::string& path) {
    if (path.find('/') != std::string::npos) {
        return path;
    }

    std::string resolved;
    const int status = PythonEnvironment::run_capture("command -v " + path, resolved);
    if (status != 0) {
        return path;
    }
    const std::string trimmed = trim(resolved);
    return trimmed.empty() ? path : trimmed;
}

bool probe_package(const std::string& python_path,
                   const std::string& import_name,
                   std::string* out_version,
                   const std::string& version_expr = "__version__") {
    std::string output;
    const std::string command = shell_quote(python_path) + " -c \"import " + import_name
        + "; print(getattr(" + import_name + ", '" + version_expr + "', 'unknown'))\"";
    const int status = PythonEnvironment::run_capture(command, output);
    if (status != 0) {
        return false;
    }
    if (out_version != nullptr) {
        *out_version = trim(output);
    }
    return true;
}

}  // namespace

PythonEnvironment::PythonEnvironment() = default;

void PythonEnvironment::detect() {
    info_ = PythonInfo{};

    std::vector<std::string> candidates;
    if (const char* env_python = std::getenv("SBOX_PYTHON")) {
        candidates.emplace_back(env_python);
    }

    const std::filesystem::path pref_path = home_directory() / ".sbox_python";
    if (std::filesystem::exists(pref_path)) {
        std::ifstream in(pref_path);
        std::string saved;
        std::getline(in, saved);
        if (!trim(saved).empty()) {
            candidates.push_back(trim(saved));
        }
    }

    candidates.push_back("python3");
    candidates.push_back("python");

    const std::filesystem::path home = home_directory();
    candidates.push_back("/usr/local/bin/python3");
    candidates.push_back("/opt/homebrew/bin/python3");
    if (!home.empty()) {
        candidates.push_back((home / "miniconda3/bin/python3").string());
        candidates.push_back((home / "anaconda3/bin/python3").string());
    }

    std::vector<std::string> deduped;
    for (const std::string& candidate : candidates) {
        if (candidate.empty()) {
            continue;
        }
        if (std::find(deduped.begin(), deduped.end(), candidate) == deduped.end()) {
            deduped.push_back(candidate);
        }
    }

    for (const std::string& candidate : deduped) {
        if (test_python(candidate)) {
            check_packages();
            return;
        }
    }
}

void PythonEnvironment::check_packages() {
    if (!info_.valid) {
        return;
    }

    info_.has_pyscf = probe_package(info_.python_path, "pyscf", &info_.pyscf_version);
    info_.has_tblite = probe_package(info_.python_path, "tblite", &info_.tblite_version);
    info_.has_xtb = probe_package(info_.python_path, "xtb", &info_.xtb_version);
    info_.has_geometric = probe_package(info_.python_path, "geometric", nullptr);
    info_.has_ase = probe_package(info_.python_path, "ase", nullptr);
}

const PythonInfo& PythonEnvironment::info() const {
    return info_;
}

bool PythonEnvironment::is_valid() const {
    return info_.valid;
}

bool PythonEnvironment::has_pyscf() const {
    return info_.has_pyscf;
}

bool PythonEnvironment::has_tblite() const {
    return info_.has_tblite;
}

void PythonEnvironment::save_preference() const {
    if (!info_.valid || info_.python_path.empty()) {
        return;
    }

    const std::filesystem::path pref_path = home_directory() / ".sbox_python";
    std::ofstream out(pref_path);
    if (out) {
        out << info_.python_path << '\n';
    }
}

void PythonEnvironment::set_python_path(const std::string& path) {
    info_ = PythonInfo{};
    (void)test_python(path);
    if (info_.valid) {
        check_packages();
    }
}

int PythonEnvironment::run_capture(const std::string& command, std::string& stdout_out) {
    stdout_out.clear();

    // TODO: replace popen with a timeout-capable fork/exec implementation.
    FILE* pipe = popen(command.c_str(), "r");
    if (pipe == nullptr) {
        return -1;
    }

    std::array<char, 256> buffer{};
    while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) != nullptr) {
        stdout_out += buffer.data();
    }

    const int status = pclose(pipe);
    if (status == -1) {
        return -1;
    }
    return status;
}

bool PythonEnvironment::test_python(const std::string& path) {
    const std::string resolved_path = resolve_python_path(path);
    std::string output;
    const int status = run_capture(python_version_command(resolved_path), output);
    if (status != 0) {
        return false;
    }

    const std::string version = trim(output);
    if (version.empty()) {
        return false;
    }

    info_.python_path = resolved_path;
    info_.version = version;
    info_.valid = true;
    return true;
}

}  // namespace sbox::backend
