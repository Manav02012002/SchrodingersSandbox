#include "core/paths.h"

#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#else
#include <unistd.h>
#endif

namespace sbox {

namespace {

std::filesystem::path executable_path() {
#ifdef _WIN32
    std::wstring buffer(MAX_PATH, L'\0');
    DWORD len = GetModuleFileNameW(nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
    while (len == buffer.size()) {
        buffer.resize(buffer.size() * 2);
        len = GetModuleFileNameW(nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
    }
    if (len == 0) {
        throw std::runtime_error("Could not determine executable path");
    }
    buffer.resize(len);
    return std::filesystem::path(buffer);
#elif defined(__APPLE__)
    uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);
    std::string buffer(size, '\0');
    if (_NSGetExecutablePath(buffer.data(), &size) != 0) {
        throw std::runtime_error("Could not determine executable path");
    }
    return std::filesystem::weakly_canonical(std::filesystem::path(buffer.c_str()));
#else
    std::vector<char> buffer(1024, '\0');
    ssize_t len = readlink("/proc/self/exe", buffer.data(), buffer.size());
    while (len >= 0 && static_cast<std::size_t>(len) == buffer.size()) {
        buffer.resize(buffer.size() * 2);
        len = readlink("/proc/self/exe", buffer.data(), buffer.size());
    }
    if (len < 0) {
        throw std::runtime_error("Could not determine executable path");
    }
    return std::filesystem::path(std::string(buffer.data(), static_cast<std::size_t>(len)));
#endif
}

bool is_valid_data_dir(const std::filesystem::path& dir) {
    std::error_code ec;
    return std::filesystem::exists(dir / "shaders", ec) && std::filesystem::is_directory(dir / "shaders", ec);
}

std::filesystem::path normalized(const std::filesystem::path& path) {
    std::error_code ec;
    const auto canonical = std::filesystem::weakly_canonical(path, ec);
    return ec ? path.lexically_normal() : canonical;
}

}  // namespace

std::string get_data_dir() {
    std::vector<std::filesystem::path> candidates;

    if (const char* env_dir = std::getenv("SBOX_DATA_DIR")) {
        if (*env_dir != '\0') {
            candidates.emplace_back(env_dir);
        }
    }

    const std::filesystem::path exe_dir = executable_path().parent_path();
    candidates.push_back(std::filesystem::current_path() / "data");
    candidates.push_back(exe_dir / "data");
    candidates.push_back(exe_dir / "../Resources/data");
    candidates.push_back(exe_dir / "../share/schrodingers_sandbox/data");
    candidates.push_back(exe_dir / "data");

    for (const auto& candidate : candidates) {
        const std::filesystem::path dir = normalized(candidate);
        if (is_valid_data_dir(dir)) {
            return dir.string();
        }
    }

    std::string message = "Could not locate application data directory. Checked:";
    for (const auto& candidate : candidates) {
        message += "\n  - " + normalized(candidate).string();
    }
    message += "\nSet SBOX_DATA_DIR to override the search path.";
    throw std::runtime_error(message);
}

std::string get_shader_path(const std::string& shader_name) {
    return (std::filesystem::path(get_data_dir()) / "shaders" / shader_name).string();
}

std::string get_script_path(const std::string& script_name) {
    return (std::filesystem::path(get_data_dir()) / "scripts" / script_name).string();
}

}  // namespace sbox
