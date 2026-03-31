#include "core/settings.h"
#include "core/logging.h"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace sbox {

namespace {

using json = nlohmann::json;

template <typename T>
void load_if_present(const json& j, const char* key, T& value) {
    const auto it = j.find(key);
    if (it != j.end() && !it->is_null()) {
        value = it->get<T>();
    }
}

std::filesystem::path app_data_base() {
    if (const char* override_dir = std::getenv("SBOX_APP_DATA_DIR")) {
        return std::filesystem::path(override_dir);
    }
#ifdef _WIN32
    if (const char* appdata = std::getenv("APPDATA")) {
        return std::filesystem::path(appdata) / "SchrodingersSandbox";
    }
    if (const char* userprofile = std::getenv("USERPROFILE")) {
        return std::filesystem::path(userprofile) / "AppData/Roaming/SchrodingersSandbox";
    }
#elif defined(__APPLE__)
    if (const char* home = std::getenv("HOME")) {
        return std::filesystem::path(home) / "Library/Application Support/SchrodingersSandbox";
    }
#else
    if (const char* xdg = std::getenv("XDG_CONFIG_HOME")) {
        return std::filesystem::path(xdg) / "SchrodingersSandbox";
    }
    if (const char* home = std::getenv("HOME")) {
        return std::filesystem::path(home) / ".config/SchrodingersSandbox";
    }
#endif
    return std::filesystem::temp_directory_path() / "SchrodingersSandbox";
}

std::filesystem::path settings_path_fs() {
    return std::filesystem::path(get_app_data_dir()) / "settings.json";
}

std::filesystem::path recent_files_path_fs() {
    return std::filesystem::path(get_app_data_dir()) / "recent_files.json";
}

}  // namespace

std::string get_app_data_dir() {
    const std::filesystem::path dir = app_data_base();
    std::filesystem::create_directories(dir);
    return dir.string() + std::string(1, std::filesystem::path::preferred_separator);
}

std::string get_settings_path() {
    return settings_path_fs().string();
}

std::string get_recent_files_path() {
    return recent_files_path_fs().string();
}

nlohmann::json Settings::to_json() const {
    return json{
        {"window_width", window_width},
        {"window_height", window_height},
        {"window_maximized", window_maximized},
        {"volume_steps", volume_steps},
        {"isosurface_steps", isosurface_steps},
        {"default_iso_value", default_iso_value},
        {"default_gamma", default_gamma},
        {"mol_render_mode", mol_render_mode},
        {"color_mode", color_mode},
        {"show_atom_labels", show_atom_labels},
        {"show_hydrogen_labels", show_hydrogen_labels},
        {"atom_scale", atom_scale},
        {"bond_scale", bond_scale},
        {"default_method", default_method},
        {"default_basis", default_basis},
        {"default_charge", default_charge},
        {"default_multiplicity", default_multiplicity},
        {"max_scf_cycles", max_scf_cycles},
        {"scf_convergence", scf_convergence},
        {"cube_resolution", cube_resolution},
        {"auto_optimize_xTB", auto_optimize_xTB},
        {"python_path", python_path},
        {"python_auto_detect", python_auto_detect},
        {"font_size", font_size},
        {"dark_mode", dark_mode},
        {"ui_scale", ui_scale},
        {"show_status_bar", show_status_bar},
        {"show_fps", show_fps},
        {"last_open_directory", last_open_directory},
        {"last_save_directory", last_save_directory},
        {"check_for_updates", check_for_updates},
        {"update_notify_only", update_notify_only},
        {"skipped_update_version", skipped_update_version},
        {"lod_threshold_atoms", lod_threshold_atoms},
        {"enable_antialiasing", enable_antialiasing},
        {"enable_vsync", enable_vsync},
    };
}

Settings Settings::from_json(const nlohmann::json& j) {
    Settings settings;
    if (!j.is_object()) {
        return settings;
    }

    load_if_present(j, "window_width", settings.window_width);
    load_if_present(j, "window_height", settings.window_height);
    load_if_present(j, "window_maximized", settings.window_maximized);
    load_if_present(j, "volume_steps", settings.volume_steps);
    load_if_present(j, "isosurface_steps", settings.isosurface_steps);
    load_if_present(j, "default_iso_value", settings.default_iso_value);
    load_if_present(j, "default_gamma", settings.default_gamma);
    load_if_present(j, "mol_render_mode", settings.mol_render_mode);
    load_if_present(j, "color_mode", settings.color_mode);
    load_if_present(j, "show_atom_labels", settings.show_atom_labels);
    load_if_present(j, "show_hydrogen_labels", settings.show_hydrogen_labels);
    load_if_present(j, "atom_scale", settings.atom_scale);
    load_if_present(j, "bond_scale", settings.bond_scale);
    load_if_present(j, "default_method", settings.default_method);
    load_if_present(j, "default_basis", settings.default_basis);
    load_if_present(j, "default_charge", settings.default_charge);
    load_if_present(j, "default_multiplicity", settings.default_multiplicity);
    load_if_present(j, "max_scf_cycles", settings.max_scf_cycles);
    load_if_present(j, "scf_convergence", settings.scf_convergence);
    load_if_present(j, "cube_resolution", settings.cube_resolution);
    load_if_present(j, "auto_optimize_xTB", settings.auto_optimize_xTB);
    load_if_present(j, "python_path", settings.python_path);
    load_if_present(j, "python_auto_detect", settings.python_auto_detect);
    load_if_present(j, "font_size", settings.font_size);
    load_if_present(j, "dark_mode", settings.dark_mode);
    load_if_present(j, "ui_scale", settings.ui_scale);
    load_if_present(j, "show_status_bar", settings.show_status_bar);
    load_if_present(j, "show_fps", settings.show_fps);
    load_if_present(j, "last_open_directory", settings.last_open_directory);
    load_if_present(j, "last_save_directory", settings.last_save_directory);
    load_if_present(j, "check_for_updates", settings.check_for_updates);
    load_if_present(j, "update_notify_only", settings.update_notify_only);
    load_if_present(j, "skipped_update_version", settings.skipped_update_version);
    load_if_present(j, "lod_threshold_atoms", settings.lod_threshold_atoms);
    load_if_present(j, "enable_antialiasing", settings.enable_antialiasing);
    load_if_present(j, "enable_vsync", settings.enable_vsync);
    return settings;
}

SettingsManager::SettingsManager() = default;

void SettingsManager::load() {
    settings_ = Settings{};
    recent_files_.clear();

    const auto settings_path = settings_path_fs();
    if (std::filesystem::exists(settings_path)) {
        try {
            std::ifstream in(settings_path);
            json j;
            in >> j;
            settings_ = Settings::from_json(j);
        } catch (const std::exception& ex) {
            SBOX_LOG_WARN("Failed to load settings from %s: %s", settings_path.string().c_str(), ex.what());
            settings_ = Settings{};
        }
    }

    const auto recent_path = recent_files_path_fs();
    if (std::filesystem::exists(recent_path)) {
        try {
            std::ifstream in(recent_path);
            json j;
            in >> j;
            if (j.is_array()) {
                for (const auto& item : j) {
                    if (item.is_string()) {
                        recent_files_.push_back(item.get<std::string>());
                    }
                }
            }
        } catch (const std::exception& ex) {
            SBOX_LOG_WARN("Failed to load recent files from %s: %s", recent_path.string().c_str(), ex.what());
            recent_files_.clear();
        }
    }

    if (recent_files_.size() > 20) {
        recent_files_.resize(20);
    }
    loaded_ = true;
}

void SettingsManager::save() {
    std::filesystem::create_directories(app_data_base());

    const auto settings_path = settings_path_fs();
    std::ofstream settings_out(settings_path);
    if (!settings_out) {
        throw std::runtime_error("Could not open settings file for writing: " + settings_path.string());
    }
    settings_out << settings_.to_json().dump(2);

    const auto recent_path = recent_files_path_fs();
    std::ofstream recent_out(recent_path);
    if (!recent_out) {
        throw std::runtime_error("Could not open recent files file for writing: " + recent_path.string());
    }
    recent_out << json(recent_files_).dump(2);
}

Settings& SettingsManager::settings() {
    return settings_;
}

const Settings& SettingsManager::settings() const {
    return settings_;
}

void SettingsManager::add_recent_file(const std::string& path) {
    recent_files_.erase(std::remove(recent_files_.begin(), recent_files_.end(), path), recent_files_.end());
    recent_files_.insert(recent_files_.begin(), path);
    if (recent_files_.size() > 20) {
        recent_files_.resize(20);
    }
    if (!loaded_) {
        loaded_ = true;
    }
    save();
}

const std::vector<std::string>& SettingsManager::recent_files() const {
    return recent_files_;
}

void SettingsManager::clear_recent_files() {
    recent_files_.clear();
    if (!loaded_) {
        loaded_ = true;
    }
    save();
}

}  // namespace sbox
