#pragma once

#include <json.hpp>

#include <string>
#include <vector>

namespace sbox {

std::string get_app_data_dir();
std::string get_settings_path();
std::string get_recent_files_path();

struct Settings {
    int window_width = 1600;
    int window_height = 1000;
    bool window_maximized = false;

    int volume_steps = 192;
    int isosurface_steps = 256;
    float default_iso_value = 0.01f;
    float default_gamma = 0.4f;
    int mol_render_mode = 0;
    int color_mode = 0;
    bool show_atom_labels = false;
    bool show_hydrogen_labels = false;
    float atom_scale = 1.0f;
    float bond_scale = 1.0f;

    int default_method = 0;
    int default_basis = 0;
    int default_charge = 0;
    int default_multiplicity = 1;
    int max_scf_cycles = 200;
    double scf_convergence = 1e-8;
    int cube_resolution = 80;
    bool auto_optimize_xTB = true;

    std::string python_path;
    bool python_auto_detect = true;

    float font_size = 15.0f;
    bool dark_mode = true;
    float ui_scale = 1.0f;
    bool show_status_bar = true;
    bool show_fps = true;

    std::string last_open_directory;
    std::string last_save_directory;

    bool check_for_updates = true;
    bool update_notify_only = true;
    std::string skipped_update_version;

    int lod_threshold_atoms = 200;
    bool enable_antialiasing = true;
    bool enable_vsync = true;

    nlohmann::json to_json() const;
    static Settings from_json(const nlohmann::json& j);
};

class SettingsManager {
public:
    SettingsManager();

    void load();
    void save();

    Settings& settings();
    const Settings& settings() const;

    void add_recent_file(const std::string& path);
    const std::vector<std::string>& recent_files() const;
    void clear_recent_files();

private:
    Settings settings_;
    std::vector<std::string> recent_files_;
    bool loaded_ = false;
};

}  // namespace sbox
