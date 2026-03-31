#include "core/settings.h"

#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <string>

namespace {

class ScopedAppDataDir {
public:
    explicit ScopedAppDataDir(const std::string& path) : path_(path) {
        std::filesystem::create_directories(path_);
        const char* old = std::getenv("SBOX_APP_DATA_DIR");
        if (old != nullptr) {
            had_old_ = true;
            old_value_ = old;
        }
#ifdef _WIN32
        _putenv_s("SBOX_APP_DATA_DIR", path_.c_str());
#else
        setenv("SBOX_APP_DATA_DIR", path_.c_str(), 1);
#endif
    }

    ~ScopedAppDataDir() {
#ifdef _WIN32
        if (had_old_) {
            _putenv_s("SBOX_APP_DATA_DIR", old_value_.c_str());
        } else {
            _putenv_s("SBOX_APP_DATA_DIR", "");
        }
#else
        if (had_old_) {
            setenv("SBOX_APP_DATA_DIR", old_value_.c_str(), 1);
        } else {
            unsetenv("SBOX_APP_DATA_DIR");
        }
#endif
    }

private:
    std::string path_;
    bool had_old_ = false;
    std::string old_value_;
};

std::string temp_dir(const char* stem) {
    static int counter = 0;
    const auto path = std::filesystem::temp_directory_path() /
                      ("sbox_settings_" + std::string(stem) + "_" + std::to_string(counter++));
    std::error_code ec;
    std::filesystem::remove_all(path, ec);
    std::filesystem::create_directories(path);
    return path.string();
}

TEST(Settings, RoundTripSerialization) {
    sbox::Settings settings;
    settings.window_width = 1920;
    settings.window_height = 1080;
    settings.window_maximized = true;
    settings.volume_steps = 321;
    settings.isosurface_steps = 654;
    settings.default_iso_value = 0.25f;
    settings.default_gamma = 0.8f;
    settings.mol_render_mode = 2;
    settings.color_mode = 3;
    settings.show_atom_labels = true;
    settings.show_hydrogen_labels = true;
    settings.atom_scale = 1.2f;
    settings.bond_scale = 0.9f;
    settings.default_method = 4;
    settings.default_basis = 7;
    settings.default_charge = -1;
    settings.default_multiplicity = 3;
    settings.max_scf_cycles = 333;
    settings.scf_convergence = 1.0e-10;
    settings.cube_resolution = 99;
    settings.auto_optimize_xTB = false;
    settings.python_path = "/usr/bin/python3";
    settings.python_auto_detect = false;
    settings.font_size = 17.5f;
    settings.dark_mode = false;
    settings.ui_scale = 1.5f;
    settings.show_status_bar = false;
    settings.show_fps = false;
    settings.last_open_directory = "/tmp/open";
    settings.last_save_directory = "/tmp/save";
    settings.check_for_updates = false;
    settings.update_notify_only = false;
    settings.skipped_update_version = "1.2.3";
    settings.lod_threshold_atoms = 512;
    settings.enable_antialiasing = false;
    settings.enable_vsync = false;

    const auto json = settings.to_json();
    const sbox::Settings loaded = sbox::Settings::from_json(json);

    EXPECT_EQ(loaded.window_width, settings.window_width);
    EXPECT_EQ(loaded.window_height, settings.window_height);
    EXPECT_EQ(loaded.window_maximized, settings.window_maximized);
    EXPECT_EQ(loaded.volume_steps, settings.volume_steps);
    EXPECT_EQ(loaded.isosurface_steps, settings.isosurface_steps);
    EXPECT_FLOAT_EQ(loaded.default_iso_value, settings.default_iso_value);
    EXPECT_FLOAT_EQ(loaded.default_gamma, settings.default_gamma);
    EXPECT_EQ(loaded.mol_render_mode, settings.mol_render_mode);
    EXPECT_EQ(loaded.color_mode, settings.color_mode);
    EXPECT_EQ(loaded.show_atom_labels, settings.show_atom_labels);
    EXPECT_EQ(loaded.show_hydrogen_labels, settings.show_hydrogen_labels);
    EXPECT_FLOAT_EQ(loaded.atom_scale, settings.atom_scale);
    EXPECT_FLOAT_EQ(loaded.bond_scale, settings.bond_scale);
    EXPECT_EQ(loaded.default_method, settings.default_method);
    EXPECT_EQ(loaded.default_basis, settings.default_basis);
    EXPECT_EQ(loaded.default_charge, settings.default_charge);
    EXPECT_EQ(loaded.default_multiplicity, settings.default_multiplicity);
    EXPECT_EQ(loaded.max_scf_cycles, settings.max_scf_cycles);
    EXPECT_DOUBLE_EQ(loaded.scf_convergence, settings.scf_convergence);
    EXPECT_EQ(loaded.cube_resolution, settings.cube_resolution);
    EXPECT_EQ(loaded.auto_optimize_xTB, settings.auto_optimize_xTB);
    EXPECT_EQ(loaded.python_path, settings.python_path);
    EXPECT_EQ(loaded.python_auto_detect, settings.python_auto_detect);
    EXPECT_FLOAT_EQ(loaded.font_size, settings.font_size);
    EXPECT_EQ(loaded.dark_mode, settings.dark_mode);
    EXPECT_FLOAT_EQ(loaded.ui_scale, settings.ui_scale);
    EXPECT_EQ(loaded.show_status_bar, settings.show_status_bar);
    EXPECT_EQ(loaded.show_fps, settings.show_fps);
    EXPECT_EQ(loaded.last_open_directory, settings.last_open_directory);
    EXPECT_EQ(loaded.last_save_directory, settings.last_save_directory);
    EXPECT_EQ(loaded.check_for_updates, settings.check_for_updates);
    EXPECT_EQ(loaded.update_notify_only, settings.update_notify_only);
    EXPECT_EQ(loaded.skipped_update_version, settings.skipped_update_version);
    EXPECT_EQ(loaded.lod_threshold_atoms, settings.lod_threshold_atoms);
    EXPECT_EQ(loaded.enable_antialiasing, settings.enable_antialiasing);
    EXPECT_EQ(loaded.enable_vsync, settings.enable_vsync);
}

TEST(Settings, MissingKeysUseDefaults) {
    const nlohmann::json j = {{"window_width", 800}};
    const sbox::Settings loaded = sbox::Settings::from_json(j);
    EXPECT_EQ(loaded.window_width, 800);
    EXPECT_EQ(loaded.window_height, 1000);
    EXPECT_EQ(loaded.default_method, 0);
    EXPECT_TRUE(loaded.enable_vsync);
}

TEST(Settings, ExtraKeysAreIgnored) {
    const nlohmann::json j = {{"window_width", 900}, {"unknown_key", 123}, {"another", "value"}};
    EXPECT_NO_THROW({
        const sbox::Settings loaded = sbox::Settings::from_json(j);
        EXPECT_EQ(loaded.window_width, 900);
    });
}

TEST(Settings, ManagerSavesAndLoads) {
    const ScopedAppDataDir scoped(temp_dir("manager"));
    sbox::SettingsManager manager;
    manager.load();
    manager.settings().window_width = 1234;
    manager.settings().python_path = "/opt/python";
    manager.save();

    sbox::SettingsManager reloaded;
    reloaded.load();
    EXPECT_EQ(reloaded.settings().window_width, 1234);
    EXPECT_EQ(reloaded.settings().python_path, "/opt/python");
}

TEST(Settings, RecentFilesBehavior) {
    const ScopedAppDataDir scoped(temp_dir("recent"));
    sbox::SettingsManager manager;
    manager.load();

    for (int i = 0; i < 5; ++i) {
        manager.add_recent_file("/tmp/file" + std::to_string(i));
    }
    ASSERT_EQ(manager.recent_files().size(), 5u);
    EXPECT_EQ(manager.recent_files().front(), "/tmp/file4");

    manager.add_recent_file("/tmp/file2");
    ASSERT_EQ(manager.recent_files().size(), 5u);
    EXPECT_EQ(manager.recent_files().front(), "/tmp/file2");

    for (int i = 0; i < 25; ++i) {
        manager.add_recent_file("/tmp/many" + std::to_string(i));
    }
    EXPECT_EQ(manager.recent_files().size(), 20u);
    EXPECT_EQ(manager.recent_files().front(), "/tmp/many24");
}

}  // namespace
