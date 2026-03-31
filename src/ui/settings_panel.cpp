#include "ui/settings_panel.h"

#include "backend/job_types.h"
#include "version.h"

#include <imgui.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <future>
#include <string>
#include <utility>

namespace sbox::ui {

namespace {

using sbox::backend::BasisSetType;
using sbox::backend::Method;

constexpr std::array<Method, 12> kMethods = {
    Method::HF,
    Method::UHF,
    Method::DFT_B3LYP,
    Method::DFT_PBE,
    Method::DFT_PBE0,
    Method::DFT_TPSS,
    Method::DFT_M06_2X,
    Method::MP2,
    Method::CCSD,
    Method::GFN2_XTB,
    Method::GFN1_XTB,
    Method::GFN_FF,
};

constexpr std::array<BasisSetType, 14> kBasisSets = {
    BasisSetType::STO_3G,
    BasisSetType::B3_21G,
    BasisSetType::B6_31G,
    BasisSetType::B6_31Gd,
    BasisSetType::B6_31Gdp,
    BasisSetType::B6_311Gdp,
    BasisSetType::cc_pVDZ,
    BasisSetType::cc_pVTZ,
    BasisSetType::cc_pVQZ,
    BasisSetType::def2_SVP,
    BasisSetType::def2_TZVP,
    BasisSetType::def2_TZVPP,
    BasisSetType::aug_cc_pVDZ,
    BasisSetType::aug_cc_pVTZ,
};

struct InstallResult {
    int status = -1;
    std::string output;
};

struct PanelState {
    sbox::Settings original;
    std::string python_path_input;
    std::string python_test_result;
    std::string update_status = "Not checked.";
    std::string install_log;
    std::future<InstallResult> install_future;
    bool install_running = false;
    bool request_close = false;
};

PanelState& panel_state() {
    static PanelState state;
    return state;
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

void help_marker(const char* text) {
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 28.0f);
        ImGui::TextUnformatted(text);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

void sync_python_path_buffer(PanelState& panel, const sbox::Settings& settings) {
    if (panel.python_path_input.empty()) {
        panel.python_path_input = settings.python_path;
    }
}

void draw_package_line(const char* label, bool available, const std::string& version) {
    ImGui::TextUnformatted(label);
    ImGui::SameLine(180.0f);
    if (available) {
        ImGui::TextColored(ImVec4(0.35f, 0.85f, 0.45f, 1.0f), "%s", version.empty() ? "found" : version.c_str());
    } else {
        ImGui::TextColored(ImVec4(0.90f, 0.30f, 0.30f, 1.0f), "not found");
    }
}

}  // namespace

bool consume_settings_panel_close_request() {
    PanelState& panel = panel_state();
    const bool value = panel.request_close;
    panel.request_close = false;
    return value;
}

void draw_settings_panel(sbox::Settings& settings,
                         sbox::backend::PythonEnvironment& python_env,
                         bool& settings_changed) {
    PanelState& panel = panel_state();
    settings_changed = false;

    if (ImGui::IsWindowAppearing()) {
        panel.original = settings;
        panel.python_path_input = settings.python_path;
        panel.python_test_result.clear();
        panel.request_close = false;
    }
    sync_python_path_buffer(panel, settings);

    if (panel.install_running &&
        panel.install_future.valid() &&
        panel.install_future.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
        const InstallResult result = panel.install_future.get();
        panel.install_running = false;
        panel.install_log = result.output;
        if (result.status == 0) {
            python_env.check_packages();
            panel.python_test_result = "PySCF installation finished successfully.";
        } else {
            panel.python_test_result = "PySCF installation failed. See log below.";
        }
    }

    if (ImGui::BeginTabBar("##settings_tabs")) {
        if (ImGui::BeginTabItem("Rendering")) {
            ImGui::SliderInt("Volume Ray March Steps", &settings.volume_steps, 64, 512);
            ImGui::SameLine();
            help_marker("Higher values improve volume quality but cost performance. Default: 192.");
            ImGui::SliderInt("Isosurface Steps", &settings.isosurface_steps, 64, 512);
            ImGui::SliderFloat("Default Iso Value", &settings.default_iso_value, 0.0001f, 0.1f, "%.4f", ImGuiSliderFlags_Logarithmic);
            ImGui::SliderFloat("Default Gamma", &settings.default_gamma, 0.1f, 1.0f, "%.2f");
            ImGui::SliderFloat("Atom Size Scale", &settings.atom_scale, 0.2f, 3.0f, "%.1f");
            ImGui::SliderFloat("Bond Size Scale", &settings.bond_scale, 0.2f, 3.0f, "%.1f");
            ImGui::Checkbox("Show Atom Labels", &settings.show_atom_labels);
            ImGui::Checkbox("Show Hydrogen Labels", &settings.show_hydrogen_labels);
            ImGui::Checkbox("Anti-aliasing", &settings.enable_antialiasing);
            ImGui::Checkbox("VSync", &settings.enable_vsync);
            ImGui::SliderInt("LOD Threshold (atoms)", &settings.lod_threshold_atoms, 50, 1000);
            ImGui::SameLine();
            help_marker("Molecules above this size use simplified rendering for far atoms.");
            ImGui::EndTabItem();
        }

        if (ImGui::BeginTabItem("Computation")) {
            int method_index = 0;
            for (int i = 0; i < static_cast<int>(kMethods.size()); ++i) {
                if (static_cast<int>(kMethods[static_cast<std::size_t>(i)]) == settings.default_method) {
                    method_index = i;
                    break;
                }
            }
            if (ImGui::BeginCombo("Default Method", sbox::backend::method_display_name(kMethods[static_cast<std::size_t>(method_index)]))) {
                for (int i = 0; i < static_cast<int>(kMethods.size()); ++i) {
                    const bool selected = i == method_index;
                    if (ImGui::Selectable(sbox::backend::method_display_name(kMethods[static_cast<std::size_t>(i)]), selected)) {
                        settings.default_method = static_cast<int>(kMethods[static_cast<std::size_t>(i)]);
                    }
                    if (selected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }

            int basis_index = 0;
            for (int i = 0; i < static_cast<int>(kBasisSets.size()); ++i) {
                if (static_cast<int>(kBasisSets[static_cast<std::size_t>(i)]) == settings.default_basis) {
                    basis_index = i;
                    break;
                }
            }
            if (ImGui::BeginCombo("Default Basis Set", sbox::backend::basis_display_name(kBasisSets[static_cast<std::size_t>(basis_index)]))) {
                for (int i = 0; i < static_cast<int>(kBasisSets.size()); ++i) {
                    const bool selected = i == basis_index;
                    if (ImGui::Selectable(sbox::backend::basis_display_name(kBasisSets[static_cast<std::size_t>(i)]), selected)) {
                        settings.default_basis = static_cast<int>(kBasisSets[static_cast<std::size_t>(i)]);
                    }
                    if (selected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }

            ImGui::InputInt("Default Charge", &settings.default_charge);
            ImGui::InputInt("Default Multiplicity", &settings.default_multiplicity);
            settings.default_multiplicity = std::max(1, settings.default_multiplicity);
            ImGui::InputInt("Max SCF Cycles", &settings.max_scf_cycles);
            settings.max_scf_cycles = std::max(1, settings.max_scf_cycles);
            ImGui::InputDouble("SCF Convergence", &settings.scf_convergence, 0.0, 0.0, "%.1e");
            settings.scf_convergence = std::max(1.0e-14, settings.scf_convergence);
            ImGui::SliderInt("Cube Grid Resolution", &settings.cube_resolution, 40, 200);
            ImGui::Checkbox("Auto-optimize with xTB after building", &settings.auto_optimize_xTB);
            ImGui::EndTabItem();
        }

        if (ImGui::BeginTabItem("Python Environment")) {
            const auto& info = python_env.info();
            if (info.version.empty()) {
                ImGui::TextWrapped("Current Python: %s",
                                   info.python_path.empty() ? "(none)" : info.python_path.c_str());
            } else {
                ImGui::TextWrapped("Current Python: %s (version %s)",
                                   info.python_path.empty() ? "(none)" : info.python_path.c_str(),
                                   info.version.c_str());
            }
            draw_package_line("PySCF", info.has_pyscf, info.pyscf_version);
            draw_package_line("tblite", info.has_tblite, info.tblite_version);
            draw_package_line("xtb", info.has_xtb, info.xtb_version);
            draw_package_line("geomeTRIC", info.has_geometric, info.has_geometric ? "found" : "");

            ImGui::Separator();
            ImGui::Checkbox("Auto-detect Python", &settings.python_auto_detect);
            if (!settings.python_auto_detect) {
                std::array<char, 1024> buf{};
                std::snprintf(buf.data(), buf.size(), "%s", panel.python_path_input.c_str());
                if (ImGui::InputText("Python Path", buf.data(), buf.size())) {
                    panel.python_path_input = buf.data();
                    settings.python_path = panel.python_path_input;
                }
                if (ImGui::Button("Test")) {
                    if (panel.python_path_input.empty()) {
                        panel.python_test_result = "Enter a Python path before testing.";
                    } else {
                        python_env.set_python_path(panel.python_path_input);
                        settings.python_path = panel.python_path_input;
                        panel.python_test_result = python_env.is_valid()
                            ? "Python environment looks valid."
                            : "Python path test failed.";
                    }
                }
            }

            ImGui::SameLine();
            if (ImGui::Button("Re-detect")) {
                python_env.detect();
                if (python_env.is_valid()) {
                    settings.python_path = python_env.info().python_path;
                    panel.python_path_input = settings.python_path;
                    panel.python_test_result = "Detected Python successfully.";
                } else {
                    panel.python_test_result = "Failed to detect a usable Python installation.";
                }
            }

            ImGui::SameLine();
            if (!panel.install_running) {
                if (ImGui::Button("Install PySCF")) {
                    const std::string python_path = settings.python_auto_detect ? python_env.info().python_path : panel.python_path_input;
                    if (python_path.empty()) {
                        panel.python_test_result = "Choose or detect a Python interpreter first.";
                    } else {
                        panel.install_running = true;
                        panel.install_log.clear();
                        panel.install_future = std::async(std::launch::async, [python_path]() {
                            InstallResult result;
                            result.status = sbox::backend::PythonEnvironment::run_capture(
                                shell_quote(python_path) + " -m pip install pyscf", result.output);
                            return result;
                        });
                    }
                }
            } else {
                ImGui::BeginDisabled();
                ImGui::Button("Installing PySCF...");
                ImGui::EndDisabled();
            }

            if (!panel.python_test_result.empty()) {
                ImGui::Spacing();
                ImGui::TextWrapped("%s", panel.python_test_result.c_str());
            }
            if (!panel.install_log.empty()) {
                ImGui::Separator();
                ImGui::TextUnformatted("Install Log");
                if (ImGui::BeginChild("##install_log", ImVec2(-1.0f, 140.0f), true)) {
                    ImGui::TextUnformatted(panel.install_log.c_str());
                }
                ImGui::EndChild();
            }
            ImGui::EndTabItem();
        }

        if (ImGui::BeginTabItem("UI")) {
            ImGui::SliderFloat("Font Size", &settings.font_size, 10.0f, 24.0f, "%.0f");
            ImGui::TextDisabled("Requires restart to take effect.");
            ImGui::SliderFloat("UI Scale", &settings.ui_scale, 0.75f, 2.0f, "%.2f");
            ImGui::TextDisabled("For HiDPI displays. Requires restart.");
            ImGui::Checkbox("Show Status Bar", &settings.show_status_bar);
            ImGui::Checkbox("Show FPS Counter", &settings.show_fps);
            ImGui::EndTabItem();
        }

        if (ImGui::BeginTabItem("Updates")) {
            ImGui::Checkbox("Check for updates on startup", &settings.check_for_updates);
            ImGui::Checkbox("Notify only (don't auto-download)", &settings.update_notify_only);
            if (ImGui::Button("Check Now")) {
                panel.update_status = "Update checks are not configured for local builds.";
            }
            ImGui::Text("Current version: %s", sbox::VERSION);
            ImGui::TextWrapped("Latest version: %s", panel.update_status.c_str());
            ImGui::EndTabItem();
        }

        ImGui::EndTabBar();
    }

    if (ImGui::BeginPopupModal("##reset_settings_confirm", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::TextWrapped("Reset all preferences to their default values?");
        if (ImGui::Button("Reset", ImVec2(120.0f, 0.0f))) {
            settings = sbox::Settings{};
            panel.python_path_input = settings.python_path;
            panel.python_test_result = "Defaults restored. Click Apply or OK to save them.";
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120.0f, 0.0f))) {
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }

    ImGui::Separator();
    if (ImGui::Button("Apply", ImVec2(110.0f, 0.0f))) {
        settings.python_path = panel.python_path_input;
        settings_changed = true;
        panel.original = settings;
    }
    ImGui::SameLine();
    if (ImGui::Button("Reset to Defaults", ImVec2(150.0f, 0.0f))) {
        ImGui::OpenPopup("##reset_settings_confirm");
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel", ImVec2(110.0f, 0.0f))) {
        settings = panel.original;
        panel.python_path_input = settings.python_path;
        panel.request_close = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("OK", ImVec2(110.0f, 0.0f))) {
        settings.python_path = panel.python_path_input;
        settings_changed = true;
        panel.original = settings;
        panel.request_close = true;
    }
}

}  // namespace sbox::ui
