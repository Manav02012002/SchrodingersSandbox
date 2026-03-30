#include "ui/setup_wizard.h"

#include <imgui.h>

#include <chrono>
#include <cstdio>
#include <cstring>
#include <future>

namespace sbox::ui {

namespace {

using Page = SetupWizardState::Page;

template <typename T>
bool future_ready(std::future<T>& future) {
    return future.valid() &&
           future.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

void begin_detection(SetupWizardState& state, sbox::backend::PythonEnvironment& env) {
    state.detection_running = true;
    state.detection_done = false;
    state.current_page = Page::Detecting;
    state.detection_future = std::async(std::launch::async, [env]() mutable {
        env.detect();
        env.check_packages();
        return env.info();
    });
}

void begin_install(SetupWizardState& state, sbox::backend::PythonEnvironment& env) {
    state.installing = true;
    state.install_complete = false;
    state.install_success = false;
    state.install_log.clear();
    state.current_page = Page::Installing;
    const std::string python_path = env.info().python_path;
    state.install_future = std::async(std::launch::async, [python_path]() {
        std::string output;
        const std::string command = "'" + python_path + "' -m pip install pyscf 2>&1";
        const int status = sbox::backend::PythonEnvironment::run_capture(command, output);
        return std::make_pair(status, output);
    });
}

void show_package_line(const char* name, bool available, const std::string& version) {
    if (available) {
        ImGui::Text("%s: %s [ok]", name, version.empty() ? "installed" : version.c_str());
    } else {
        ImGui::Text("%s: not installed [missing]", name);
    }
}

}  // namespace

bool draw_setup_wizard(SetupWizardState& state, sbox::backend::PythonEnvironment& env) {
    if (!state.show_wizard) {
        return false;
    }

    if (!ImGui::IsPopupOpen("Python Setup")) {
        ImGui::OpenPopup("Python Setup");
    }

    if (state.current_page == Page::Detecting && future_ready(state.detection_future)) {
        state.detected_info = state.detection_future.get();
        state.detection_done = true;
        state.detection_running = false;
        if (state.detected_info.valid) {
            env.set_python_path(state.detected_info.python_path);
        }
        state.current_page = (env.is_valid() ? Page::Found : Page::NotFound);
    }

    if (state.current_page == Page::Installing && future_ready(state.install_future)) {
        const auto [status, log] = state.install_future.get();
        state.install_log = log;
        env.check_packages();
        state.detected_info = env.info();
        state.install_success = (status == 0) && env.has_pyscf();
        state.install_complete = true;
        state.installing = false;
        state.current_page = state.install_success ? Page::Done : Page::Found;
        if (!state.install_success && state.install_log.empty()) {
            state.install_log = "PySCF installation failed.";
        }
    }

    ImGui::SetNextWindowSize(ImVec2(620.0f, 360.0f), ImGuiCond_Appearing);
    const bool opened = ImGui::BeginPopupModal("Python Setup", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
    if (!opened) {
        return true;
    }

    switch (state.current_page) {
    case Page::Welcome:
        ImGui::TextWrapped("Schrodinger's Sandbox needs Python with PySCF to run electronic structure calculations.");
        ImGui::Spacing();
        ImGui::TextWrapped("The atomic orbital viewer works without Python - only molecular calculations require it.");
        ImGui::Spacing();
        if (ImGui::Button("Detect Automatically")) {
            begin_detection(state, env);
        }
        ImGui::SameLine();
        if (ImGui::Button("Enter Path Manually")) {
            state.current_page = Page::ManualPath;
            std::snprintf(state.python_path_input, sizeof(state.python_path_input), "%s", env.info().python_path.c_str());
        }
        ImGui::SameLine();
        if (ImGui::Button("Skip for Now")) {
            state.show_wizard = false;
            ImGui::CloseCurrentPopup();
        }
        break;

    case Page::Detecting:
        ImGui::TextUnformatted("Searching for Python...");
        ImGui::Spacing();
        ImGui::TextDisabled("Checking SBOX_PYTHON, saved preferences, PATH, and common install locations.");
        ImGui::TextUnformatted("Please wait...");
        break;

    case Page::NotFound:
        ImGui::TextWrapped("Could not find a suitable Python environment.");
        ImGui::Spacing();
        ImGui::TextWrapped("You need Python 3.8+ with PySCF installed.");
        ImGui::Spacing();
        ImGui::TextUnformatted("Quick setup (recommended):");
        ImGui::TextUnformatted("  pip install pyscf");
        ImGui::Spacing();
        ImGui::TextUnformatted("Or with conda:");
        ImGui::TextUnformatted("  conda install -c conda-forge pyscf");
        ImGui::Spacing();
        if (ImGui::Button("Try Again")) {
            begin_detection(state, env);
        }
        ImGui::SameLine();
        if (ImGui::Button("Enter Path Manually")) {
            state.current_page = Page::ManualPath;
        }
        ImGui::SameLine();
        if (ImGui::Button("Skip")) {
            state.show_wizard = false;
            ImGui::CloseCurrentPopup();
        }
        break;

    case Page::Found:
        ImGui::TextWrapped("Found Python at: %s", env.info().python_path.c_str());
        ImGui::Text("Version: %s", env.info().version.c_str());
        show_package_line("PySCF", env.info().has_pyscf, env.info().pyscf_version);
        show_package_line("tblite", env.info().has_tblite, env.info().tblite_version);
        ImGui::Spacing();
        if (env.info().has_pyscf) {
            if (ImGui::Button("Use This Python")) {
                env.save_preference();
                state.show_wizard = false;
                ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            if (ImGui::Button("Search Again")) {
                begin_detection(state, env);
            }
        } else {
            ImGui::TextWrapped("PySCF is not installed in this Python environment.");
            if (ImGui::Button("Install PySCF")) {
                begin_install(state, env);
            }
            ImGui::SameLine();
            if (ImGui::Button("Enter Different Path")) {
                state.current_page = Page::ManualPath;
            }
            ImGui::SameLine();
            if (ImGui::Button("Skip")) {
                state.show_wizard = false;
                ImGui::CloseCurrentPopup();
            }
        }
        break;

    case Page::ManualPath:
        ImGui::TextWrapped("Enter the full path to a Python executable.");
        ImGui::InputText("Python Path", state.python_path_input, sizeof(state.python_path_input));
        if (ImGui::Button("Test This Path")) {
            env.set_python_path(state.python_path_input);
            env.check_packages();
            state.detected_info = env.info();
            state.current_page = env.is_valid() ? Page::Found : Page::NotFound;
        }
        ImGui::SameLine();
        if (ImGui::Button("Back")) {
            state.current_page = Page::Welcome;
        }
        break;

    case Page::Installing:
        ImGui::TextUnformatted("Installing PySCF...");
        ImGui::Spacing();
        ImGui::BeginChild("InstallLog", ImVec2(560.0f, 220.0f), true, ImGuiWindowFlags_HorizontalScrollbar);
        ImGui::TextUnformatted(state.install_log.empty() ? "Running pip install pyscf..." : state.install_log.c_str());
        ImGui::EndChild();
        break;

    case Page::Done:
        ImGui::TextUnformatted("Setup complete! PySCF is ready.");
        ImGui::Spacing();
        if (ImGui::Button("Start Using Schrodinger's Sandbox")) {
            env.save_preference();
            state.show_wizard = false;
            ImGui::CloseCurrentPopup();
        }
        break;
    }

    ImGui::EndPopup();
    return true;
}

}  // namespace sbox::ui
