#include "ui/about_dialog.h"

#include "core/paths.h"
#include "core/settings.h"
#include "version.h"

#include <glad/gl.h>
#include <imgui.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace sbox::ui {

namespace {

const sbox::backend::PythonEnvironment* g_python_env = nullptr;

std::string safe_gl_string(unsigned int name) {
    const GLubyte* value = glGetString(name);
    return value != nullptr ? reinterpret_cast<const char*>(value) : "Unavailable";
}

std::string get_os_info() {
#if defined(__APPLE__)
    std::string output;
    FILE* pipe = popen("sw_vers 2>/dev/null", "r");
    if (pipe != nullptr) {
        char buffer[256];
        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            output += buffer;
        }
        pclose(pipe);
    }
    if (!output.empty()) {
        std::istringstream iss(output);
        std::string line;
        std::string product;
        std::string version;
        while (std::getline(iss, line)) {
            if (line.rfind("ProductName:", 0) == 0) {
                product = line.substr(std::string("ProductName:").size());
            } else if (line.rfind("ProductVersion:", 0) == 0) {
                version = line.substr(std::string("ProductVersion:").size());
            }
        }
        if (!product.empty() || !version.empty()) {
            return product + (product.empty() || version.empty() ? "" : " ") + version;
        }
    }
    return "macOS";
#elif defined(__linux__)
    std::ifstream in("/etc/os-release");
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("PRETTY_NAME=", 0) == 0) {
            std::string value = line.substr(std::string("PRETTY_NAME=").size());
            if (!value.empty() && value.front() == '"' && value.back() == '"') {
                value = value.substr(1, value.size() - 2);
            }
            return value;
        }
    }
    return "Linux";
#elif defined(_WIN32)
    return "Windows";
#else
    return "Unknown";
#endif
}

std::string get_data_dir_safe() {
    try {
        return sbox::get_data_dir();
    } catch (...) {
        return "Unavailable";
    }
}

std::string package_line(const char* name, bool available, const std::string& version = {}) {
    if (!available) {
        return std::string(name) + ": not installed";
    }
    if (!version.empty()) {
        return std::string(name) + ": " + version;
    }
    return std::string(name) + ": available";
}

std::string build_system_info_text() {
    std::ostringstream oss;
    oss << "Application: " << sbox::APP_NAME << "\n";
    oss << "Version: " << sbox::VERSION << "\n";
    oss << "Built: " << __DATE__ << " " << __TIME__ << "\n";
    oss << "OS: " << get_os_info() << "\n";
    oss << "GPU: " << safe_gl_string(GL_RENDERER) << "\n";
    oss << "OpenGL: " << safe_gl_string(GL_VERSION) << "\n";
    oss << "Vendor: " << safe_gl_string(GL_VENDOR) << "\n";
    if (g_python_env != nullptr) {
        const auto& info = g_python_env->info();
        oss << "Python: " << (info.valid ? info.python_path + " (" + info.version + ")" : std::string("not detected")) << "\n";
        oss << package_line("PySCF", info.has_pyscf, info.pyscf_version) << "\n";
        oss << package_line("tblite", info.has_tblite, info.tblite_version) << "\n";
        oss << package_line("xtb", info.has_xtb, info.xtb_version) << "\n";
        oss << package_line("geomeTRIC", info.has_geometric) << "\n";
    }
    oss << "Data directory: " << get_data_dir_safe() << "\n";
    oss << "Settings directory: " << sbox::get_app_data_dir() << "\n";
    oss << "Log file: " << (std::filesystem::path(sbox::get_app_data_dir()) / "schrodingers_sandbox.log").string() << "\n";
    return oss.str();
}

void draw_shortcuts_table_section(const char* title, const std::initializer_list<std::pair<const char*, const char*>>& items) {
    ImGui::SeparatorText(title);
    if (ImGui::BeginTable(title, 2, ImGuiTableFlags_SizingStretchProp)) {
        for (const auto& [shortcut, description] : items) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::TextUnformatted(shortcut);
            ImGui::TableSetColumnIndex(1);
            ImGui::TextUnformatted(description);
        }
        ImGui::EndTable();
    }
}

}  // namespace

void set_about_dialog_context(const sbox::backend::PythonEnvironment* python_env) {
    g_python_env = python_env;
}

void draw_about_dialog(bool& show) {
    if (show) {
        ImGui::OpenPopup("About Schrödinger's Sandbox");
    }

    bool open = show;
    if (ImGui::BeginPopupModal("About Schrödinger's Sandbox", &open, ImGuiWindowFlags_AlwaysAutoResize)) {
        show = open;

        const float avail = ImGui::GetContentRegionAvail().x;
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + std::max(0.0f, (avail - 220.0f) * 0.5f));
        ImGui::SetWindowFontScale(1.5f);
        ImGui::TextUnformatted("Schrödinger's Sandbox");
        ImGui::SetWindowFontScale(1.0f);

        ImGui::Separator();
        ImGui::Text("Version: %s", sbox::VERSION);
        ImGui::Text("Built: %s %s", __DATE__, __TIME__);
        ImGui::Spacing();
        ImGui::TextUnformatted("Interactive quantum chemistry visualisation");
        ImGui::TextUnformatted("and molecular orbital exploration.");
        ImGui::Spacing();
        ImGui::TextUnformatted("Author: Manav Rawal");
        ImGui::TextUnformatted("License: See repository LICENSE");
        ImGui::Spacing();

        const std::string system_info = build_system_info_text();
        if (ImGui::TreeNode("System Information")) {
            ImGui::TextWrapped("%s", system_info.c_str());
            ImGui::TreePop();
        }

        if (ImGui::Button("Copy System Info")) {
            ImGui::SetClipboardText(system_info.c_str());
        }

        if (ImGui::TreeNode("Open Source Libraries")) {
            ImGui::TextUnformatted("GLFW 3.4 - zlib/libpng license");
            ImGui::TextUnformatted("Dear ImGui - MIT license");
            ImGui::TextUnformatted("ImPlot - MIT license");
            ImGui::TextUnformatted("Eigen - MPL2 license");
            ImGui::TextUnformatted("GoogleTest - BSD 3-Clause");
            ImGui::TextUnformatted("nlohmann/json - MIT license");
            ImGui::TextUnformatted("nativefiledialog-extended - zlib license");
            ImGui::TextUnformatted("PySCF - Apache 2.0");
            ImGui::TreePop();
        }

        ImGui::Spacing();
        if (ImGui::Button("OK", ImVec2(120.0f, 0.0f))) {
            show = false;
            ImGui::CloseCurrentPopup();
        }

        ImGui::EndPopup();
    } else {
        show = false;
    }
}

void draw_shortcuts_dialog(bool& show) {
    if (show) {
        ImGui::OpenPopup("Keyboard Shortcuts");
    }

    bool open = show;
    if (ImGui::BeginPopupModal("Keyboard Shortcuts", &open, ImGuiWindowFlags_AlwaysAutoResize)) {
        show = open;

        draw_shortcuts_table_section("General", {
            {"Ctrl+O", "Open file"},
            {"Ctrl+S", "Save project"},
            {"Ctrl+Z", "Undo"},
            {"Ctrl+Shift+Z", "Redo"},
            {"Ctrl+A", "Select all"},
            {"Escape", "Clear selection / cancel"},
            {"Ctrl+,", "Preferences"},
        });
        draw_shortcuts_table_section("Camera", {
            {"Left drag", "Rotate"},
            {"Middle drag", "Pan"},
            {"Scroll", "Zoom"},
        });
        draw_shortcuts_table_section("Editor", {
            {"1", "Select mode"},
            {"2", "Draw mode"},
            {"3", "Erase mode"},
            {"4", "Measure mode"},
            {"5", "Fragment mode"},
            {"Delete", "Delete selected"},
            {"H", "Add hydrogens to selected"},
        });
        draw_shortcuts_table_section("Viewport", {
            {"V", "Toggle volume/isosurface"},
            {"P", "Toggle phase coloring"},
            {"L", "Toggle atom labels"},
        });

        if (ImGui::Button("OK", ImVec2(120.0f, 0.0f))) {
            show = false;
            ImGui::CloseCurrentPopup();
        }

        ImGui::EndPopup();
    } else {
        show = false;
    }
}

}  // namespace sbox::ui
