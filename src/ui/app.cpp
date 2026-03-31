#include "ui/app.h"

#include "core/logging.h"
#include "core/settings.h"
#include "core/update_checker.h"
#include "ui/about_dialog.h"
#include "ui/bond_order_panel.h"
#include "core/gaussian_eval.h"
#include "core/hydrogen.h"
#include "core/molden_parser.h"
#include "core/paths.h"
#include "io/project_io.h"
#include "ui/charge_overlay.h"
#include "ui/complex_builder.h"
#include "ui/computation_panel.h"
#include "ui/constraint_editor.h"
#include "ui/context_menu.h"
#include "ui/crystal_field_panel.h"
#include "ui/d_orbital_viewer.h"
#include "ui/dos_panel.h"
#include "ui/editor_toolbar.h"
#include "ui/esp_controls.h"
#include "ui/file_dialog.h"
#include "ui/ir_spectrum_panel.h"
#include "ui/mo_diagram.h"
#include "ui/mo_diagram_plot.h"
#include "ui/molecule_info.h"
#include "ui/optimization_panel.h"
#include "ui/orbital_composition_panel.h"
#include "ui/panels.h"
#include "ui/pes_panel.h"
#include "ui/plot_utils.h"
#include "ui/population_panel.h"
#include "ui/property_dashboard.h"
#include "ui/settings_panel.h"
#include "ui/results_panel.h"
#include "ui/setup_wizard.h"
#include "ui/solvent_panel.h"
#include "ui/spectrochemical_panel.h"

#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <implot.h>
#include <imgui_internal.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <vector>

namespace sbox {

namespace {

std::unique_ptr<Shader> try_load_shader(const std::string& vertex,
                                        const std::string& fragment,
                                        const char* label) {
    try {
        auto shader = std::make_unique<Shader>(vertex, fragment);
        int info_log_len = 0;
        glGetProgramiv(shader->id(), GL_INFO_LOG_LENGTH, &info_log_len);
        if (info_log_len > 1) {
            std::vector<char> info_log(static_cast<std::size_t>(info_log_len));
            glGetProgramInfoLog(shader->id(), info_log_len, nullptr, info_log.data());
            SBOX_LOG_DEBUG("%s shader compile/link log: %s", label, info_log.data());
        } else {
            SBOX_LOG_DEBUG("%s shader compile/link log: <empty>", label);
        }
        return shader;
    } catch (const std::exception& ex) {
        SBOX_LOG_WARN("%s shader unavailable: %s", label, ex.what());
        return nullptr;
    }
}

void clear_mo_summary(ui::AppState& state) {
    state.num_mo = 0;
    state.homo_index = -1;
    state.selected_mo = -1;
    state.mol_num_basis = 0;
    state.mol_total_energy_h = 0.0;
    state.mol_homo_lumo_gap_ev = 0.0;
    state.mol_has_mo_summary = false;
}

const std::vector<double>* current_charges_for_render(const std::optional<sbox::backend::JobResult>& result) {
    if (result.has_value() && !result->mulliken_charges.empty()) {
        return &result->mulliken_charges;
    }
    return nullptr;
}

bool is_transition_metal_Z(int Z) {
    return (Z >= 21 && Z <= 30) || (Z >= 39 && Z <= 48) || (Z >= 57 && Z <= 80) || (Z >= 89 && Z <= 112);
}

int first_transition_metal_index(const sbox::chem::MolecularSystem& mol) {
    for (int i = 0; i < mol.num_atoms(); ++i) {
        if (is_transition_metal_Z(mol.atom(i).Z)) {
            return i;
        }
    }
    return -1;
}

void open_in_system_viewer(const std::string& path) {
#ifdef __APPLE__
    const std::string command = "open \"" + path + "\"";
    std::system(command.c_str());
#elif defined(__linux__)
    const std::string command = "xdg-open \"" + path + "\"";
    std::system(command.c_str());
#elif defined(_WIN32)
    const std::string command = "start \"\" \"" + path + "\"";
    std::system(command.c_str());
#endif
}

void open_url_in_system_browser(const std::string& url) {
#ifdef __APPLE__
    const std::string command = "open \"" + url + "\"";
    std::system(command.c_str());
#elif defined(__linux__)
    const std::string command = "xdg-open \"" + url + "\"";
    std::system(command.c_str());
#elif defined(_WIN32)
    const std::string command = "start \"\" \"" + url + "\"";
    std::system(command.c_str());
#endif
}

void draw_update_dialog(bool& show, const sbox::UpdateInfo& info, sbox::SettingsManager& settings_manager) {
    if (show) {
        ImGui::OpenPopup("Update Available");
    }
    bool open = show;
    if (ImGui::BeginPopupModal("Update Available", &open, ImGuiWindowFlags_AlwaysAutoResize)) {
        show = open;
        ImGui::TextUnformatted("A new version of Schrödinger's Sandbox is available!");
        ImGui::Spacing();
        ImGui::Text("Current version: %s", info.current_version.c_str());
        ImGui::Text("Latest version: %s", info.latest_version.c_str());
        ImGui::Text("Released: %s", info.published_date.c_str());
        ImGui::Spacing();
        ImGui::TextUnformatted("What's new:");
        ImGui::TextWrapped("%s", info.changelog.empty() ? "No release notes provided." : info.changelog.c_str());
        ImGui::Spacing();
        if (ImGui::Button("Download")) {
            open_url_in_system_browser(!info.download_url.empty() ? info.download_url : info.release_url);
        }
        ImGui::SameLine();
        if (ImGui::Button("View Release")) {
            open_url_in_system_browser(info.release_url);
        }
        if (ImGui::Button("Remind Me Later")) {
            show = false;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Skip This Version")) {
            settings_manager.settings().skipped_update_version = info.latest_version;
            settings_manager.save();
            show = false;
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    } else {
        show = false;
    }
}

std::string default_screenshot_path() {
    const char* home = std::getenv("HOME");
    std::filesystem::path base = home != nullptr ? std::filesystem::path(home) : std::filesystem::current_path();
    const std::filesystem::path desktop = base / "Desktop";
    const std::filesystem::path pictures = base / "Pictures";
    if (std::filesystem::exists(desktop)) {
        base = desktop;
    } else if (std::filesystem::exists(pictures)) {
        base = pictures;
    }

    std::time_t now = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &now);
#else
    localtime_r(&now, &tm);
#endif
    char buffer[64];
    std::strftime(buffer, sizeof(buffer), "SchrodingersSandbox_%Y%m%d_%H%M%S.png", &tm);
    return (base / buffer).string();
}

void draw_color_legend(const ui::AppState& state,
                       const sbox::io::PDBData& pdb_data,
                       const std::optional<sbox::backend::JobResult>& result,
                       const ImVec2& viewport_pos) {
    if (state.view_mode != ui::ViewMode::MolecularOrbital) {
        return;
    }
    ImDrawList* draw = ImGui::GetForegroundDrawList();
    const ImVec2 box_min(viewport_pos.x + 12.0f, viewport_pos.y + 12.0f);
    const ImVec2 box_max(box_min.x + 220.0f, box_min.y + 94.0f);
    draw->AddRectFilled(box_min, box_max, IM_COL32(18, 20, 28, 210), 6.0f);
    draw->AddRect(box_min, box_max, IM_COL32(90, 96, 110, 210), 6.0f);

    auto draw_gradient = [&](float y, const Eigen::Vector3f& left, const Eigen::Vector3f& mid, const Eigen::Vector3f& right, const char* l, const char* c, const char* r) {
        const float x0 = box_min.x + 12.0f;
        const float width = 160.0f;
        for (int i = 0; i < 80; ++i) {
            const float t0 = static_cast<float>(i) / 80.0f;
            const float t1 = static_cast<float>(i + 1) / 80.0f;
            Eigen::Vector3f col = t0 < 0.5f
                                      ? left + (mid - left) * (t0 * 2.0f)
                                      : mid + (right - mid) * ((t0 - 0.5f) * 2.0f);
            draw->AddRectFilled(ImVec2(x0 + width * t0, y), ImVec2(x0 + width * t1, y + 12.0f),
                                IM_COL32(static_cast<int>(col.x() * 255.0f),
                                         static_cast<int>(col.y() * 255.0f),
                                         static_cast<int>(col.z() * 255.0f), 255));
        }
        draw->AddText(ImVec2(x0, y + 16.0f), IM_COL32(225, 228, 235, 230), l);
        draw->AddText(ImVec2(x0 + 68.0f, y + 16.0f), IM_COL32(225, 228, 235, 230), c);
        draw->AddText(ImVec2(x0 + 130.0f, y + 16.0f), IM_COL32(225, 228, 235, 230), r);
    };

    switch (static_cast<sbox::render::ColorMode>(state.color_mode)) {
    case sbox::render::ColorMode::ByChain: {
        draw->AddText(ImVec2(box_min.x + 12.0f, box_min.y + 8.0f), IM_COL32(230, 234, 242, 240), "Chains");
        const int limit = std::min<int>(static_cast<int>(pdb_data.chains.size()), 5);
        for (int i = 0; i < limit; ++i) {
            const Eigen::Vector3f col = sbox::render::chain_color(i);
            const ImU32 c = IM_COL32(static_cast<int>(col.x() * 255.0f),
                                     static_cast<int>(col.y() * 255.0f),
                                     static_cast<int>(col.z() * 255.0f), 255);
            const float y = box_min.y + 30.0f + i * 12.0f;
            draw->AddRectFilled(ImVec2(box_min.x + 12.0f, y), ImVec2(box_min.x + 24.0f, y + 8.0f), c, 2.0f);
            const std::string label = "Chain " + (pdb_data.chains[i].id.empty() ? std::string("?") : pdb_data.chains[i].id);
            draw->AddText(ImVec2(box_min.x + 30.0f, y - 3.0f), IM_COL32(225, 228, 235, 230), label.c_str());
        }
        break;
    }
    case sbox::render::ColorMode::ByResidue: {
        draw->AddText(ImVec2(box_min.x + 12.0f, box_min.y + 8.0f), IM_COL32(230, 234, 242, 240), "Residues");
        const std::array<const char*, 5> residues = {"ALA", "ASP", "LYS", "SER", "HOH"};
        for (int i = 0; i < static_cast<int>(residues.size()); ++i) {
            const Eigen::Vector3f col = sbox::render::residue_color(residues[static_cast<std::size_t>(i)]);
            const ImU32 c = IM_COL32(static_cast<int>(col.x() * 255.0f),
                                     static_cast<int>(col.y() * 255.0f),
                                     static_cast<int>(col.z() * 255.0f), 255);
            const float y = box_min.y + 30.0f + i * 12.0f;
            draw->AddRectFilled(ImVec2(box_min.x + 12.0f, y), ImVec2(box_min.x + 24.0f, y + 8.0f), c, 2.0f);
            draw->AddText(ImVec2(box_min.x + 30.0f, y - 3.0f), IM_COL32(225, 228, 235, 230), residues[static_cast<std::size_t>(i)]);
        }
        break;
    }
    case sbox::render::ColorMode::ByBFactor:
        draw->AddText(ImVec2(box_min.x + 12.0f, box_min.y + 8.0f), IM_COL32(230, 234, 242, 240), "B-Factor");
        draw_gradient(box_min.y + 32.0f,
                      Eigen::Vector3f(0.0f, 0.0f, 0.8f),
                      Eigen::Vector3f(1.0f, 1.0f, 1.0f),
                      Eigen::Vector3f(0.8f, 0.0f, 0.0f),
                      "low", "mid", "high");
        break;
    case sbox::render::ColorMode::ByCharge:
        draw->AddText(ImVec2(box_min.x + 12.0f, box_min.y + 8.0f), IM_COL32(230, 234, 242, 240), "Charge");
        draw_gradient(box_min.y + 32.0f,
                      Eigen::Vector3f(0.15f, 0.35f, 0.95f),
                      Eigen::Vector3f(1.0f, 1.0f, 1.0f),
                      Eigen::Vector3f(0.95f, 0.20f, 0.15f),
                      "-", "0", "+");
        if (!result.has_value() || result->mulliken_charges.empty()) {
            draw->AddText(ImVec2(box_min.x + 12.0f, box_min.y + 70.0f), IM_COL32(225, 160, 120, 230), "No charges loaded");
        }
        break;
    default:
        break;
    }
}

std::optional<sbox::basis::MOData> mo_data_from_fchk(const sbox::io::FchkData& fchk) {
    if (fchk.shell_types.empty() || fchk.shell_to_atom_map.empty() || fchk.primitives_per_shell.empty() ||
        fchk.primitive_exponents.empty() || fchk.contraction_coefficients.empty() || fchk.mo_coefficients.size() == 0) {
        return std::nullopt;
    }

    if (fchk.shell_types.size() != fchk.shell_to_atom_map.size() ||
        fchk.shell_types.size() != fchk.primitives_per_shell.size()) {
        return std::nullopt;
    }

    sbox::basis::MOData mo_data;
    mo_data.atomic_numbers = fchk.atomic_numbers;
    mo_data.energies = fchk.mo_energies;
    mo_data.occupations = fchk.occupations;
    mo_data.coefficients = fchk.mo_coefficients;
    mo_data.total_energy = fchk.total_energy;
    mo_data.atom_positions.reserve(static_cast<std::size_t>(fchk.num_atoms));
    for (int i = 0; i < fchk.num_atoms; ++i) {
        mo_data.atom_positions.emplace_back(fchk.coordinates[static_cast<std::size_t>(3 * i + 0)],
                                            fchk.coordinates[static_cast<std::size_t>(3 * i + 1)],
                                            fchk.coordinates[static_cast<std::size_t>(3 * i + 2)]);
    }

    bool spherical = true;
    for (int shell_type : fchk.shell_types) {
        if (shell_type > 1) {
            spherical = false;
        }
        if (shell_type == -1) {
            continue;
        }
    }
    mo_data.basis.spherical = spherical;

    std::size_t primitive_offset = 0;
    for (std::size_t shell_index = 0; shell_index < fchk.shell_types.size(); ++shell_index) {
        const int shell_type = fchk.shell_types[shell_index];
        const int atom_index = fchk.shell_to_atom_map[shell_index] - 1;
        const int nprim = fchk.primitives_per_shell[shell_index];
        if (atom_index < 0 || atom_index >= fchk.num_atoms || nprim < 0) {
            return std::nullopt;
        }
        if (primitive_offset + static_cast<std::size_t>(nprim) > fchk.primitive_exponents.size() ||
            primitive_offset + static_cast<std::size_t>(nprim) > fchk.contraction_coefficients.size()) {
            return std::nullopt;
        }

        auto fill_shell = [&](sbox::basis::BasisShell& shell, const std::vector<double>& coeffs) {
            shell.atom_index = atom_index;
            shell.primitives.reserve(static_cast<std::size_t>(nprim));
            for (int p = 0; p < nprim; ++p) {
                shell.primitives.push_back({fchk.primitive_exponents[primitive_offset + static_cast<std::size_t>(p)],
                                            coeffs[primitive_offset + static_cast<std::size_t>(p)]});
            }
        };

        if (shell_type == -1) {
            if (primitive_offset + static_cast<std::size_t>(nprim) > fchk.sp_contraction_coefficients.size()) {
                return std::nullopt;
            }
            sbox::basis::BasisShell shell_s;
            shell_s.angular_momentum = 0;
            fill_shell(shell_s, fchk.contraction_coefficients);
            mo_data.basis.shells.push_back(shell_s);

            sbox::basis::BasisShell shell_p;
            shell_p.angular_momentum = 1;
            fill_shell(shell_p, fchk.sp_contraction_coefficients);
            mo_data.basis.shells.push_back(shell_p);
        } else {
            sbox::basis::BasisShell shell;
            shell.angular_momentum = std::abs(shell_type);
            if (shell.angular_momentum < 0 || shell.angular_momentum > 3) {
                return std::nullopt;
            }
            fill_shell(shell, fchk.contraction_coefficients);
            mo_data.basis.shells.push_back(shell);
        }

        primitive_offset += static_cast<std::size_t>(nprim);
    }

    if (mo_data.basis.num_basis_functions() != mo_data.coefficients.rows()) {
        return std::nullopt;
    }

    return mo_data;
}

}  // namespace

App::App() {
    settings_manager_.load();
    const auto& settings = settings_manager_.settings();

    window_ = std::make_unique<Window>(settings.window_width, settings.window_height, "Schrödinger's Sandbox");
    if (settings.window_maximized) {
        glfwMaximizeWindow(window_->handle());
    }
    glfwSwapInterval(settings.enable_vsync ? 1 : 0);
    gradient_shader_ = std::make_unique<Shader>(sbox::get_shader_path("fullscreen_quad.vert"),
                                                sbox::get_shader_path("test_gradient.frag"));
    orbital_shader_ = try_load_shader(sbox::get_shader_path("orbital_raymarch.vert"),
                                      sbox::get_shader_path("orbital_raymarch.frag"),
                                      "Orbital");
    mo_shader_ = try_load_shader(sbox::get_shader_path("mo_raymarch.vert"),
                                 sbox::get_shader_path("mo_raymarch.frag"),
                                 "MO");
    cube_shader_ = try_load_shader(sbox::get_shader_path("cube_raymarch.vert"),
                                   sbox::get_shader_path("cube_raymarch.frag"),
                                   "Cube");
    esp_shader_ = try_load_shader(sbox::get_shader_path("mo_raymarch.vert"),
                                  sbox::get_shader_path("esp_surface.frag"),
                                  "ESP Surface");

    glGenVertexArrays(1, &fullscreen_vao_);

    python_env_.detect();
    python_env_.check_packages();
    backend_.init(python_env_);
    if (!python_env_.has_pyscf() && !python_env_.has_tblite()) {
        wizard_state_.show_wizard = true;
    }

    glfwSetWindowUserPointer(window_->handle(), this);
    glfwSetScrollCallback(window_->handle(), &App::ScrollCallback);

    editor_state_.select_mode = std::make_unique<sbox::editor::SelectMode>();
    editor_state_.draw_mode = std::make_unique<sbox::editor::DrawMode>();
    editor_state_.erase_mode = std::make_unique<sbox::editor::EraseMode>();
    editor_state_.measure_mode = std::make_unique<sbox::editor::MeasureMode>();
    editor_state_.fragment_mode = std::make_unique<sbox::editor::FragmentMode>(&editor_state_.fragment_library);
    editor_state_.select_mode->set_context_menu_state(&editor_state_.context_menu);

    state_.iso_value = settings.default_iso_value;
    state_.gamma = settings.default_gamma;
    state_.mol_render_mode = settings.mol_render_mode;
    state_.color_mode = settings.color_mode;
    state_.computation.method = static_cast<sbox::backend::Method>(std::clamp(settings.default_method, 0, static_cast<int>(sbox::backend::Method::GFN_FF)));
    state_.computation.basis = static_cast<sbox::backend::BasisSetType>(std::clamp(settings.default_basis, 0, static_cast<int>(sbox::backend::BasisSetType::aug_cc_pVTZ)));
    state_.computation.charge = settings.default_charge;
    state_.computation.multiplicity = settings.default_multiplicity;
    sbox::render::set_atom_radius_scale(settings.atom_scale);
    sbox::render::set_bond_radius_scale(settings.bond_scale);
    if (!settings.python_path.empty()) {
        python_env_.set_python_path(settings.python_path);
        python_env_.detect();
        python_env_.check_packages();
        backend_.init(python_env_);
    }

    ui::set_about_dialog_context(&python_env_);
    if (settings.check_for_updates) {
        update_checker_ = std::make_unique<sbox::UpdateChecker>("Manav02012002/SchrodingersSandbox");
        update_checker_->check_async();
    }

    initImGui();
    ensureViewportTarget(viewport_width_, viewport_height_);
}

App::~App() {
    if (window_ != nullptr) {
        int win_w = settings_manager_.settings().window_width;
        int win_h = settings_manager_.settings().window_height;
        glfwGetWindowSize(window_->handle(), &win_w, &win_h);
        settings_manager_.settings().window_width = win_w;
        settings_manager_.settings().window_height = win_h;
        settings_manager_.settings().window_maximized =
            glfwGetWindowAttrib(window_->handle(), GLFW_MAXIMIZED) == GLFW_TRUE;
        settings_manager_.settings().python_path = python_env_.info().python_path;
        settings_manager_.settings().default_method = static_cast<int>(state_.computation.method);
        settings_manager_.settings().default_basis = static_cast<int>(state_.computation.basis);
        settings_manager_.settings().default_charge = state_.computation.charge;
        settings_manager_.settings().default_multiplicity = state_.computation.multiplicity;
        settings_manager_.settings().default_iso_value = state_.iso_value;
        settings_manager_.settings().default_gamma = state_.gamma;
        settings_manager_.settings().mol_render_mode = state_.mol_render_mode;
        settings_manager_.settings().color_mode = state_.color_mode;
    }
    settings_manager_.save();

    shutdownImGui();

    if (fullscreen_vao_ != 0U) {
        glDeleteVertexArrays(1, &fullscreen_vao_);
        fullscreen_vao_ = 0;
    }

    if (viewport_depth_rbo_ != 0U) {
        glDeleteRenderbuffers(1, &viewport_depth_rbo_);
        viewport_depth_rbo_ = 0;
    }

    if (viewport_color_tex_ != 0U) {
        glDeleteTextures(1, &viewport_color_tex_);
        viewport_color_tex_ = 0;
    }

    if (viewport_fbo_ != 0U) {
        glDeleteFramebuffers(1, &viewport_fbo_);
        viewport_fbo_ = 0;
    }
}

void App::run() {
    static bool dock_layout_built = false;

    while (!window_->shouldClose()) {
        window_->pollEvents();
        ++update_poll_frame_counter_;
        if (update_checker_ && update_poll_frame_counter_ % 60 == 0) {
            const std::optional<sbox::UpdateInfo> result = update_checker_->get_result();
            if (result.has_value() &&
                result->update_available &&
                result->latest_version != settings_manager_.settings().skipped_update_version &&
                !update_notification_shown_) {
                update_notification_shown_ = true;
                pending_update_ = *result;
            }
        }

        const std::vector<int> completed_jobs = backend_.poll_completed();
        for (int job_id : completed_jobs) {
            const sbox::backend::JobResult* job_result = backend_.result(job_id);
            if (job_id == state_.solvent.gas_job_id && job_result != nullptr) {
                state_.solvent.gas_result = *job_result;
                state_.solvent.gas_done = true;
                continue;
            }
            if (job_id == state_.pes.scan_job_id && job_result != nullptr) {
                state_.pes.result = job_result->scan_result;
                state_.pes.scan_running = false;
                state_.pes.scan_complete = job_result->has_scan;
                if (!state_.pes.result.geometries.empty()) {
                    state_.pes.viewed_point = std::clamp(
                        state_.pes.viewed_point,
                        0,
                        static_cast<int>(state_.pes.result.geometries.size()) - 1);
                }
                continue;
            }
            if (job_id == state_.solvent.solvent_job_id && job_result != nullptr) {
                state_.solvent.solvent_result = *job_result;
                state_.solvent.solvent_done = true;
                if (state_.solvent.gas_done && state_.solvent.solvent_result.converged() && state_.solvent.gas_result.converged()) {
                    const auto& solvent = state_.solvent.solvent_result;
                    const auto& gas = state_.solvent.gas_result;
                    const auto& opt = state_.solvent.selected_solvent_index;
                    const std::array<std::pair<const char*, double>, 10> solvents = {{
                        {"Gas Phase", 1.0}, {"Water", 78.4}, {"DMSO", 46.7}, {"Acetonitrile", 35.7},
                        {"Methanol", 32.7}, {"Ethanol", 24.3}, {"DCM", 8.9}, {"THF", 7.6},
                        {"Toluene", 2.4}, {"Hexane", 1.9}
                    }};
                    const auto& sel = solvents[static_cast<std::size_t>(std::clamp(opt, 0, static_cast<int>(solvents.size()) - 1))];
                    state_.solvent.history.push_back({
                        sel.first,
                        sel.second,
                        solvent.total_energy,
                        (solvent.total_energy - gas.total_energy) * 627.509,
                        solvent.dipole_moment.norm(),
                        solvent.homo_index() >= 0 ? solvent.mo_data.energies(solvent.homo_index()) * 27.2114 : 0.0,
                        solvent.lumo_index() >= 0 ? solvent.mo_data.energies(solvent.lumo_index()) * 27.2114 : 0.0
                    });
                }
                continue;
            }
            state_.computation.active_job_id = job_id;
            state_.computation.job_running = false;
            state_.computation.job_completed = true;
            state_.computation.last_progress_iteration = 0;
            if (job_result != nullptr) {
                latest_result_ = *job_result;
                state_.computation.last_error = job_result->error_message;
                if (job_result->status == sbox::backend::JobStatus::Converged) {
                    try {
                        applyBackendResult(*job_result);
                        if (job_result->has_mo_data && first_transition_metal_index(current_molecule_) >= 0) {
                            run_crystal_field_analysis();
                            if (has_d_orbital_analysis_) {
                                message_popup_title_ = "Crystal Field Analysis";
                                message_popup_text_ = "Crystal field analysis available — see Crystal Field Diagram and d-Orbital Viewer.";
                                open_message_popup_ = true;
                            }
                        }
                    } catch (const std::exception& ex) {
                        state_.computation.last_error = ex.what();
                    }
                }
            }
        }

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        const bool wizard_active = ui::draw_setup_wizard(wizard_state_, python_env_);
        backend_.init(python_env_);
        (void)wizard_active;

        if (ImGui::BeginMainMenuBar()) {
            if (ImGui::BeginMenu("File")) {
                if (ImGui::MenuItem("New Molecule")) {
                    current_molecule_.clear();
                    current_pdb_data_ = sbox::io::PDBData{};
                    current_molecule_.set_name("Untitled");
                    editor_state_.commands.clear();
                    editor_state_.selection.clear();
                    uploadCurrentMoleculeToRenderers();
                    state_.view_mode = ui::ViewMode::MolecularOrbital;
                    has_mo_data_ = false;
                    has_cube_data_ = false;
                    use_cube_fallback_ = false;
                    clear_mo_summary(state_);
                    state_.molecule_loaded = false;
                    state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);
                }
                if (ImGui::MenuItem("Open Molden...")) {
                    try {
                        const std::string path = ui::open_file_dialog("Open Molden File", "molden");
                        if (!path.empty()) {
                            loadMoldenFile(path);
                        }
                    } catch (const std::exception& ex) {
                        SBOX_LOG_ERROR("Failed to load Molden file: %s", ex.what());
                    }
                }
                if (ImGui::MenuItem("Open XYZ...")) {
                    try {
                        const std::string path = ui::open_file_dialog("Open XYZ File", "xyz");
                        if (!path.empty()) {
                            loadXYZFile(path);
                        }
                    } catch (const std::exception& ex) {
                        SBOX_LOG_ERROR("Failed to load XYZ file: %s", ex.what());
                    }
                }
                if (ImGui::MenuItem("Open Trajectory...")) {
                    try {
                        const std::string path = ui::open_file_dialog("Open Trajectory File", "xyz");
                        if (!path.empty()) {
                            loadTrajectoryFile(path);
                        }
                    } catch (const std::exception& ex) {
                        SBOX_LOG_ERROR("Failed to load trajectory file: %s", ex.what());
                    }
                }
                if (ImGui::MenuItem("Open SDF...")) {
                    try {
                        const std::string path = ui::open_file_dialog("Open SDF File", "sdf,mol");
                        if (!path.empty()) {
                            loadSDFFile(path);
                        }
                    } catch (const std::exception& ex) {
                        SBOX_LOG_ERROR("Failed to load SDF file: %s", ex.what());
                    }
                }
                if (ImGui::MenuItem("Open PDB...")) {
                    try {
                        const std::string path = ui::open_file_dialog("Open PDB File", "pdb,ent");
                        if (!path.empty()) {
                            loadPDBFile(path);
                        }
                    } catch (const std::exception& ex) {
                        SBOX_LOG_ERROR("Failed to load PDB file: %s", ex.what());
                    }
                }
                if (ImGui::MenuItem("Open Cube...")) {
                    try {
                        const std::string path = ui::open_file_dialog("Open Cube File", "cube");
                        if (!path.empty()) {
                            loadCubeFile(path);
                        }
                    } catch (const std::exception& ex) {
                        SBOX_LOG_ERROR("Failed to load Cube file: %s", ex.what());
                    }
                }
                if (ImGui::MenuItem("Open FCHK...")) {
                    try {
                        const std::string path = ui::open_file_dialog("Open FCHK File", "fchk,fch");
                        if (!path.empty()) {
                            loadFchkFile(path);
                        }
                    } catch (const std::exception& ex) {
                        SBOX_LOG_ERROR("Failed to load FCHK file: %s", ex.what());
                    }
                }
                if (ImGui::MenuItem("Save Screenshot...", "Ctrl+Shift+S")) {
                    const std::string path = ui::save_file_dialog("Save Screenshot", "png,jpg,bmp", "SchrodingersSandbox.png");
                    if (!path.empty()) {
                        renderViewportToTarget(viewport_fbo_, viewport_width_, viewport_height_, false);
                        if (sbox::render::save_screenshot(path, viewport_fbo_, viewport_width_, viewport_height_)) {
                            SBOX_LOG_INFO("Screenshot saved to %s", path.c_str());
                        } else {
                            SBOX_LOG_ERROR("Failed to save screenshot to %s", path.c_str());
                        }
                    }
                }
                if (ImGui::BeginMenu("Export Image")) {
                    if (ImGui::MenuItem("1080p (1920x1080)")) {
                        const std::string path = ui::save_file_dialog("Save Screenshot", "png", "SchrodingersSandbox_1080p.png");
                        if (!path.empty()) {
                            if (sbox::render::save_screenshot_highres(path, 1920, 1080, [this](unsigned int fbo, int w, int h) {
                                    renderViewportToTarget(fbo, w, h, false);
                                })) {
                                SBOX_LOG_INFO("High-resolution screenshot saved to %s", path.c_str());
                            } else {
                                SBOX_LOG_ERROR("Failed to save high-resolution screenshot to %s", path.c_str());
                            }
                        }
                    }
                    if (ImGui::MenuItem("4K (3840x2160)")) {
                        const std::string path = ui::save_file_dialog("Save Screenshot", "png", "SchrodingersSandbox_4k.png");
                        if (!path.empty()) {
                            if (sbox::render::save_screenshot_highres(path, 3840, 2160, [this](unsigned int fbo, int w, int h) {
                                    renderViewportToTarget(fbo, w, h, false);
                                })) {
                                SBOX_LOG_INFO("4K screenshot saved to %s", path.c_str());
                            } else {
                                SBOX_LOG_ERROR("Failed to save 4K screenshot to %s", path.c_str());
                            }
                            ensureViewportTarget(viewport_width_, viewport_height_);
                        }
                    }
                    if (ImGui::MenuItem("Transparent Background (PNG)")) {
                        const std::string path = ui::save_file_dialog("Save Screenshot", "png", "SchrodingersSandbox_transparent.png");
                        if (!path.empty()) {
                            renderViewportToTarget(viewport_fbo_, viewport_width_, viewport_height_, true);
                            if (sbox::render::save_screenshot_transparent(path, viewport_fbo_, viewport_width_, viewport_height_)) {
                                SBOX_LOG_INFO("Transparent screenshot saved to %s", path.c_str());
                            } else {
                                SBOX_LOG_ERROR("Failed to save transparent screenshot to %s", path.c_str());
                            }
                            renderViewportToTarget(viewport_fbo_, viewport_width_, viewport_height_, false);
                        }
                    }
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Recent Files")) {
                    for (const auto& path : settings_manager_.recent_files()) {
                        const std::string filename = std::filesystem::path(path).filename().string();
                        if (ImGui::MenuItem(filename.c_str())) {
                            try {
                                load_file_by_extension(path);
                            } catch (const std::exception& ex) {
                                SBOX_LOG_ERROR("Failed to load recent file: %s", ex.what());
                            }
                        }
                    }
                    if (!settings_manager_.recent_files().empty()) {
                        ImGui::Separator();
                    }
                    if (ImGui::MenuItem("Clear Recent Files")) {
                        settings_manager_.clear_recent_files();
                    }
                    ImGui::EndMenu();
                }
                ImGui::Separator();
                if (ImGui::MenuItem("Save Project...")) {
                    (void)ui::open_file_dialog("Save Project", "sbox");
                }
                ImGui::Separator();
                if (ImGui::MenuItem("Quit")) {
                    glfwSetWindowShouldClose(window_->handle(), GLFW_TRUE);
                }
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Edit")) {
                if (ImGui::MenuItem("Undo", "Ctrl+Z", false, editor_state_.commands.can_undo())) {
                    editor_state_.commands.undo(current_molecule_);
                }
                if (ImGui::MenuItem("Redo", "Ctrl+Shift+Z", false, editor_state_.commands.can_redo())) {
                    editor_state_.commands.redo(current_molecule_);
                }
                ImGui::Separator();
                if (ImGui::MenuItem("Select All", "Ctrl+A")) {
                    editor_state_.selection.clear();
                    for (int i = 0; i < current_molecule_.num_atoms(); ++i) {
                        editor_state_.selection.toggle_atom(i);
                    }
                }
                if (ImGui::MenuItem("Deselect All", "Escape")) {
                    editor_state_.selection.clear();
                }
                ImGui::Separator();
                if (ImGui::MenuItem("Add Hydrogens")) {
                    editor_state_.commands.execute(
                        std::make_unique<sbox::editor::AddHydrogensCommand>(-1), current_molecule_);
                }
                if (ImGui::MenuItem("Remove Hydrogens")) {
                    editor_state_.commands.execute(
                        std::make_unique<sbox::editor::RemoveHydrogensCommand>(), current_molecule_);
                }
                ImGui::Separator();
                if (ImGui::MenuItem("Preferences...", "Cmd+,")) {
                    state_.show_settings = true;
                }
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("View")) {
                if (ImGui::MenuItem("Reset Layout")) {
                    dock_layout_built = false;
                    const ImGuiID reset_dockspace_id = ImGui::GetID("MainDockSpace");
                    ImGui::DockBuilderRemoveNode(reset_dockspace_id);
                }
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Chemistry")) {
                if (ImGui::MenuItem("Complex Builder...")) {
                    state_.show_complex_builder = true;
                }
                if (ImGui::MenuItem("Crystal Field Analysis")) {
                    run_crystal_field_analysis();
                }
                if (ImGui::MenuItem("Spectrochemical Series...")) {
                    state_.show_spectrochemical = true;
                }
                ImGui::Separator();
                if (ImGui::MenuItem("Detect Metal Center")) {
                    detect_metal_center();
                }
                ImGui::EndMenu();
            }

            if (ImGui::BeginMenu("Help")) {
                if (ImGui::MenuItem("Keyboard Shortcuts...")) {
                    show_shortcuts_ = true;
                }
                ImGui::Separator();
                if (ImGui::MenuItem("View Log File")) {
                    open_in_system_viewer((std::filesystem::path(get_app_data_dir()) / "schrodingers_sandbox.log").string());
                }
                if (ImGui::MenuItem("Open Data Directory")) {
                    open_in_system_viewer(get_app_data_dir());
                }
                ImGui::Separator();
                if (ImGui::MenuItem("About...")) {
                    show_about_ = true;
                }
                ImGui::EndMenu();
            }
            ImGui::EndMainMenuBar();
        }

        if ((ImGui::GetIO().KeySuper || ImGui::GetIO().KeyCtrl) && ImGui::IsKeyPressed(ImGuiKey_Comma)) {
            state_.show_settings = true;
        }
        if ((ImGui::GetIO().KeySuper || ImGui::GetIO().KeyCtrl) &&
            ImGui::GetIO().KeyShift &&
            ImGui::IsKeyPressed(ImGuiKey_S)) {
            const std::string path = default_screenshot_path();
            renderViewportToTarget(viewport_fbo_, viewport_width_, viewport_height_, false);
            if (sbox::render::save_screenshot(path, viewport_fbo_, viewport_width_, viewport_height_)) {
                SBOX_LOG_INFO("Screenshot saved to %s", path.c_str());
            } else {
                SBOX_LOG_ERROR("Failed to save screenshot to %s", path.c_str());
            }
        }

        if (open_message_popup_) {
            ImGui::OpenPopup("##chem_message_popup");
            open_message_popup_ = false;
        }
        if (ImGui::BeginPopupModal("##chem_message_popup", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::Text("%s", message_popup_title_.c_str());
            ImGui::Separator();
            ImGui::TextWrapped("%s", message_popup_text_.c_str());
            if (ImGui::Button("OK", ImVec2(120.0f, 0.0f))) {
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }

        if (open_metal_picker_popup_) {
            ImGui::OpenPopup("##metal_picker_popup");
            open_metal_picker_popup_ = false;
        }
        if (ImGui::BeginPopupModal("##metal_picker_popup", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::TextUnformatted("Multiple transition metals detected.");
            ImGui::TextUnformatted("Choose the metal center to analyze:");
            ImGui::Separator();
            for (int atom_index : detected_metal_choices_) {
                const auto& atom = current_molecule_.atom(atom_index);
                char label[64];
                std::snprintf(label, sizeof(label), "%s%d", sbox::elements::get_element(atom.Z).symbol, atom_index + 1);
                if (ImGui::Selectable(label)) {
                    current_metal_index_ = atom_index;
                    ImGui::CloseCurrentPopup();
                }
            }
            if (ImGui::Button("Use First")) {
                if (!detected_metal_choices_.empty()) {
                    current_metal_index_ = detected_metal_choices_.front();
                }
                ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            if (ImGui::Button("Cancel")) {
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }

        if (state_.show_settings) {
            bool changed = false;
            ImGui::SetNextWindowSize(ImVec2(600.0f, 500.0f), ImGuiCond_FirstUseEver);
            if (ImGui::Begin("Preferences", &state_.show_settings)) {
                ui::draw_settings_panel(settings_manager_.settings(), python_env_, changed);
                if (changed) {
                    apply_settings();
                }
            }
            ImGui::End();
            if (ui::consume_settings_panel_close_request()) {
                state_.show_settings = false;
            }
        }

        ui::draw_about_dialog(show_about_);
        ui::draw_shortcuts_dialog(show_shortcuts_);
        if (pending_update_.has_value() && show_update_dialog_) {
            draw_update_dialog(show_update_dialog_, *pending_update_, settings_manager_);
        }

        const ImGuiID dockspace_id = ImGui::GetID("MainDockSpace");
        ImGui::DockSpaceOverViewport(dockspace_id, ImGui::GetMainViewport());

        if (!dock_layout_built) {
            ImGuiDockNode* dockspace_node = ImGui::DockBuilderGetNode(dockspace_id);
            const bool has_existing_layout =
                dockspace_node != nullptr &&
                (dockspace_node->ChildNodes[0] != nullptr ||
                 dockspace_node->ChildNodes[1] != nullptr ||
                 dockspace_node->Windows.Size > 0);

            if (!has_existing_layout) {
                ImGui::DockBuilderRemoveNode(dockspace_id);
                ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
                ImGui::DockBuilderSetNodeSize(dockspace_id, ImGui::GetMainViewport()->WorkSize);

                ImGuiID left_node = dockspace_id;
                ImGuiID center_right_node = 0;
                ImGui::DockBuilderSplitNode(left_node, ImGuiDir_Left, 0.20f, &left_node, &center_right_node);

                ImGuiID top_center_right_node = center_right_node;
                ImGuiID bottom_center_right_node = 0;
                ImGui::DockBuilderSplitNode(top_center_right_node,
                                            ImGuiDir_Down,
                                            0.30f,
                                            &bottom_center_right_node,
                                            &top_center_right_node);

                ImGuiID center_node = top_center_right_node;
                ImGuiID right_node = 0;
                ImGui::DockBuilderSplitNode(center_node, ImGuiDir_Left, 0.6875f, &center_node, &right_node);

                ImGuiID orbital_node = left_node;
                ImGuiID properties_node = 0;
                ImGui::DockBuilderSplitNode(orbital_node, ImGuiDir_Down, 0.45f, &properties_node, &orbital_node);
                ImGuiID editor_node = left_node;
                ImGuiID browser_stack_node = 0;
                ImGui::DockBuilderSplitNode(editor_node, ImGuiDir_Down, 0.32f, &browser_stack_node, &editor_node);

                ImGui::DockBuilderDockWindow("Editor", editor_node);
                ImGui::DockBuilderDockWindow("Orbital Browser", browser_stack_node);
                ImGui::DockBuilderDockWindow("Atomic Properties", properties_node);
                ImGui::DockBuilderDockWindow("Energy Levels", properties_node);
                ImGui::DockBuilderDockWindow("MO Energy Diagram", properties_node);
                ImGui::DockBuilderDockWindow("Crystal Field Diagram", properties_node);
                ImGui::DockBuilderDockWindow("d-Orbital Viewer", properties_node);
                ImGui::DockBuilderDockWindow("Density of States", properties_node);
                ImGui::DockBuilderDockWindow("Orbital Composition", properties_node);
                ImGui::DockBuilderDockWindow("3D Viewport", center_node);
                ImGui::DockBuilderDockWindow("Computation", right_node);
                ImGui::DockBuilderDockWindow("Results", right_node);
                ImGui::DockBuilderDockWindow("Properties", right_node);
                ImGui::DockBuilderDockWindow("Population Analysis", right_node);
                ImGui::DockBuilderDockWindow("Bond Orders", right_node);
                ImGui::DockBuilderDockWindow("Electrostatic Properties", right_node);
                ImGui::DockBuilderDockWindow("IR Spectrum", right_node);
                ImGui::DockBuilderDockWindow("Solvent Effects", right_node);
                ImGui::DockBuilderDockWindow("Potential Energy Surface", right_node);
                ImGui::DockBuilderDockWindow("Complex Builder", right_node);
                ImGui::DockBuilderDockWindow("Periodic Table", bottom_center_right_node);
                ImGui::DockBuilderFinish(dockspace_id);
            }

            dock_layout_built = true;
        }

        if (state_.view_mode == ui::ViewMode::AtomicOrbital && state_.needs_update) {
            state_.update();
        }

        ui::draw_periodic_table(state_);
        if (state_.view_mode == ui::ViewMode::MolecularOrbital) {
            ui::draw_editor_toolbar(editor_state_, current_molecule_);
            ui::draw_constraint_editor(state_, current_molecule_, editor_state_.selection);
        }
        if (state_.show_complex_builder) {
            ui::draw_complex_builder(state_, current_molecule_, editor_state_.commands, ligand_library_);
        }
        ui::draw_orbital_browser(state_);
        if (state_.view_mode == ui::ViewMode::MolecularOrbital) {
            ui::draw_molecule_info(state_, current_molecule_);
        } else {
            ui::draw_properties(state_);
        }
        if (state_.view_mode == ui::ViewMode::MolecularOrbital && has_mo_data_ &&
            (!latest_result_.has_value() || !latest_result_->has_mo_data)) {
            ui::draw_mo_diagram(state_, current_mo_data_);
        } else {
            if (state_.view_mode != ui::ViewMode::MolecularOrbital) {
                ui::draw_energy_diagram(state_);
            }
        }
        ui::draw_computation_panel(state_, backend_);
        ui::draw_pes_panel(state_, backend_, current_molecule_, editor_state_.selection);
        ui::draw_solvent_panel(state_, backend_);
        if (state_.show_spectrochemical) {
            ui::draw_spectrochemical_panel(state_, backend_, ligand_library_);
        }
        if (latest_result_.has_value() && (!latest_result_->opt_history.empty() || state_.computation.job_running)) {
            ui::draw_optimization_panel(state_, *latest_result_, current_molecule_, mol_renderer_);
        }
        if (latest_result_.has_value() && latest_result_->converged()) {
            ui::draw_results_panel(state_, *latest_result_, current_molecule_);
            if (has_d_orbital_analysis_) {
                ui::draw_crystal_field_panel(state_, *latest_result_, current_molecule_);
                ui::draw_d_orbital_viewer(state_, *latest_result_, current_molecule_, current_d_orbitals_);
            }
            switch (state_.property_view) {
            case ui::PropertyView::Dashboard:
                ui::draw_property_dashboard(state_, *latest_result_, current_molecule_);
                break;
            case ui::PropertyView::Population:
                ui::draw_population_panel(state_, *latest_result_, current_molecule_);
                break;
            case ui::PropertyView::BondOrders:
                ui::draw_bond_order_panel(state_, *latest_result_, current_molecule_);
                break;
            case ui::PropertyView::MoDiagram:
                ui::draw_mo_diagram_plot(state_, *latest_result_);
                break;
            case ui::PropertyView::DOS:
                ui::draw_dos_panel(state_, *latest_result_);
                break;
            case ui::PropertyView::OrbitalComposition:
                ui::draw_orbital_composition_panel(state_, *latest_result_, current_molecule_);
                break;
            case ui::PropertyView::IRSpectrum:
                ui::draw_ir_spectrum_panel(state_, *latest_result_, current_molecule_);
                break;
            case ui::PropertyView::ESPControls:
                ui::draw_esp_controls(state_, *latest_result_);
                break;
            }
        }

        if (state_.computation.run_requested && !state_.computation.job_running && state_.molecule_loaded) {
            state_.computation.run_requested = false;
            try {
                const sbox::backend::JobSpec spec = makeJobSpecFromState();
                state_.computation.active_job_id = backend_.submit(spec);
                state_.computation.job_running = true;
                state_.computation.job_completed = false;
                state_.computation.last_error.clear();
                state_.computation.last_progress_iteration = 0;
                state_.computation.scf_plot_energies.clear();
            } catch (const std::exception& ex) {
                state_.computation.last_error = ex.what();
                state_.computation.job_running = false;
            }
        }

        if (state_.solvent.run_requested && state_.molecule_loaded) {
            state_.solvent.run_requested = false;
            const std::array<const char*, 10> solvent_values = {"", "water", "dmso", "acetonitrile", "methanol", "ethanol", "dcm", "thf", "toluene", "hexane"};
            try {
                sbox::backend::JobSpec gas_spec = makeJobSpecFromState();
                gas_spec.solvent.clear();
                gas_spec.properties = {
                    sbox::backend::PropertyRequest::MullikenCharges,
                    sbox::backend::PropertyRequest::DipoleMoment,
                    sbox::backend::PropertyRequest::MoldenFile,
                };
                sbox::backend::JobSpec solv_spec = gas_spec;
                solv_spec.solvent = solvent_values[static_cast<std::size_t>(std::clamp(state_.solvent.selected_solvent_index, 0, 9))];

                state_.solvent.gas_job_id = backend_.submit(gas_spec);
                state_.solvent.solvent_job_id = backend_.submit(solv_spec);
                state_.solvent.gas_done = false;
                state_.solvent.solvent_done = false;
            } catch (const std::exception& ex) {
                state_.computation.last_error = ex.what();
            }
        }

        if (state_.computation.apply_results_requested) {
            state_.computation.apply_results_requested = false;
            if (const sbox::backend::JobResult* job_result = backend_.result(state_.computation.active_job_id)) {
                applyBackendResult(*job_result);
            }
        }

        const ui::ViewportPanelState viewport_state = ui::draw_viewport(state_, viewport_color_tex_);

        if (state_.view_mode == ui::ViewMode::MolecularOrbital && viewport_state.hovered) {
            const ImVec2 mouse_pos = ImGui::GetIO().MousePos;
            const sbox::editor::Ray ray = sbox::editor::screen_to_ray(
                mouse_pos.x,
                mouse_pos.y,
                viewport_state.pos,
                viewport_state.size,
                camera_.inv_view_projection());

            const bool shift = ImGui::GetIO().KeyShift;
            const bool ctrl = ImGui::GetIO().KeyCtrl;
            sbox::editor::EditorMode* mode = editor_state_.active_mode();

            if (mode != nullptr) {
                if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                    if (editor_state_.current_mode == ui::EditorState::Mode::Select) {
                        const sbox::editor::PickResult hit = sbox::editor::pick(ray, current_molecule_);
                        editor_consuming_left_drag_ = (hit.type != sbox::editor::PickResult::Type::None);
                    } else {
                        editor_consuming_left_drag_ = true;
                    }
                    mode->on_mouse_down(ray, 0, shift, current_molecule_, editor_state_.selection, editor_state_.commands);
                }
                if (ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
                    mode->on_mouse_down(ray, 1, shift, current_molecule_, editor_state_.selection, editor_state_.commands);
                }

                if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                    mode->on_mouse_up(ray, 0, current_molecule_, editor_state_.selection, editor_state_.commands);
                    editor_consuming_left_drag_ = false;
                }
                if (ImGui::IsMouseReleased(ImGuiMouseButton_Right)) {
                    mode->on_mouse_up(ray, 1, current_molecule_, editor_state_.selection, editor_state_.commands);
                }

                const bool dragging = ImGui::IsMouseDown(ImGuiMouseButton_Left);
                const ImVec2 delta = ImGui::GetIO().MouseDelta;
                mode->on_mouse_move(ray, delta.x, delta.y, dragging, current_molecule_, editor_state_.selection, editor_state_.commands);

                struct KeyMap { ImGuiKey imgui; int glfw; };
                const std::array<KeyMap, 11> keys = {{
                    {ImGuiKey_A, GLFW_KEY_A},
                    {ImGuiKey_Escape, GLFW_KEY_ESCAPE},
                    {ImGuiKey_Delete, GLFW_KEY_DELETE},
                    {ImGuiKey_Backspace, GLFW_KEY_BACKSPACE},
                    {ImGuiKey_Z, GLFW_KEY_Z},
                    {ImGuiKey_Y, GLFW_KEY_Y},
                    {ImGuiKey_1, GLFW_KEY_1},
                    {ImGuiKey_6, GLFW_KEY_6},
                    {ImGuiKey_7, GLFW_KEY_7},
                    {ImGuiKey_8, GLFW_KEY_8},
                    {ImGuiKey_9, GLFW_KEY_9},
                }};
                for (const KeyMap& key : keys) {
                    if (ImGui::IsKeyPressed(key.imgui)) {
                        mode->on_key(key.glfw, ctrl, shift, current_molecule_, editor_state_.selection, editor_state_.commands);
                    }
                }
            }

            if (ctrl && ImGui::IsKeyPressed(ImGuiKey_Z)) {
                if (shift) {
                    editor_state_.commands.redo(current_molecule_);
                } else {
                    editor_state_.commands.undo(current_molecule_);
                }
            }
            if (ctrl && ImGui::IsKeyPressed(ImGuiKey_Y)) {
                editor_state_.commands.redo(current_molecule_);
            }
            if (ctrl && ImGui::IsKeyPressed(ImGuiKey_A)) {
                editor_state_.selection.clear();
                for (int i = 0; i < current_molecule_.num_atoms(); ++i) {
                    editor_state_.selection.toggle_atom(i);
                }
            }
        }

        {
            static std::string last_editor_signature;
            const std::string editor_signature =
                std::to_string(editor_state_.commands.size()) + "|" +
                editor_state_.commands.undo_description() + "|" +
                editor_state_.commands.redo_description() + "|" +
                std::to_string(current_molecule_.num_atoms()) + "|" +
                std::to_string(current_molecule_.num_bonds());
            if (editor_signature != last_editor_signature) {
                last_editor_signature = editor_signature;
                uploadCurrentMoleculeToRenderers();
                state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);
                state_.molecule_loaded = current_molecule_.num_atoms() > 0;
            }
        }

        static int last_color_mode = -1;
        if (last_color_mode != state_.color_mode) {
            last_color_mode = state_.color_mode;
            uploadCurrentMoleculeToRenderers();
        }

        {
            static int last_pes_job_id = -1;
            static int last_pes_point = -1;
            if (state_.pes.scan_complete && !state_.pes.result.geometries.empty()) {
                const int point = std::clamp(
                    state_.pes.viewed_point < 0 ? 0 : state_.pes.viewed_point,
                    0,
                    static_cast<int>(state_.pes.result.geometries.size()) - 1);
                if (state_.pes.scan_job_id != last_pes_job_id || point != last_pes_point) {
                    last_pes_job_id = state_.pes.scan_job_id;
                    last_pes_point = point;
                    state_.pes.viewed_point = point;
                    current_molecule_ = state_.pes.result.geometries[static_cast<std::size_t>(point)];
                    uploadCurrentMoleculeToRenderers();
                    state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);
                    state_.molecule_loaded = current_molecule_.num_atoms() > 0;
                }
            } else {
                last_pes_job_id = -1;
                last_pes_point = -1;
            }
        }

        if (state_.view_mode == ui::ViewMode::MolecularOrbital) {
            sbox::editor::EditorMode* mode = editor_state_.active_mode();
            if (mode != nullptr) {
                ImDrawList* draw_list = ImGui::GetForegroundDrawList();
                const Eigen::Matrix4f vp = camera_.projectionMatrix() * camera_.viewMatrix();
                mode->draw_overlay(draw_list,
                                   current_molecule_,
                                   editor_state_.selection,
                                   vp,
                                   viewport_state.pos,
                                   viewport_state.size);
            }
            ui::draw_context_menu(editor_state_.context_menu,
                                  current_molecule_,
                                  editor_state_.selection,
                                  editor_state_.commands);
            if (settings_manager_.settings().show_atom_labels) {
                const Eigen::Matrix4f vp = camera_.projectionMatrix() * camera_.viewMatrix();
                mol_renderer_.render_atom_labels(current_molecule_,
                                                 vp,
                                                 viewport_state.pos,
                                                 viewport_state.size,
                                                 true,
                                                 true,
                                                 settings_manager_.settings().show_hydrogen_labels);
            }
            if (editor_state_.context_menu.center_view_requested) {
                editor_state_.context_menu.center_view_requested = false;
                Eigen::Vector3f target = Eigen::Vector3f::Zero();
                if (current_molecule_.num_atoms() > 0) {
                    target = current_molecule_.center_of_mass().cast<float>();
                }
                camera_.setTarget(target);
            }
            if (editor_state_.context_menu.fit_view_requested) {
                editor_state_.context_menu.fit_view_requested = false;
                Eigen::Vector3f target = Eigen::Vector3f::Zero();
                if (current_molecule_.num_atoms() > 0) {
                    target = current_molecule_.center_of_mass().cast<float>();
                }
                camera_.setTarget(target);
                camera_.setDistance(std::max(10.0f, state_.mol_bound_radius * 1.8f));
            }
        }

        if ((state_.show_charges || state_.show_bond_orders || state_.show_dipole) &&
            state_.computation.active_job_id >= 0) {
            if (const sbox::backend::JobResult* job_result = backend_.result(state_.computation.active_job_id)) {
                const Eigen::Matrix4f view = camera_.viewMatrix();
                const Eigen::Matrix4f proj = camera_.projectionMatrix();
                if (state_.show_charges && !job_result->mulliken_charges.empty()) {
                    ui::draw_charge_labels(current_molecule_,
                                           job_result->mulliken_charges,
                                           view,
                                           proj,
                                           viewport_state.pos,
                                           viewport_state.size);
                }
                if (state_.show_bond_orders && job_result->mayer_bond_orders.size() > 0) {
                    ui::draw_bond_order_labels(current_molecule_,
                                               job_result->mayer_bond_orders,
                                               view,
                                               proj,
                                               viewport_state.pos,
                                               viewport_state.size);
                }
                if (state_.show_dipole && job_result->dipole_moment.norm() > 1.0e-6) {
                    ui::draw_dipole_arrow(current_molecule_,
                                          job_result->dipole_moment,
                                          view,
                                          proj,
                                          viewport_state.pos,
                                          viewport_state.size);
                }
            }
        }
        if (!state_.constraints.empty()) {
            const Eigen::Matrix4f vp = camera_.projectionMatrix() * camera_.viewMatrix();
            ui::draw_constraint_overlays(state_.constraints, current_molecule_, vp, viewport_state.pos, viewport_state.size);
        }
        draw_color_legend(state_, current_pdb_data_, latest_result_, viewport_state.pos);
        if (settings_manager_.settings().show_status_bar) {
            bool update_clicked = false;
            ui::draw_status_bar(state_, settings_manager_.settings().show_fps, pending_update_ ? &*pending_update_ : nullptr, &update_clicked);
            if (update_clicked) {
                show_update_dialog_ = true;
            }
        }
        const int new_width = std::max(1, static_cast<int>(std::floor(viewport_state.size.x)));
        const int new_height = std::max(1, static_cast<int>(std::floor(viewport_state.size.y)));
        ensureViewportTarget(new_width, new_height);

        camera_.setViewportSize(static_cast<float>(viewport_width_), static_cast<float>(viewport_height_));

        const bool left_drag = glfwGetMouseButton(window_->handle(), GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        const bool middle_drag = glfwGetMouseButton(window_->handle(), GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS;

        const ImVec2 mouse_delta = ImGui::GetIO().MouseDelta;
        bool allow_camera_rotation = true;
        if (state_.view_mode == ui::ViewMode::MolecularOrbital) {
            if (editor_state_.current_mode != ui::EditorState::Mode::Select || editor_consuming_left_drag_) {
                allow_camera_rotation = false;
            }
        }
        camera_.handleRotationDrag(mouse_delta.x, mouse_delta.y, viewport_state.hovered && left_drag && allow_camera_rotation);
        camera_.handlePanDrag(mouse_delta.x, mouse_delta.y, viewport_state.hovered && middle_drag);
        camera_.handleScroll(scroll_delta_, viewport_state.hovered);
        scroll_delta_ = 0.0f;

        if (state_.view_mode == ui::ViewMode::AtomicOrbital && state_.needs_update) {
            state_.update();
        }
        updateMaxDensityEstimate();

        renderViewportToTexture();

        ImGui::Render();

        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0, 0, window_->framebufferWidth(), window_->framebufferHeight());
        glClearColor(0.07f, 0.08f, 0.10f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        window_->swapBuffers();
    }
}

void App::ScrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    (void)xoffset;
    auto* self = static_cast<App*>(glfwGetWindowUserPointer(window));
    if (self != nullptr) {
        self->scroll_delta_ += static_cast<float>(yoffset);
    }
}

void App::initImGui() {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;

    ImGuiStyle& style = ImGui::GetStyle();
    ImGui::StyleColorsDark(&style);
    style.WindowRounding = 4.0f;
    style.FrameRounding = 3.0f;
    style.GrabRounding = 3.0f;
    style.FramePadding = ImVec2(6.0f, 4.0f);
    style.ChildBorderSize = 1.0f;
    style.ScrollbarSize = 12.0f;
    style.ScrollbarRounding = 6.0f;

    ImVec4* colors = style.Colors;
    const ImVec4 accent_active(0.15f, 0.55f, 0.65f, 1.0f);
    const ImVec4 accent(0.13f, 0.44f, 0.52f, 1.0f);
    const ImVec4 accent_hovered(0.18f, 0.62f, 0.74f, 1.0f);

    colors[ImGuiCol_WindowBg] = ImVec4(0.10f, 0.10f, 0.14f, 1.0f);
    colors[ImGuiCol_TitleBg] = ImVec4(0.08f, 0.08f, 0.12f, 1.0f);
    colors[ImGuiCol_TitleBgActive] = ImVec4(0.12f, 0.12f, 0.18f, 1.0f);
    colors[ImGuiCol_FrameBgActive] = accent_active;
    colors[ImGuiCol_ButtonActive] = accent_active;
    colors[ImGuiCol_HeaderActive] = accent_active;
    colors[ImGuiCol_Button] = accent;
    colors[ImGuiCol_ButtonHovered] = accent_hovered;
    colors[ImGuiCol_Header] = accent;
    colors[ImGuiCol_HeaderHovered] = accent_hovered;
    colors[ImGuiCol_CheckMark] = accent_hovered;
    colors[ImGuiCol_SliderGrab] = accent;
    colors[ImGuiCol_SliderGrabActive] = accent_active;
    colors[ImGuiCol_ResizeGrip] = accent;
    colors[ImGuiCol_ResizeGripHovered] = accent_hovered;
    colors[ImGuiCol_ResizeGripActive] = accent_active;
    colors[ImGuiCol_Tab] = ImVec4(0.10f, 0.26f, 0.30f, 1.0f);
    colors[ImGuiCol_TabHovered] = accent_hovered;
    colors[ImGuiCol_TabActive] = accent_active;
    colors[ImGuiCol_TabUnfocused] = ImVec4(0.08f, 0.18f, 0.22f, 1.0f);
    colors[ImGuiCol_TabUnfocusedActive] = ImVec4(0.12f, 0.38f, 0.45f, 1.0f);
    ui::setup_dark_plot_style();

    io.Fonts->Clear();
    const float kFontSize = settings_manager_.settings().font_size;
    const std::array<const char*, 8> font_candidates = {
        "/Library/Fonts/Inter-Regular.ttf",
        "/System/Library/Fonts/Supplemental/Inter.ttc",
        "/System/Library/Fonts/Inter.ttc",
        "/usr/share/fonts/truetype/inter/Inter-Regular.ttf",
        "/usr/share/fonts/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/local/share/fonts/Inter-Regular.ttf",
        "/usr/local/share/fonts/dejavu/DejaVuSans.ttf",
    };

    ImFont* ui_font = nullptr;
    for (const char* path : font_candidates) {
        if (std::filesystem::exists(path)) {
            ui_font = io.Fonts->AddFontFromFileTTF(path, kFontSize);
            if (ui_font != nullptr) {
                break;
            }
        }
    }

    if (ui_font == nullptr) {
        ImFontConfig font_cfg;
        font_cfg.SizePixels = kFontSize;
        ui_font = io.Fonts->AddFontDefault(&font_cfg);
    }
    io.FontDefault = ui_font;

    if (!ImGui_ImplGlfw_InitForOpenGL(window_->handle(), true)) {
        throw std::runtime_error("Failed to initialize ImGui GLFW backend");
    }
    if (!ImGui_ImplOpenGL3_Init("#version 410")) {
        throw std::runtime_error("Failed to initialize ImGui OpenGL3 backend");
    }
}

void App::shutdownImGui() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
}

void App::ensureViewportTarget(int width, int height) {
    width = std::max(width, 1);
    height = std::max(height, 1);

    if (viewport_fbo_ != 0U && width == viewport_width_ && height == viewport_height_) {
        return;
    }

    viewport_width_ = width;
    viewport_height_ = height;

    if (viewport_fbo_ == 0U) {
        glGenFramebuffers(1, &viewport_fbo_);
    }
    if (viewport_color_tex_ == 0U) {
        glGenTextures(1, &viewport_color_tex_);
    }
    if (viewport_depth_rbo_ == 0U) {
        glGenRenderbuffers(1, &viewport_depth_rbo_);
    }

    glBindTexture(GL_TEXTURE_2D, viewport_color_tex_);
    glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_RGBA8,
                 viewport_width_,
                 viewport_height_,
                 0,
                 GL_RGBA,
                 GL_UNSIGNED_BYTE,
                 nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindRenderbuffer(GL_RENDERBUFFER, viewport_depth_rbo_);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, viewport_width_, viewport_height_);

    glBindFramebuffer(GL_FRAMEBUFFER, viewport_fbo_);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, viewport_color_tex_, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, viewport_depth_rbo_);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        throw std::runtime_error("Viewport framebuffer is incomplete");
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

int App::find_homo_index() const {
    int homo = -1;
    for (int i = 0; i < current_mo_data_.occupations.size(); ++i) {
        if (current_mo_data_.occupations(i) > 0.5) {
            homo = i;
        }
    }
    if (homo < 0 && current_mo_data_.coefficients.cols() > 0) {
        homo = 0;
    }
    return homo;
}

float App::compute_mol_bound_radius(const sbox::chem::MolecularSystem& mol) const {
    float max_r = 5.0f;
    const Eigen::Vector3d com = mol.num_atoms() > 0 ? mol.center_of_mass() : Eigen::Vector3d::Zero();
    for (int i = 0; i < mol.num_atoms(); ++i) {
        const float r = static_cast<float>((mol.atom(i).position - com).norm()) + 10.0f;
        max_r = std::max(max_r, r);
    }
    return max_r;
}

void App::apply_settings() {
    sbox::Settings& settings = settings_manager_.settings();

    state_.iso_value = settings.default_iso_value;
    state_.gamma = settings.default_gamma;
    state_.mol_render_mode = settings.mol_render_mode;
    state_.color_mode = settings.color_mode;
    state_.computation.method = static_cast<sbox::backend::Method>(
        std::clamp(settings.default_method, 0, static_cast<int>(sbox::backend::Method::GFN_FF)));
    state_.computation.basis = static_cast<sbox::backend::BasisSetType>(
        std::clamp(settings.default_basis, 0, static_cast<int>(sbox::backend::BasisSetType::aug_cc_pVTZ)));
    state_.computation.charge = settings.default_charge;
    state_.computation.multiplicity = std::max(1, settings.default_multiplicity);

    sbox::render::set_atom_radius_scale(settings.atom_scale);
    sbox::render::set_bond_radius_scale(settings.bond_scale);
    glfwSwapInterval(settings.enable_vsync ? 1 : 0);

    if (settings.python_auto_detect) {
        python_env_.detect();
        if (python_env_.is_valid()) {
            settings.python_path = python_env_.info().python_path;
        }
    } else if (!settings.python_path.empty()) {
        python_env_.set_python_path(settings.python_path);
    }

    if (python_env_.is_valid()) {
        python_env_.check_packages();
        backend_.init(python_env_);
    }

    if (current_molecule_.num_atoms() > 0) {
        uploadCurrentMoleculeToRenderers();
    }
    settings_manager_.save();
}

void App::applyMOData(const sbox::basis::MOData& mo_data, const std::string& name_hint) {
    current_mo_data_ = mo_data;
    current_molecule_.clear();
    current_pdb_data_ = sbox::io::PDBData{};
    current_molecule_.set_name(name_hint);

    const std::size_t atom_count = std::min(current_mo_data_.atom_positions.size(), current_mo_data_.atomic_numbers.size());
    for (std::size_t i = 0; i < atom_count; ++i) {
        current_molecule_.add_atom({current_mo_data_.atomic_numbers[i], current_mo_data_.atom_positions[i], "", 0});
    }
    current_molecule_.set_charge(state_.computation.charge);
    current_molecule_.set_multiplicity(state_.computation.multiplicity);
    current_molecule_.perceive_bonds();
    uploadCurrentMoleculeToRenderers();

    use_cube_fallback_ = !basis_textures_.upload(current_mo_data_);
    has_mo_data_ = !use_cube_fallback_;
    has_cube_data_ = false;
    state_.molecule_loaded = current_molecule_.num_atoms() > 0;
    state_.mol_has_mo_summary = true;
    state_.mol_num_basis = current_mo_data_.basis.num_basis_functions();
    state_.mol_total_energy_h = current_mo_data_.total_energy;

    state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);

    state_.view_mode = ui::ViewMode::MolecularOrbital;
    state_.num_mo = static_cast<int>(current_mo_data_.coefficients.cols());
    state_.homo_index = find_homo_index();
    state_.selected_mo = state_.homo_index;
    if (state_.homo_index >= 0 && state_.homo_index + 1 < current_mo_data_.energies.size()) {
        state_.mol_homo_lumo_gap_ev =
            (current_mo_data_.energies(state_.homo_index + 1) - current_mo_data_.energies(state_.homo_index)) * 27.2114;
    } else {
        state_.mol_homo_lumo_gap_ev = 0.0;
    }
}

void App::applyBackendResult(const sbox::backend::JobResult& result) {
    if (result.has_optimized_geometry) {
        current_molecule_ = result.optimized_geometry;
        current_pdb_data_ = sbox::io::PDBData{};
        uploadCurrentMoleculeToRenderers();
        state_.molecule_loaded = current_molecule_.num_atoms() > 0;
        state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);
    }

    if (result.has_mo_data) {
        const std::string name_hint = current_molecule_.name().empty() ? "Backend Result" : current_molecule_.name();
        applyMOData(result.mo_data, name_hint);
        if (result.has_density_cube && result.has_esp_cube) {
            loadESPSurface(result);
        }
        return;
    }

    if (result.has_homo_cube) {
        if (!volume_texture_.upload(result.homo_cube)) {
            throw std::runtime_error("Failed to upload HOMO cube result");
        }
        current_molecule_.clear();
        for (std::size_t i = 0; i < result.homo_cube.atom_Z.size() && i < result.homo_cube.atom_pos.size(); ++i) {
            current_molecule_.add_atom({result.homo_cube.atom_Z[i], result.homo_cube.atom_pos[i], "", 0});
        }
        current_molecule_.perceive_bonds();
        current_pdb_data_ = sbox::io::PDBData{};
        uploadCurrentMoleculeToRenderers();
        has_cube_data_ = true;
        has_mo_data_ = false;
        use_cube_fallback_ = true;
        state_.molecule_loaded = current_molecule_.num_atoms() > 0;
        state_.view_mode = ui::ViewMode::MolecularOrbital;
        state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);
        clear_mo_summary(state_);
        return;
    }

    if (result.has_density_cube) {
        if (!volume_texture_.upload(result.density_cube)) {
            throw std::runtime_error("Failed to upload density cube result");
        }
        current_molecule_.clear();
        for (std::size_t i = 0; i < result.density_cube.atom_Z.size() && i < result.density_cube.atom_pos.size(); ++i) {
            current_molecule_.add_atom({result.density_cube.atom_Z[i], result.density_cube.atom_pos[i], "", 0});
        }
        current_molecule_.perceive_bonds();
        current_pdb_data_ = sbox::io::PDBData{};
        uploadCurrentMoleculeToRenderers();
        has_cube_data_ = true;
        has_mo_data_ = false;
        use_cube_fallback_ = true;
        state_.molecule_loaded = current_molecule_.num_atoms() > 0;
        state_.view_mode = ui::ViewMode::MolecularOrbital;
        state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);
        clear_mo_summary(state_);
        if (result.has_esp_cube) {
            loadESPSurface(result);
        }
    }
}

void App::loadESPSurface(const sbox::backend::JobResult& result) {
    if (result.has_density_cube && result.has_esp_cube) {
        esp_surface_.upload(result.density_cube, result.esp_cube);
    }
}

void App::detect_metal_center() {
    detected_metal_choices_.clear();
    for (int i = 0; i < current_molecule_.num_atoms(); ++i) {
        if (is_transition_metal_Z(current_molecule_.atom(i).Z)) {
            detected_metal_choices_.push_back(i);
        }
    }

    if (detected_metal_choices_.empty()) {
        current_metal_index_ = -1;
        message_popup_title_ = "Metal Detection";
        message_popup_text_ = "No transition metal found in the molecule.";
        open_message_popup_ = true;
        return;
    }

    current_metal_index_ = detected_metal_choices_.front();
    if (detected_metal_choices_.size() > 1) {
        open_metal_picker_popup_ = true;
    }
}

void App::run_crystal_field_analysis() {
    if (current_molecule_.num_atoms() == 0) {
        has_d_orbital_analysis_ = false;
        message_popup_title_ = "Crystal Field Analysis";
        message_popup_text_ = "No molecule loaded.";
        open_message_popup_ = true;
        return;
    }

    if (current_metal_index_ < 0 || current_metal_index_ >= current_molecule_.num_atoms() ||
        !is_transition_metal_Z(current_molecule_.atom(current_metal_index_).Z)) {
        detect_metal_center();
        if (current_metal_index_ < 0) {
            message_popup_title_ = "Crystal Field Analysis";
            message_popup_text_ = "No transition metal found in the molecule.";
            open_message_popup_ = true;
            return;
        }
    }

    if ((!latest_result_.has_value() || !latest_result_->has_mo_data) && current_mo_data_.coefficients.cols() == 0) {
        has_d_orbital_analysis_ = false;
        message_popup_title_ = "Crystal Field Analysis";
        message_popup_text_ = "Run a calculation first (HF, DFT, or xTB).";
        open_message_popup_ = true;
        return;
    }

    const sbox::basis::MOData& mo_data =
        (latest_result_.has_value() && latest_result_->has_mo_data) ? latest_result_->mo_data : current_mo_data_;
    current_d_orbitals_ = sbox::analysis::extract_d_orbitals(mo_data, current_molecule_, current_metal_index_);
    const sbox::chem::CoordinationGeometry geometry =
        sbox::chem::detect_coordination_geometry(current_molecule_, current_metal_index_);
    sbox::analysis::identify_splitting(current_d_orbitals_, geometry);
    has_d_orbital_analysis_ = true;
    state_.selected_crystal_field_metal = 0;
}

sbox::backend::JobSpec App::makeJobSpecFromState() const {
    if (current_molecule_.num_atoms() == 0) {
        throw std::runtime_error("No molecule loaded");
    }

    sbox::backend::JobSpec spec;
    spec.geometry = current_molecule_;
    spec.geometry.set_charge(state_.computation.charge);
    spec.geometry.set_multiplicity(state_.computation.multiplicity);
    spec.method = state_.computation.method;
    spec.basis = state_.computation.basis;
    spec.charge = state_.computation.charge;
    spec.multiplicity = state_.computation.multiplicity;
    spec.optimize_geometry = state_.computation.optimize;
    spec.solvent = state_.computation.solvent;
    spec.properties = {
        sbox::backend::PropertyRequest::MullikenCharges,
        sbox::backend::PropertyRequest::DipoleMoment,
    };

    if (sbox::backend::method_is_xtb(spec.method)) {
        spec.properties.push_back(sbox::backend::PropertyRequest::MoldenFile);
    } else {
        spec.properties.push_back(sbox::backend::PropertyRequest::MoldenFile);
        spec.properties.push_back(sbox::backend::PropertyRequest::CubeHOMO);
        spec.properties.push_back(sbox::backend::PropertyRequest::CubeLUMO);
    }

    if (spec.optimize_geometry) {
        spec.properties.push_back(sbox::backend::PropertyRequest::Optimization);
    }
    for (const auto& constraint : state_.constraints) {
        if (!constraint.active) {
            continue;
        }
        switch (constraint.type) {
        case ui::AppState::GeometricConstraint::Type::FreezeAtom:
            if (!constraint.atom_indices.empty()) {
                spec.constraints.freeze_atoms.push_back(constraint.atom_indices[0]);
            }
            break;
        case ui::AppState::GeometricConstraint::Type::FixDistance:
            if (constraint.atom_indices.size() >= 2) {
                spec.constraints.fixed_distances.emplace_back(
                    constraint.atom_indices[0], constraint.atom_indices[1], constraint.value);
            }
            break;
        case ui::AppState::GeometricConstraint::Type::FixAngle:
            if (constraint.atom_indices.size() >= 3) {
                spec.constraints.fixed_angles.emplace_back(
                    constraint.atom_indices[0], constraint.atom_indices[1], constraint.atom_indices[2], constraint.value);
            }
            break;
        case ui::AppState::GeometricConstraint::Type::FixDihedral:
            if (constraint.atom_indices.size() >= 4) {
                spec.constraints.fixed_dihedrals.emplace_back(
                    constraint.atom_indices[0], constraint.atom_indices[1], constraint.atom_indices[2], constraint.atom_indices[3], constraint.value);
            }
            break;
        }
    }
    return spec;
}

void App::loadMoldenFile(const std::string& path) {
    applyMOData(sbox::molden::parse_molden_file(path), std::filesystem::path(path).filename().string());
    state_.computation.charge = current_molecule_.charge();
    state_.computation.multiplicity = current_molecule_.multiplicity();
    settings_manager_.settings().last_open_directory = std::filesystem::path(path).parent_path().string();
    settings_manager_.add_recent_file(path);
}

void App::loadCubeFile(const std::string& path) {
    const sbox::io::CubeData cube = sbox::io::read_cube(path);
    current_molecule_.clear();
    current_pdb_data_ = sbox::io::PDBData{};
    current_molecule_.set_name(std::filesystem::path(path).filename().string());
    for (std::size_t i = 0; i < cube.atom_Z.size() && i < cube.atom_pos.size(); ++i) {
        current_molecule_.add_atom({cube.atom_Z[i], cube.atom_pos[i], "", 0});
    }
    current_molecule_.perceive_bonds();
    uploadCurrentMoleculeToRenderers();

    if (!volume_texture_.upload(cube)) {
        throw std::runtime_error("Failed to upload cube volume texture");
    }

    has_cube_data_ = true;
    use_cube_fallback_ = true;
    has_mo_data_ = false;
    state_.molecule_loaded = current_molecule_.num_atoms() > 0;
    state_.view_mode = ui::ViewMode::MolecularOrbital;
    state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);
    clear_mo_summary(state_);
    state_.computation.charge = current_molecule_.charge();
    state_.computation.multiplicity = current_molecule_.multiplicity();
    settings_manager_.settings().last_open_directory = std::filesystem::path(path).parent_path().string();
    settings_manager_.add_recent_file(path);
}

void App::loadXYZFile(const std::string& path) {
    current_trajectory_ = sbox::io::Trajectory{};
    has_trajectory_ = false;
    current_molecule_ = sbox::io::read_xyz(path);
    current_pdb_data_ = sbox::io::PDBData{};
    uploadCurrentMoleculeToRenderers();
    current_mo_data_ = sbox::basis::MOData{};
    has_mo_data_ = false;
    has_cube_data_ = false;
    use_cube_fallback_ = false;
    state_.molecule_loaded = current_molecule_.num_atoms() > 0;
    state_.view_mode = ui::ViewMode::MolecularOrbital;
    state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);
    clear_mo_summary(state_);
    state_.computation.charge = current_molecule_.charge();
    state_.computation.multiplicity = current_molecule_.multiplicity();
    settings_manager_.settings().last_open_directory = std::filesystem::path(path).parent_path().string();
    settings_manager_.add_recent_file(path);
}

void App::loadTrajectoryFile(const std::string& path) {
    current_trajectory_ = sbox::io::read_trajectory_xyz(path);
    has_trajectory_ = !current_trajectory_.empty();
    if (!has_trajectory_) {
        throw std::runtime_error("Trajectory file did not contain any frames");
    }

    current_molecule_ = current_trajectory_.frames.front().geometry;
    current_molecule_.set_name(std::filesystem::path(path).filename().string());
    current_pdb_data_ = sbox::io::PDBData{};
    uploadCurrentMoleculeToRenderers();
    state_.view_mode = ui::ViewMode::MolecularOrbital;
    has_mo_data_ = false;
    has_cube_data_ = false;
    use_cube_fallback_ = false;
    clear_mo_summary(state_);
    state_.molecule_loaded = current_molecule_.num_atoms() > 0;
    state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);
    state_.computation.charge = current_molecule_.charge();
    state_.computation.multiplicity = current_molecule_.multiplicity();
    state_.optimization_player = {};
    state_.optimization_player.total_frames = current_trajectory_.num_frames();
    settings_manager_.settings().last_open_directory = std::filesystem::path(path).parent_path().string();
    settings_manager_.add_recent_file(path);
}

void App::loadSDFFile(const std::string& path) {
    current_molecule_ = sbox::io::read_sdf(path);
    current_pdb_data_ = sbox::io::PDBData{};
    uploadCurrentMoleculeToRenderers();
    current_mo_data_ = sbox::basis::MOData{};
    has_mo_data_ = false;
    has_cube_data_ = false;
    use_cube_fallback_ = false;
    state_.molecule_loaded = current_molecule_.num_atoms() > 0;
    state_.view_mode = ui::ViewMode::MolecularOrbital;
    state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);
    clear_mo_summary(state_);
    state_.computation.charge = current_molecule_.charge();
    state_.computation.multiplicity = current_molecule_.multiplicity();
    settings_manager_.settings().last_open_directory = std::filesystem::path(path).parent_path().string();
    settings_manager_.add_recent_file(path);
}

void App::loadPDBFile(const std::string& path) {
    current_pdb_data_ = sbox::io::read_pdb(path);
    current_molecule_ = current_pdb_data_.to_molecular_system();
    current_molecule_.set_name(current_pdb_data_.title.empty() ? std::filesystem::path(path).filename().string() : current_pdb_data_.title);
    uploadCurrentMoleculeToRenderers();
    current_mo_data_ = sbox::basis::MOData{};
    has_mo_data_ = false;
    has_cube_data_ = false;
    use_cube_fallback_ = false;
    state_.molecule_loaded = current_molecule_.num_atoms() > 0;
    state_.view_mode = ui::ViewMode::MolecularOrbital;
    state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);
    clear_mo_summary(state_);
    state_.computation.charge = current_molecule_.charge();
    state_.computation.multiplicity = current_molecule_.multiplicity();
    settings_manager_.settings().last_open_directory = std::filesystem::path(path).parent_path().string();
    settings_manager_.add_recent_file(path);
}

void App::loadFchkFile(const std::string& path) {
    const sbox::io::FchkData fchk = sbox::io::read_fchk(path);
    current_molecule_.clear();
    current_pdb_data_ = sbox::io::PDBData{};
    current_molecule_.set_name(fchk.title);
    current_molecule_.set_charge(fchk.charge);
    current_molecule_.set_multiplicity(fchk.multiplicity);
    for (int i = 0; i < fchk.num_atoms; ++i) {
        current_molecule_.add_atom({fchk.atomic_numbers[static_cast<std::size_t>(i)],
                                    Eigen::Vector3d(fchk.coordinates[static_cast<std::size_t>(3 * i + 0)],
                                                    fchk.coordinates[static_cast<std::size_t>(3 * i + 1)],
                                                    fchk.coordinates[static_cast<std::size_t>(3 * i + 2)]),
                                    "",
                                    0});
    }
    current_molecule_.perceive_bonds();
    uploadCurrentMoleculeToRenderers();

    current_mo_data_ = sbox::basis::MOData{};
    const std::optional<sbox::basis::MOData> maybe_mo = mo_data_from_fchk(fchk);
    if (maybe_mo.has_value()) {
        current_mo_data_ = *maybe_mo;
        use_cube_fallback_ = !basis_textures_.upload(current_mo_data_);
        has_mo_data_ = !use_cube_fallback_;
    } else {
        has_mo_data_ = false;
        use_cube_fallback_ = false;
    }
    has_cube_data_ = false;
    state_.view_mode = ui::ViewMode::MolecularOrbital;
    state_.molecule_loaded = current_molecule_.num_atoms() > 0;
    state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);
    state_.num_mo = static_cast<int>(fchk.mo_energies.size());
    state_.homo_index = -1;
    for (int i = 0; i < fchk.occupations.size(); ++i) {
        if (fchk.occupations(i) > 0.5) {
            state_.homo_index = i;
        }
    }
    state_.selected_mo = state_.homo_index;
    state_.mol_num_basis = fchk.num_basis;
    state_.mol_total_energy_h = fchk.total_energy;
    state_.mol_has_mo_summary = true;
    if (state_.homo_index >= 0 && state_.homo_index + 1 < fchk.mo_energies.size()) {
        state_.mol_homo_lumo_gap_ev =
            (fchk.mo_energies(state_.homo_index + 1) - fchk.mo_energies(state_.homo_index)) * 27.2114;
    } else {
        state_.mol_homo_lumo_gap_ev = 0.0;
    }
    state_.computation.charge = current_molecule_.charge();
    state_.computation.multiplicity = current_molecule_.multiplicity();
    settings_manager_.settings().last_open_directory = std::filesystem::path(path).parent_path().string();
    settings_manager_.add_recent_file(path);
}

void App::loadProjectFile(const std::string& path) {
    nlohmann::json extra_state;
    current_molecule_ = sbox::io::load_project(path, &extra_state);
    current_pdb_data_ = sbox::io::PDBData{};
    current_trajectory_ = sbox::io::Trajectory{};
    has_trajectory_ = false;
    uploadCurrentMoleculeToRenderers();
    current_mo_data_ = sbox::basis::MOData{};
    has_mo_data_ = false;
    has_cube_data_ = false;
    use_cube_fallback_ = false;
    clear_mo_summary(state_);
    state_.molecule_loaded = current_molecule_.num_atoms() > 0;
    state_.view_mode = ui::ViewMode::MolecularOrbital;
    state_.mol_bound_radius = compute_mol_bound_radius(current_molecule_);
    state_.computation.charge = current_molecule_.charge();
    state_.computation.multiplicity = current_molecule_.multiplicity();
    settings_manager_.settings().last_open_directory = std::filesystem::path(path).parent_path().string();
    settings_manager_.add_recent_file(path);
}

void App::load_file_by_extension(const std::string& path) {
    std::string ext = std::filesystem::path(path).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    if (ext == ".xyz") {
        loadXYZFile(path);
    } else if (ext == ".sdf" || ext == ".mol") {
        loadSDFFile(path);
    } else if (ext == ".molden") {
        loadMoldenFile(path);
    } else if (ext == ".cube") {
        loadCubeFile(path);
    } else if (ext == ".fchk" || ext == ".fch") {
        loadFchkFile(path);
    } else if (ext == ".pdb" || ext == ".ent") {
        loadPDBFile(path);
    } else if (ext == ".sbox") {
        loadProjectFile(path);
    } else {
        throw std::runtime_error("Unsupported file type: " + ext);
    }
}

float App::computeMaxDensityEstimate() const {
    if (state_.view_mode == ui::ViewMode::MolecularOrbital) {
        if (use_cube_fallback_ && has_cube_data_ && volume_texture_.is_uploaded()) {
            const float max_psi = volume_texture_.max_abs_value();
            return std::max(max_psi * max_psi, 1e-6f);
        }

        if (has_mo_data_ && current_mo_data_.coefficients.cols() > 0) {
            const int mo_index = state_.selected_mo >= 0 ? state_.selected_mo : find_homo_index();
            double max_density = 0.0;
            Eigen::Vector3d center = current_molecule_.num_atoms() > 0 ? current_molecule_.center_of_mass() : Eigen::Vector3d::Zero();
            const std::array<Eigen::Vector3d, 7> offsets = {
                Eigen::Vector3d::Zero(),
                Eigen::Vector3d(1.0, 0.0, 0.0),
                Eigen::Vector3d(-1.0, 0.0, 0.0),
                Eigen::Vector3d(0.0, 1.0, 0.0),
                Eigen::Vector3d(0.0, -1.0, 0.0),
                Eigen::Vector3d(0.0, 0.0, 1.0),
                Eigen::Vector3d(0.0, 0.0, -1.0),
            };

            max_density = std::max(max_density,
                                   sbox::basis::evaluate_mo_density_at_point(current_mo_data_, mo_index, center));
            for (const sbox::chem::Atom& atom : current_molecule_.atoms()) {
                for (const Eigen::Vector3d& offset : offsets) {
                    max_density = std::max(max_density,
                                           sbox::basis::evaluate_mo_density_at_point(current_mo_data_,
                                                                                     mo_index,
                                                                                     atom.position + offset));
                }
            }

            if (max_density <= 0.0 || std::isnan(max_density)) {
                return 1.0f;
            }
            return static_cast<float>(std::max(max_density, 1e-6));
        }

        return 1.0f;
    }

    const int n = state_.current_n;
    const int l = state_.current_l;
    const int m = state_.selected_m;
    const double zeff = static_cast<double>(state_.current_Zeff);

    if (n <= 0 || zeff <= 0.0) {
        return 1.0f;
    }

    double max_density = 0.0;
    if (l == 0) {
        max_density = sbox::hydrogen::probability_density(n, l, m, zeff, 0.0, 0.0, 0.0);
    } else {
        const int samples = 100;
        const double max_r = 4.0 * static_cast<double>(n * n) / zeff;
        for (int i = 0; i < samples; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(samples - 1);
            const double r = max_r * t;
            const double density_z = sbox::hydrogen::probability_density(n, l, m, zeff, 0.0, 0.0, r);
            const double density_x = sbox::hydrogen::probability_density(n, l, m, zeff, r, 0.0, 0.0);
            max_density = std::max(max_density, std::max(density_z, density_x));
        }
    }

    if (max_density <= 0.0 || std::isnan(max_density)) {
        max_density = 1.0;
    } else if (max_density <= 1e-9) {
        max_density = 1e-6;
    }
    return static_cast<float>(max_density);
}

void App::updateMaxDensityEstimate() {
    if (state_.view_mode == ui::ViewMode::MolecularOrbital) {
        max_density_estimate_ = computeMaxDensityEstimate();
        return;
    }

    const bool changed = (density_n_ != state_.current_n) ||
                         (density_l_ != state_.current_l) ||
                         (density_m_ != state_.selected_m) ||
                         (std::abs(density_zeff_ - state_.current_Zeff) > 1e-5f);

    if (!changed) {
        return;
    }

    max_density_estimate_ = computeMaxDensityEstimate();
    density_n_ = state_.current_n;
    density_l_ = state_.current_l;
    density_m_ = state_.selected_m;
    density_zeff_ = state_.current_Zeff;
}

void App::renderViewportToTexture() {
    renderViewportToTarget(viewport_fbo_, viewport_width_, viewport_height_, false);
}

void App::render_single_frame() {
    window_->pollEvents();
    updateMaxDensityEstimate();
    renderViewportToTexture();
}

void App::render_to_fbo(unsigned int fbo, int w, int h) {
    renderViewportToTarget(fbo, w, h, false);
}

ui::AppState& App::state() {
    return state_;
}

const ui::AppState& App::state() const {
    return state_;
}

void App::renderViewportToTarget(unsigned int fbo, int width, int height, bool transparent_background) {
    camera_.setViewportSize(static_cast<float>(width), static_cast<float>(height));
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glViewport(0, 0, width, height);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.04f, 0.055f, 0.09f, transparent_background ? 0.0f : 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (state_.view_mode == ui::ViewMode::AtomicOrbital) {
        Shader* active_shader = orbital_shader_ ? orbital_shader_.get() : gradient_shader_.get();
        if (active_shader == nullptr) {
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            camera_.setViewportSize(static_cast<float>(viewport_width_), static_cast<float>(viewport_height_));
            return;
        }

        active_shader->bind();
        if (orbital_shader_) {
            active_shader->setUniform("u_n", state_.current_n);
            active_shader->setUniform("u_l", state_.current_l);
            active_shader->setUniform("u_m", state_.selected_m);
            active_shader->setUniform("u_Zeff", state_.current_Zeff);
            active_shader->setUniform("u_inv_vp", camera_.inv_view_projection());
            active_shader->setUniform("u_camera_pos", camera_.camera_position());
            active_shader->setUniform("u_iso_value", state_.iso_value);
            active_shader->setUniform("u_render_mode", state_.render_mode);
            active_shader->setUniform("u_gamma", state_.gamma);
            active_shader->setUniform("u_volume_steps", settings_manager_.settings().volume_steps);
            active_shader->setUniform("u_isosurface_steps", settings_manager_.settings().isosurface_steps);

            const int res_loc = glGetUniformLocation(active_shader->id(), "u_resolution");
            if (res_loc >= 0) {
                glUniform2f(res_loc, static_cast<float>(width), static_cast<float>(height));
            }

            const int density_loc = glGetUniformLocation(active_shader->id(), "u_max_density");
            if (density_loc >= 0) {
                glUniform1f(density_loc, max_density_estimate_);
            }
        }

        glBindVertexArray(fullscreen_vao_);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glBindVertexArray(0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        camera_.setViewportSize(static_cast<float>(viewport_width_), static_cast<float>(viewport_height_));
        return;
    }

    state_.lod_atoms_rendered = 0;
    state_.lod_atoms_culled = 0;
    state_.lod_bonds_rendered = 0;

    if (mol_renderer_.has_data()) {
        if (current_molecule_.num_atoms() > settings_manager_.settings().lod_threshold_atoms) {
            lod_renderer_.render(camera_.viewMatrix(),
                                 camera_.projectionMatrix(),
                                 camera_.camera_position(),
                                 static_cast<sbox::render::MolRenderMode>(state_.mol_render_mode),
                                 current_molecule_,
                                 static_cast<sbox::render::ColorMode>(state_.color_mode),
                                 current_pdb_data_.atoms.empty() ? nullptr : &current_pdb_data_,
                                 current_charges_for_render(latest_result_));
            state_.lod_atoms_rendered = lod_renderer_.atoms_rendered();
            state_.lod_atoms_culled = lod_renderer_.atoms_culled();
            state_.lod_bonds_rendered = lod_renderer_.bonds_rendered();
        } else {
            mol_renderer_.render(camera_.viewMatrix(),
                                 camera_.projectionMatrix(),
                                 camera_.camera_position(),
                                 static_cast<sbox::render::MolRenderMode>(state_.mol_render_mode));
            state_.lod_atoms_rendered = current_molecule_.num_atoms();
            state_.lod_atoms_culled = 0;
            state_.lod_bonds_rendered = current_molecule_.num_bonds();
        }
        if (!editor_state_.selection.empty()) {
            mol_renderer_.render_selection(camera_.viewMatrix(),
                                           camera_.projectionMatrix(),
                                           camera_.camera_position(),
                                           current_molecule_,
                                           editor_state_.selection,
                                           static_cast<sbox::render::MolRenderMode>(state_.mol_render_mode));
        }
    }

    Shader* active = nullptr;
    if (state_.show_esp_surface && esp_surface_.is_uploaded() && esp_shader_) {
        active = esp_shader_.get();
    } else if (use_cube_fallback_ && has_cube_data_ && cube_shader_) {
        active = cube_shader_.get();
    } else if (has_mo_data_ && mo_shader_) {
        active = mo_shader_.get();
    }

    if (active != nullptr) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        active->bind();
        active->setUniform("u_inv_vp", camera_.inv_view_projection());
        active->setUniform("u_camera_pos", camera_.camera_position());
        active->setUniform("u_render_mode", state_.render_mode);
        active->setUniform("u_iso_value", active == esp_shader_.get() ? state_.esp_density_iso : state_.iso_value);
        active->setUniform("u_gamma", state_.gamma);
        active->setUniform("u_max_density", max_density_estimate_);
        active->setUniform("u_bound_radius", state_.mol_bound_radius);
        active->setUniform("u_volume_steps", settings_manager_.settings().volume_steps);
        active->setUniform("u_isosurface_steps", settings_manager_.settings().isosurface_steps);

        const int res_loc = glGetUniformLocation(active->id(), "u_resolution");
        if (res_loc >= 0) {
            glUniform2f(res_loc, static_cast<float>(width), static_cast<float>(height));
        }

        if (active == esp_shader_.get()) {
            esp_surface_.bind(active->id(), 6, 7);
            const Eigen::Vector3f orig = esp_surface_.origin();
            const Eigen::Matrix3f w2g = esp_surface_.world_to_grid();
            const int origin_loc = glGetUniformLocation(active->id(), "u_grid_origin");
            if (origin_loc >= 0) {
                glUniform3f(origin_loc, orig.x(), orig.y(), orig.z());
            }
            const int w2g_loc = glGetUniformLocation(active->id(), "u_world_to_grid");
            if (w2g_loc >= 0) {
                glUniformMatrix3fv(w2g_loc, 1, GL_FALSE, w2g.data());
            }
            if (latest_result_.has_value() && latest_result_->has_density_cube) {
                const int dims_loc = glGetUniformLocation(active->id(), "u_grid_dims");
                if (dims_loc >= 0) {
                    glUniform3i(dims_loc,
                                latest_result_->density_cube.nx,
                                latest_result_->density_cube.ny,
                                latest_result_->density_cube.nz);
                }
            }
            const int min_loc = glGetUniformLocation(active->id(), "u_esp_min");
            if (min_loc >= 0) {
                glUniform1f(min_loc, state_.esp_auto_range ? esp_surface_.esp_min() : state_.esp_color_min);
            }
            const int max_loc = glGetUniformLocation(active->id(), "u_esp_max");
            if (max_loc >= 0) {
                glUniform1f(max_loc, state_.esp_auto_range ? esp_surface_.esp_max() : state_.esp_color_max);
            }
        } else if (!use_cube_fallback_) {
            const int mo_idx = state_.selected_mo >= 0 ? state_.selected_mo : find_homo_index();
            active->setUniform("u_mo_index", mo_idx);
            active->setUniform("u_num_shells", basis_textures_.num_shells());
            active->setUniform("u_num_basis", basis_textures_.num_basis());
            active->setUniform("u_num_mo", basis_textures_.num_mo());
            basis_textures_.bind(active->id(), 1);
        } else {
            volume_texture_.bind(active->id(), 5);
            const Eigen::Vector3f orig = volume_texture_.origin();
            const Eigen::Matrix3f w2g = volume_texture_.world_to_grid();
            const int origin_loc = glGetUniformLocation(active->id(), "u_grid_origin");
            if (origin_loc >= 0) {
                glUniform3f(origin_loc, orig.x(), orig.y(), orig.z());
            }
            const int w2g_loc = glGetUniformLocation(active->id(), "u_world_to_grid");
            if (w2g_loc >= 0) {
                glUniformMatrix3fv(w2g_loc, 1, GL_FALSE, w2g.data());
            }
            const int dims_loc = glGetUniformLocation(active->id(), "u_grid_dims");
            if (dims_loc >= 0) {
                glUniform3i(dims_loc, volume_texture_.nx(), volume_texture_.ny(), volume_texture_.nz());
            }
        }

        glBindVertexArray(fullscreen_vao_);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glBindVertexArray(0);

        if (active == esp_shader_.get()) {
            esp_surface_.unbind();
        } else if (!use_cube_fallback_) {
            basis_textures_.unbind();
        } else {
            volume_texture_.unbind();
        }

        glDisable(GL_BLEND);
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    camera_.setViewportSize(static_cast<float>(viewport_width_), static_cast<float>(viewport_height_));
}

void App::uploadCurrentMoleculeToRenderers() {
    has_d_orbital_analysis_ = false;
    current_d_orbitals_ = {};
    current_metal_index_ = -1;
    sbox::render::set_atom_radius_scale(settings_manager_.settings().atom_scale);
    sbox::render::set_bond_radius_scale(settings_manager_.settings().bond_scale);
    const auto color_mode = static_cast<sbox::render::ColorMode>(state_.color_mode);
    const sbox::io::PDBData* pdb_ptr = current_pdb_data_.atoms.empty() ? nullptr : &current_pdb_data_;
    const std::vector<double>* charges = current_charges_for_render(latest_result_);
    mol_renderer_.upload(current_molecule_, color_mode, pdb_ptr, charges);
    lod_renderer_.rebuild(current_molecule_, color_mode, pdb_ptr, charges);
}

}  // namespace sbox
