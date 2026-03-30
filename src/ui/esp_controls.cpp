#include "ui/esp_controls.h"

#include "core/elements.h"

#include <imgui.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>

namespace sbox::ui {

namespace {

constexpr double kHartreeToKcalMol = 627.509;

ImVec4 esp_color(float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    const float v = t * 2.0f - 1.0f;
    if (v < -0.5f) {
        const float u = (v + 1.0f) * 2.0f;
        return ImVec4(0.0f + 0.2f * u, 0.0f + 0.4f * u, 0.6f + 0.4f * u, 1.0f);
    }
    if (v < 0.0f) {
        const float u = (v + 0.5f) * 2.0f;
        return ImVec4(0.2f + 0.7f * u, 0.4f + 0.5f * u, 1.0f - 0.05f * u, 1.0f);
    }
    if (v < 0.5f) {
        const float u = v * 2.0f;
        return ImVec4(0.9f + 0.1f * u, 0.9f - 0.5f * u, 0.95f - 0.75f * u, 1.0f);
    }
    const float u = (v - 0.5f) * 2.0f;
    return ImVec4(1.0f - 0.4f * u, 0.4f * (1.0f - u), 0.2f * (1.0f - u), 1.0f);
}

ImVec4 component_color(double value) {
    return value >= 0.0 ? ImVec4(0.90f, 0.35f, 0.35f, 1.0f) : ImVec4(0.30f, 0.55f, 0.95f, 1.0f);
}

const char* polarity_label(double dipole_mag) {
    if (dipole_mag < 0.5) {
        return "nonpolar";
    }
    if (dipole_mag < 1.5) {
        return "weakly polar";
    }
    if (dipole_mag < 3.0) {
        return "polar";
    }
    return "highly polar";
}

const char* dominant_axis_label(const Eigen::Vector3d& dipole) {
    const std::array<double, 3> mags = {std::abs(dipole.x()), std::abs(dipole.y()), std::abs(dipole.z())};
    const auto it = std::max_element(mags.begin(), mags.end());
    const int idx = static_cast<int>(std::distance(mags.begin(), it));
    switch (idx) {
    case 0: return "X";
    case 1: return "Y";
    default: return "Z";
    }
}

void draw_component_bars(const Eigen::Vector3d& dipole) {
    const float width = ImGui::GetContentRegionAvail().x;
    const float bar_height = 14.0f;
    const float label_width = 28.0f;
    const double max_abs = std::max({std::abs(dipole.x()), std::abs(dipole.y()), std::abs(dipole.z()), 1.0e-6});
    const std::array<const char*, 3> labels = {"mx", "my", "mz"};
    const std::array<double, 3> values = {dipole.x(), dipole.y(), dipole.z()};

    for (int i = 0; i < 3; ++i) {
        ImGui::TextUnformatted(labels[static_cast<std::size_t>(i)]);
        ImGui::SameLine(label_width);
        const ImVec2 pos = ImGui::GetCursorScreenPos();
        const float avail = width - label_width;
        ImDrawList* draw_list = ImGui::GetWindowDrawList();
        draw_list->AddRectFilled(pos, ImVec2(pos.x + avail, pos.y + bar_height), IM_COL32(38, 44, 54, 180), 3.0f);
        const float center_x = pos.x + avail * 0.5f;
        draw_list->AddLine(ImVec2(center_x, pos.y), ImVec2(center_x, pos.y + bar_height), IM_COL32(220, 224, 232, 90), 1.0f);
        const float extent = static_cast<float>(std::abs(values[static_cast<std::size_t>(i)]) / max_abs) * (avail * 0.5f);
        const ImU32 col = ImGui::GetColorU32(component_color(values[static_cast<std::size_t>(i)]));
        if (values[static_cast<std::size_t>(i)] >= 0.0) {
            draw_list->AddRectFilled(ImVec2(center_x, pos.y), ImVec2(center_x + extent, pos.y + bar_height), col, 3.0f);
        } else {
            draw_list->AddRectFilled(ImVec2(center_x - extent, pos.y), ImVec2(center_x, pos.y + bar_height), col, 3.0f);
        }
        ImGui::Dummy(ImVec2(avail, bar_height));
        ImGui::SameLine();
        ImGui::Text("%+.3f", values[static_cast<std::size_t>(i)]);
    }
}

void draw_legend(float min_h, float max_h) {
    const float width = ImGui::GetContentRegionAvail().x;
    const float height = 22.0f;
    const ImVec2 pos = ImGui::GetCursorScreenPos();
    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    constexpr int steps = 100;
    for (int i = 0; i < steps; ++i) {
        const float t0 = static_cast<float>(i) / static_cast<float>(steps);
        const float t1 = static_cast<float>(i + 1) / static_cast<float>(steps);
        draw_list->AddRectFilled(ImVec2(pos.x + width * t0, pos.y),
                                 ImVec2(pos.x + width * t1, pos.y + height),
                                 ImGui::GetColorU32(esp_color(t0)));
    }
    draw_list->AddRect(pos, ImVec2(pos.x + width, pos.y + height), IM_COL32(215, 220, 228, 120), 3.0f);
    ImGui::Dummy(ImVec2(width, height + 6.0f));
    ImGui::TextUnformatted("d- (nucleophilic)");
    ImGui::SameLine();
    ImGui::SetCursorPosX(std::max(ImGui::GetCursorPosX(), ImGui::GetWindowWidth() * 0.42f));
    ImGui::TextUnformatted("neutral");
    ImGui::SameLine();
    ImGui::SetCursorPosX(std::max(ImGui::GetCursorPosX(), ImGui::GetWindowWidth() - 170.0f));
    ImGui::TextUnformatted("d+ (electrophilic)");
    ImGui::Text("%.1f kcal/mol                        %.1f kcal/mol",
                static_cast<double>(min_h) * kHartreeToKcalMol,
                static_cast<double>(max_h) * kHartreeToKcalMol);
}

}  // namespace

void draw_esp_controls(AppState& state, const sbox::backend::JobResult& result) {
    if (!result.converged()) {
        return;
    }

    if (!ImGui::Begin("Electrostatic Properties")) {
        ImGui::End();
        return;
    }

    if (ImGui::Button("<- Back to Dashboard")) {
        state.property_view = PropertyView::Dashboard;
    }
    ImGui::Separator();

    const Eigen::Vector3d dipole = result.dipole_moment;
    const double dipole_mag = dipole.norm();

    ImGui::TextUnformatted("Dipole Moment");
    ImGui::Text("|mu| = %.4f Debye", dipole_mag);
    ImGui::Text("  mux = %.4f D", dipole.x());
    ImGui::Text("  muy = %.4f D", dipole.y());
    ImGui::Text("  muz = %.4f D", dipole.z());
    ImGui::Checkbox("Show dipole arrow", &state.show_dipole);

    if (dipole_mag > 0.1) {
        ImGui::Text("Direction: primarily along %s axis", dominant_axis_label(dipole));
        draw_component_bars(dipole);
    }

    ImGui::Separator();
    ImGui::TextUnformatted("ESP Surface Controls");
    const bool has_esp = result.has_density_cube && result.has_esp_cube;
    if (!has_esp) {
        ImGui::TextDisabled("ESP surface unavailable for this result.");
    } else {
        ImGui::Checkbox("Show ESP Surface", &state.show_esp_surface);
        ImGui::SliderFloat("Density Iso Value",
                           &state.esp_density_iso,
                           0.0001f,
                           0.1f,
                           "%.4f",
                           ImGuiSliderFlags_Logarithmic);
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Controls the shape of the molecular surface. Lower values = larger surface.");
        }
        ImGui::Checkbox("Auto Color Range", &state.esp_auto_range);
        if (!state.esp_auto_range) {
            ImGui::SliderFloat("ESP Min (Hartree)", &state.esp_color_min, -0.2f, 0.0f, "%.4f");
            ImGui::SliderFloat("ESP Max (Hartree)", &state.esp_color_max, 0.0f, 0.2f, "%.4f");
            if (state.esp_color_min > state.esp_color_max) {
                std::swap(state.esp_color_min, state.esp_color_max);
            }
        }
        ImGui::Text("Color range: %.1f to %.1f kcal/mol",
                    static_cast<double>(state.esp_color_min) * kHartreeToKcalMol,
                    static_cast<double>(state.esp_color_max) * kHartreeToKcalMol);
    }

    ImGui::Separator();
    ImGui::TextUnformatted("ESP Colormap");
    draw_legend(state.esp_color_min, state.esp_color_max);

    ImGui::Separator();
    ImGui::TextUnformatted("Quick Info");
    ImGui::Text("Polarity: %s", polarity_label(dipole_mag));

    if (!result.mulliken_charges.empty()) {
        auto max_it = std::max_element(result.mulliken_charges.begin(), result.mulliken_charges.end());
        auto min_it = std::min_element(result.mulliken_charges.begin(), result.mulliken_charges.end());
        const int max_idx = static_cast<int>(std::distance(result.mulliken_charges.begin(), max_it));
        const int min_idx = static_cast<int>(std::distance(result.mulliken_charges.begin(), min_it));
        auto atom_label = [&](int atom_index) {
            if (result.has_optimized_geometry && atom_index >= 0 && atom_index < result.optimized_geometry.num_atoms()) {
                return std::string(sbox::elements::get_element(result.optimized_geometry.atom(atom_index).Z).symbol) +
                       std::to_string(atom_index + 1);
            }
            if (result.has_mo_data && atom_index >= 0 &&
                atom_index < static_cast<int>(result.mo_data.atomic_numbers.size())) {
                return std::string(sbox::elements::get_element(result.mo_data.atomic_numbers[static_cast<std::size_t>(atom_index)]).symbol) +
                       std::to_string(atom_index + 1);
            }
            return std::string("Atom ") + std::to_string(atom_index + 1);
        };
        if (max_idx >= 0) {
            ImGui::Text("Most electrophilic atom: %s (charge %+.2f)",
                        atom_label(max_idx).c_str(),
                        *max_it);
        }
        if (min_idx >= 0) {
            ImGui::Text("Most nucleophilic atom: %s (charge %+.2f)",
                        atom_label(min_idx).c_str(),
                        *min_it);
        }
    }

    ImGui::End();
}

}  // namespace sbox::ui
