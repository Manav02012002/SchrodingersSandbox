#include "ui/energy_diagram.h"

#include "core/elements.h"
#include "core/slater.h"
#include "ui/ui_utils.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace sbox::ui {
namespace {

struct LevelVisual {
    int n = 1;
    int l = 0;
    int electrons = 0;
    float log_energy = 0.0f;
    bool selected = false;
};

void draw_arrow(ImDrawList* draw_list, ImVec2 center, float size, ImU32 color, bool up) {
    const float half_w = size * 0.26f;
    const float body_h = size * 0.52f;
    const float tip_h = size * 0.28f;

    if (up) {
        draw_list->AddLine(ImVec2(center.x, center.y + body_h * 0.5f),
                           ImVec2(center.x, center.y - body_h * 0.5f),
                           color,
                           1.2f);
        draw_list->AddTriangleFilled(ImVec2(center.x, center.y - body_h * 0.5f - tip_h),
                                     ImVec2(center.x - half_w, center.y - body_h * 0.5f),
                                     ImVec2(center.x + half_w, center.y - body_h * 0.5f),
                                     color);
    } else {
        draw_list->AddLine(ImVec2(center.x, center.y - body_h * 0.5f),
                           ImVec2(center.x, center.y + body_h * 0.5f),
                           color,
                           1.2f);
        draw_list->AddTriangleFilled(ImVec2(center.x, center.y + body_h * 0.5f + tip_h),
                                     ImVec2(center.x - half_w, center.y + body_h * 0.5f),
                                     ImVec2(center.x + half_w, center.y + body_h * 0.5f),
                                     color);
    }
}

}  // namespace

void draw_energy_diagram(AppState& state) {
    if (!ImGui::Begin("Energy Levels")) {
        ImGui::End();
        return;
    }

    const auto& element = sbox::elements::get_element(state.selected_Z);
    const auto& config = element.config;
    if (config.empty()) {
        ImGui::TextUnformatted("No orbitals available.");
        ImGui::End();
        return;
    }

    int selected_index = state.selected_orbital_index;
    if (selected_index < 0 || selected_index >= static_cast<int>(config.size())) {
        selected_index = static_cast<int>(config.size()) - 1;
    }

    std::vector<LevelVisual> levels;
    levels.reserve(config.size());
    for (std::size_t i = 0; i < config.size(); ++i) {
        const auto& subshell = config[i];
        if (subshell.electrons <= 0) {
            continue;
        }
        const float zeff = static_cast<float>(sbox::slater::compute_zeff(element.Z, config, subshell.n, subshell.l));
        const float energy = -13.6f * std::pow(zeff / static_cast<float>(subshell.n), 2.0f);
        const float minus_e = std::max(1e-6f, -energy);

        LevelVisual level;
        level.n = subshell.n;
        level.l = subshell.l;
        level.electrons = subshell.electrons;
        level.log_energy = std::log(minus_e);
        level.selected = static_cast<int>(i) == selected_index;
        levels.push_back(level);
    }

    if (levels.empty()) {
        ImGui::TextUnformatted("No occupied levels.");
        ImGui::End();
        return;
    }

    float min_log_e = levels.front().log_energy;
    float max_log_e = levels.front().log_energy;
    for (const LevelVisual& level : levels) {
        min_log_e = std::min(min_log_e, level.log_energy);
        max_log_e = std::max(max_log_e, level.log_energy);
    }
    if (std::abs(max_log_e - min_log_e) < 1e-6f) {
        max_log_e = min_log_e + 1.0f;
    }

    ImGui::Text("Energy Levels (%s)", element.symbol);
    ImGui::Separator();

    const ImVec2 canvas_origin = ImGui::GetCursorScreenPos();
    ImVec2 canvas_size = ImGui::GetContentRegionAvail();
    canvas_size.x = std::max(canvas_size.x, 120.0f);
    canvas_size.y = std::max(canvas_size.y, 120.0f);

    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    const ImVec2 canvas_max(canvas_origin.x + canvas_size.x, canvas_origin.y + canvas_size.y);

    draw_list->AddRectFilled(canvas_origin, canvas_max, ImGui::GetColorU32(ImVec4(0.07f, 0.08f, 0.11f, 0.85f)), 4.0f);
    draw_list->AddRect(canvas_origin, canvas_max, ImGui::GetColorU32(ImVec4(0.20f, 0.24f, 0.30f, 1.0f)), 4.0f);

    const float left_pad = 56.0f;
    const float right_pad = 14.0f;
    const float top_pad = 16.0f;
    const float bottom_pad = 14.0f;
    const float x0 = canvas_origin.x + left_pad;
    const float x1 = canvas_max.x - right_pad;
    const float y_top = canvas_origin.y + top_pad;
    const float y_bottom = canvas_max.y - bottom_pad;

    for (const LevelVisual& level : levels) {
        const float t = (level.log_energy - min_log_e) / (max_log_e - min_log_e);
        const float y = y_top + t * (y_bottom - y_top);

        const ImU32 level_color = level.selected
                                      ? ImGui::GetColorU32(ImVec4(0.15f, 0.55f, 0.65f, 1.0f))
                                      : ImGui::GetColorU32(ImVec4(0.74f, 0.78f, 0.88f, 0.95f));
        const float thickness = level.selected ? 2.6f : 1.5f;
        draw_list->AddLine(ImVec2(x0, y), ImVec2(x1, y), level_color, thickness);

        std::string label = std::to_string(level.n) + kLLabels[static_cast<std::size_t>(std::clamp(level.l, 0, 6))];
        draw_list->AddText(ImVec2(canvas_origin.x + 10.0f, y - 8.0f), level_color, label.c_str());

        const int orbital_count = 2 * level.l + 1;
        const int up_count = std::min(level.electrons, orbital_count);
        const int down_count = std::max(0, level.electrons - orbital_count);
        const float arrow_spacing = std::min(14.0f, (x1 - x0) / std::max(1, orbital_count));
        const float occupancy_start = x0 + 12.0f;
        const float arrow_size = 9.0f;

        for (int orb = 0; orb < orbital_count; ++orb) {
            const float cx = occupancy_start + static_cast<float>(orb) * arrow_spacing;
            if (orb < up_count) {
                draw_arrow(draw_list, ImVec2(cx, y - 5.0f), arrow_size, level_color, true);
            }
            if (orb < down_count) {
                draw_arrow(draw_list,
                           ImVec2(cx + 4.0f, y + 5.0f),
                           arrow_size,
                           ImGui::GetColorU32(ImVec4(0.95f, 0.95f, 0.95f, level.selected ? 1.0f : 0.9f)),
                           false);
            }
        }
    }

    ImGui::InvisibleButton("##energy_canvas", canvas_size);
    ImGui::End();
}

}  // namespace sbox::ui
