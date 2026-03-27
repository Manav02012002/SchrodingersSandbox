#include "ui/mo_diagram.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace sbox::ui {
namespace {

constexpr float kHartreeToEv = 27.2114f;

struct MOLine {
    int index = 0;
    float energy_h = 0.0f;
    float energy_ev = 0.0f;
    float occupation = 0.0f;
    bool occupied = false;
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

void draw_mo_diagram(AppState& state, const sbox::basis::MOData& mo_data) {
    if (!ImGui::Begin("MO Energy Levels")) {
        ImGui::End();
        return;
    }

    if (mo_data.energies.size() == 0) {
        ImGui::TextUnformatted("No molecular orbitals available.");
        ImGui::End();
        return;
    }

    std::vector<MOLine> levels;
    levels.reserve(static_cast<std::size_t>(mo_data.energies.size()));

    int homo = -1;
    int lumo = -1;
    for (int i = 0; i < mo_data.energies.size(); ++i) {
        const float occ = i < mo_data.occupations.size() ? static_cast<float>(mo_data.occupations(i)) : 0.0f;
        const bool occupied = occ > 0.5f;
        if (occupied) {
            homo = i;
        } else if (lumo < 0) {
            lumo = i;
        }

        MOLine level;
        level.index = i;
        level.energy_h = static_cast<float>(mo_data.energies(i));
        level.energy_ev = level.energy_h * kHartreeToEv;
        level.occupation = occ;
        level.occupied = occupied;
        level.selected = (state.selected_mo == i);
        levels.push_back(level);
    }
    if (lumo < 0 && homo + 1 < static_cast<int>(levels.size())) {
        lumo = homo + 1;
    }

    state.num_mo = static_cast<int>(levels.size());
    state.homo_index = homo;
    if (state.selected_mo < 0 && homo >= 0) {
        state.selected_mo = homo;
    }
    state.selected_mo = std::clamp(state.selected_mo, 0, static_cast<int>(levels.size()) - 1);

    float gap_ev = 0.0f;
    bool has_gap = homo >= 0 && lumo >= 0 && lumo < static_cast<int>(levels.size());
    if (has_gap) {
        gap_ev = levels[static_cast<std::size_t>(lumo)].energy_ev - levels[static_cast<std::size_t>(homo)].energy_ev;
    }

    ImGui::Text("Orbitals: %d", static_cast<int>(levels.size()));
    if (has_gap) {
        ImGui::SameLine();
        ImGui::Text("HOMO-LUMO gap: %.3f eV", gap_ev);
    }
    ImGui::Separator();

    static int visible_count = 15;
    visible_count = std::clamp(visible_count, 8, std::max(8, static_cast<int>(levels.size())));
    int center_index = homo >= 0 ? homo : state.selected_mo;
    if (lumo >= 0 && homo >= 0) {
        center_index = (homo + lumo) / 2;
    }
    if (state.selected_mo >= 0) {
        center_index = state.selected_mo;
    }

    int start = 0;
    int end = static_cast<int>(levels.size());
    if (static_cast<int>(levels.size()) > 30) {
        start = std::max(0, center_index - visible_count / 2);
        end = std::min(static_cast<int>(levels.size()), start + visible_count);
        start = std::max(0, end - visible_count);
    }

    const ImVec2 canvas_origin = ImGui::GetCursorScreenPos();
    ImVec2 canvas_size = ImGui::GetContentRegionAvail();
    canvas_size.x = std::max(canvas_size.x, 180.0f);
    canvas_size.y = std::max(canvas_size.y, 180.0f);
    const ImVec2 canvas_max(canvas_origin.x + canvas_size.x, canvas_origin.y + canvas_size.y);

    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    draw_list->AddRectFilled(canvas_origin, canvas_max, ImGui::GetColorU32(ImVec4(0.07f, 0.08f, 0.11f, 0.85f)), 4.0f);
    draw_list->AddRect(canvas_origin, canvas_max, ImGui::GetColorU32(ImVec4(0.20f, 0.24f, 0.30f, 1.0f)), 4.0f);

    const float left_pad = 64.0f;
    const float right_pad = 18.0f;
    const float top_pad = 18.0f;
    const float bottom_pad = 18.0f;
    const float x0 = canvas_origin.x + left_pad;
    const float x1 = canvas_max.x - right_pad;
    const float y_top = canvas_origin.y + top_pad;
    const float y_bottom = canvas_max.y - bottom_pad;

    float min_ev = levels[static_cast<std::size_t>(start)].energy_ev;
    float max_ev = min_ev;
    for (int i = start; i < end; ++i) {
        min_ev = std::min(min_ev, levels[static_cast<std::size_t>(i)].energy_ev);
        max_ev = std::max(max_ev, levels[static_cast<std::size_t>(i)].energy_ev);
    }
    if (std::abs(max_ev - min_ev) < 1e-5f) {
        max_ev = min_ev + 1.0f;
    }
    const float pad = std::max(0.5f, 0.08f * (max_ev - min_ev));
    min_ev -= pad;
    max_ev += pad;

    for (int tick = 0; tick <= 4; ++tick) {
        const float t = static_cast<float>(tick) / 4.0f;
        const float y = y_bottom - t * (y_bottom - y_top);
        const float value_ev = min_ev + t * (max_ev - min_ev);
        draw_list->AddLine(ImVec2(x0 - 6.0f, y),
                           ImVec2(x1, y),
                           ImGui::GetColorU32(ImVec4(0.18f, 0.22f, 0.28f, tick == 0 || tick == 4 ? 0.8f : 0.45f)),
                           1.0f);
        char label[32];
        std::snprintf(label, sizeof(label), "%.1f", value_ev);
        draw_list->AddText(ImVec2(canvas_origin.x + 8.0f, y - 8.0f),
                           ImGui::GetColorU32(ImVec4(0.72f, 0.76f, 0.84f, 0.9f)),
                           label);
    }

    draw_list->AddText(ImVec2(canvas_origin.x + 8.0f, canvas_origin.y + 2.0f),
                       ImGui::GetColorU32(ImVec4(0.80f, 0.84f, 0.92f, 0.9f)),
                       "eV");

    int hovered_index = -1;
    float hovered_distance = 1e9f;

    for (int i = start; i < end; ++i) {
        const MOLine& level = levels[static_cast<std::size_t>(i)];
        const float t = (level.energy_ev - min_ev) / (max_ev - min_ev);
        const float y = y_bottom - t * (y_bottom - y_top);

        ImU32 level_color = level.selected
                                ? ImGui::GetColorU32(ImVec4(0.15f, 0.55f, 0.65f, 1.0f))
                                : (level.occupied ? ImGui::GetColorU32(ImVec4(0.24f, 0.72f, 0.82f, 0.95f))
                                                  : ImGui::GetColorU32(ImVec4(0.62f, 0.66f, 0.72f, 0.85f)));
        const float thickness = level.selected ? 3.0f : 1.6f;

        draw_list->AddLine(ImVec2(x0, y), ImVec2(x1 - 64.0f, y), level_color, thickness);

        std::string index_label = std::to_string(level.index);
        draw_list->AddText(ImVec2(x1 - 56.0f, y - 8.0f), level_color, index_label.c_str());

        if (level.index == homo) {
            draw_list->AddText(ImVec2(x0 + 8.0f, y - 18.0f), level_color, "HOMO");
        } else if (level.index == lumo) {
            draw_list->AddText(ImVec2(x0 + 8.0f, y - 18.0f), level_color, "LUMO");
        }

        const float arrow_base_x = x0 + 18.0f;
        const float arrow_size = 9.0f;
        if (level.occupation >= 1.0f) {
            draw_arrow(draw_list, ImVec2(arrow_base_x, y - 5.0f), arrow_size, level_color, true);
        }
        if (level.occupation >= 2.0f) {
            draw_arrow(draw_list,
                       ImVec2(arrow_base_x + 6.0f, y + 5.0f),
                       arrow_size,
                       ImGui::GetColorU32(ImVec4(0.95f, 0.95f, 0.95f, 0.9f)),
                       false);
        }

        if (ImGui::IsMouseHoveringRect(ImVec2(x0, y - 6.0f), ImVec2(x1, y + 6.0f))) {
            const float dist = std::abs(ImGui::GetIO().MousePos.y - y);
            if (dist < hovered_distance) {
                hovered_distance = dist;
                hovered_index = level.index;
            }
        }
    }

    ImGui::InvisibleButton("##mo_canvas", canvas_size);
    const bool hovered = ImGui::IsItemHovered();
    if (hovered) {
        const float wheel = ImGui::GetIO().MouseWheel;
        if (std::abs(wheel) > 0.0f && static_cast<int>(levels.size()) > 30) {
            visible_count = std::clamp(visible_count - (wheel > 0.0f ? 1 : -1) * 2,
                                       8,
                                       static_cast<int>(levels.size()));
        }
        if (hovered_index >= 0) {
            if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                state.selected_mo = hovered_index;
            }
            const MOLine& level = levels[static_cast<std::size_t>(hovered_index)];
            if (ImGui::BeginTooltip()) {
                ImGui::Text("MO %d", level.index);
                ImGui::Text("Energy: %.6f Ha", level.energy_h);
                ImGui::Text("Energy: %.3f eV", level.energy_ev);
                ImGui::EndTooltip();
            }
        }
    }

    ImGui::End();
}

}  // namespace sbox::ui
