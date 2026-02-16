#include "ui/periodic_table.h"

#include "core/elements.h"

#include <imgui.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cfloat>
#include <string>

namespace sbox::ui {
namespace {

constexpr std::array<const char*, 7> kLLabels = {"s", "p", "d", "f", "g", "h", "i"};

std::string to_lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

std::string superscript_number(int value) {
    static const std::array<const char*, 10> kDigits = {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};

    if (value == 0) {
        return kDigits[0];
    }

    std::string out;
    if (value < 0) {
        out += "⁻";
        value = -value;
    }

    std::string digits = std::to_string(value);
    for (char ch : digits) {
        out += kDigits[static_cast<std::size_t>(ch - '0')];
    }
    return out;
}

std::string config_to_string(const sbox::slater::ElectronConfig& config) {
    std::string out;
    for (std::size_t i = 0; i < config.size(); ++i) {
        const auto& subshell = config[i];
        if (i > 0) {
            out += " ";
        }

        out += std::to_string(subshell.n);
        const int l = std::clamp(subshell.l, 0, static_cast<int>(kLLabels.size() - 1));
        out += kLLabels[static_cast<std::size_t>(l)];
        out += superscript_number(subshell.electrons);
    }
    return out;
}

bool matches_filter(const sbox::elements::Element& element, const std::string& filter) {
    if (filter.empty()) {
        return true;
    }

    const std::string name = to_lower(element.name);
    const std::string symbol = to_lower(element.symbol);
    const std::string z_text = std::to_string(element.Z);

    return name.find(filter) != std::string::npos ||
           symbol.find(filter) != std::string::npos ||
           z_text.find(filter) != std::string::npos;
}

ImVec4 category_color(const char* category) {
    const std::string value = category != nullptr ? category : "";

    if (value == "alkali") {
        return ImVec4(0.97f, 0.44f, 0.44f, 0.7f);
    }
    if (value == "alkaline") {
        return ImVec4(0.98f, 0.57f, 0.37f, 0.7f);
    }
    if (value == "transition") {
        return ImVec4(0.38f, 0.65f, 0.98f, 0.7f);
    }
    if (value == "post-transition") {
        return ImVec4(0.20f, 0.82f, 0.60f, 0.7f);
    }
    if (value == "metalloid") {
        return ImVec4(0.98f, 0.75f, 0.15f, 0.7f);
    }
    if (value == "nonmetal") {
        return ImVec4(0.18f, 0.83f, 0.75f, 0.7f);
    }
    if (value == "halogen") {
        return ImVec4(0.91f, 0.47f, 0.95f, 0.7f);
    }
    if (value == "noble") {
        return ImVec4(0.65f, 0.55f, 0.98f, 0.7f);
    }
    if (value == "lanthanide") {
        return ImVec4(0.96f, 0.44f, 0.70f, 0.7f);
    }
    if (value == "actinide") {
        return ImVec4(0.75f, 0.52f, 0.98f, 0.7f);
    }

    return ImVec4(0.30f, 0.30f, 0.30f, 0.7f);
}

ImVec4 brighten(ImVec4 color, float amount) {
    color.x = std::min(1.0f, color.x * amount);
    color.y = std::min(1.0f, color.y * amount);
    color.z = std::min(1.0f, color.z * amount);
    return color;
}

}  // namespace

void draw_periodic_table(AppState& state) {
    ImGui::SetNextWindowSize(ImVec2(400, 300), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Periodic Table")) {
        ImGui::End();
        return;
    }

    static char search[64] = "";
    ImGui::InputText("Search", search, sizeof(search));
    const std::string filter = to_lower(search);

    const float spacing = 4.0f;
    const int total_cols = 18;
    const int total_rows = 10;

    if (ImGui::BeginChild("##periodic_grid", ImVec2(0.0f, 0.0f), false, ImGuiWindowFlags_HorizontalScrollbar)) {
        const ImVec2 avail = ImGui::GetContentRegionAvail();
        float cell_size = (avail.x - spacing * static_cast<float>(total_cols - 1)) / static_cast<float>(total_cols);
        if (avail.y > 0.0f) {
            const float max_cell_by_height = (avail.y - spacing * static_cast<float>(total_rows - 1)) / static_cast<float>(total_rows);
            cell_size = std::min(cell_size, max_cell_by_height);
        }
        cell_size = std::clamp(cell_size, 24.0f, 46.0f);

        const float grid_width = static_cast<float>(total_cols) * cell_size + spacing * static_cast<float>(total_cols - 1);
        const float grid_height = static_cast<float>(total_rows) * cell_size + spacing * static_cast<float>(total_rows - 1);

        const ImVec2 origin = ImGui::GetCursorScreenPos();
        ImDrawList* draw_list = ImGui::GetWindowDrawList();
        ImFont* font = ImGui::GetFont();
        const float base_font_size = ImGui::GetFontSize();
        const float small_font_size = base_font_size * 0.7f;
        const float symbol_font_size = base_font_size * 1.2f;

        const ImU32 label_color = ImGui::GetColorU32(ImVec4(0.70f, 0.70f, 0.70f, 0.85f));
        const float row8_y = origin.y + 8.0f * (cell_size + spacing);
        const float row9_y = origin.y + 9.0f * (cell_size + spacing);
        draw_list->AddText(font, small_font_size, ImVec2(origin.x + 2.0f, row8_y + 0.2f * cell_size), label_color, "La–Lu");
        draw_list->AddText(font, small_font_size, ImVec2(origin.x + 2.0f, row9_y + 0.2f * cell_size), label_color, "Ac–Lr");

        for (const auto& pos : sbox::elements::PT_LAYOUT) {
            const auto& element = sbox::elements::get_element(pos.Z);
            const float x = origin.x + static_cast<float>(pos.col) * (cell_size + spacing);
            const float y = origin.y + static_cast<float>(pos.row) * (cell_size + spacing);
            const ImVec2 min(x, y);
            const ImVec2 max(x + cell_size, y + cell_size);

            ImGui::SetCursorScreenPos(min);
            ImGui::PushID(element.Z);
            ImGui::InvisibleButton("##element", ImVec2(cell_size, cell_size));
            const bool hovered = ImGui::IsItemHovered();
            const bool clicked = ImGui::IsItemClicked(ImGuiMouseButton_Left);
            ImGui::PopID();

            if (clicked) {
                state.selected_Z = element.Z;
                state.selected_orbital_index = -1;
                state.selected_m = 0;
                state.needs_update = true;
            }

            const bool selected = (state.selected_Z == element.Z);
            const bool match = matches_filter(element, filter);

            ImVec4 color = category_color(element.category);
            if (!match) {
                color.w = std::min(color.w, 0.2f);
            }
            if (hovered) {
                color = brighten(color, 1.15f);
            }

            draw_list->AddRectFilled(min, max, ImGui::GetColorU32(color), 4.0f);
            const ImU32 border_color = selected
                                        ? ImGui::GetColorU32(ImVec4(1.0f, 1.0f, 1.0f, 1.0f))
                                        : ImGui::GetColorU32(ImVec4(0.18f, 0.20f, 0.23f, 0.9f));
            const float border_thickness = selected ? 2.0f : 1.0f;
            draw_list->AddRect(min, max, border_color, 4.0f, 0, border_thickness);

            const ImU32 text_color = ImGui::GetColorU32(ImVec4(0.96f, 0.97f, 1.0f, match ? 0.95f : 0.50f));
            const std::string z_text = std::to_string(element.Z);
            draw_list->AddText(font,
                               small_font_size,
                               ImVec2(min.x + 3.0f, min.y + 2.0f),
                               text_color,
                               z_text.c_str());

            const ImVec2 symbol_size = font->CalcTextSizeA(symbol_font_size,
                                                           FLT_MAX,
                                                           0.0f,
                                                           element.symbol,
                                                           nullptr,
                                                           nullptr);
            draw_list->AddText(font,
                               symbol_font_size,
                               ImVec2(min.x + (cell_size - symbol_size.x) * 0.5f,
                                      min.y + (cell_size - symbol_size.y) * 0.5f + 1.0f),
                               text_color,
                               element.symbol);

            if (hovered) {
                const std::string config = config_to_string(element.config);
                ImGui::BeginTooltip();
                ImGui::Text("%s (%s)", element.name, element.symbol);
                ImGui::Text("Z = %d", element.Z);
                ImGui::TextUnformatted(config.c_str());
                ImGui::TextUnformatted(element.category);
                ImGui::EndTooltip();
            }
        }

        ImGui::SetCursorScreenPos(ImVec2(origin.x, origin.y + grid_height));
        ImGui::Dummy(ImVec2(grid_width, 1.0f));
    }
    ImGui::EndChild();

    ImGui::End();
}

}  // namespace sbox::ui
