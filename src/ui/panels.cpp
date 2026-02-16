#include "ui/panels.h"

#include <cstdint>

namespace sbox::ui {

ViewportPanelState draw_viewport(AppState& state, unsigned int texture_id) {
    ViewportPanelState viewport{};

    ImGui::SetNextWindowSize(ImVec2(400, 300), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("3D Viewport")) {
        const char* mode_items[] = {"Volume", "Isosurface", "Phase Isosurface"};
        ImGui::SetNextItemWidth(220.0f);
        ImGui::Combo("Render Mode", &state.render_mode, mode_items, 3);

        if (state.render_mode > 0) {
            ImGui::SetNextItemWidth(220.0f);
            ImGui::SliderFloat("Iso Value",
                               &state.iso_value,
                               0.0001f,
                               0.5f,
                               "%.5f",
                               ImGuiSliderFlags_Logarithmic);
        }

        if (state.render_mode == 0) {
            ImGui::SetNextItemWidth(220.0f);
            ImGui::SliderFloat("Gamma", &state.gamma, 0.1f, 1.0f, "%.2f");
        }

        ImGui::Separator();

        viewport.size = ImGui::GetContentRegionAvail();
        viewport.hovered = ImGui::IsWindowHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);

        if (viewport.size.x > 0.0f && viewport.size.y > 0.0f && texture_id != 0U) {
            ImGui::Image((ImTextureID)(intptr_t)texture_id,
                         viewport.size,
                         ImVec2(0.0f, 1.0f),
                         ImVec2(1.0f, 0.0f));
        } else {
            ImGui::TextUnformatted("Viewport unavailable.");
        }
    }
    ImGui::End();

    return viewport;
}

}  // namespace sbox::ui
