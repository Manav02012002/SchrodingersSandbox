#include "ui/status_bar.h"

#include "core/elements.h"
#include "ui/ui_utils.h"

#include <imgui.h>
#include <imgui_internal.h>

#include <algorithm>
#include <cstddef>

namespace sbox::ui {

void draw_status_bar(const AppState& state, bool show_fps, const sbox::UpdateInfo* pending_update, bool* update_clicked) {
    ImGuiViewport* viewport = ImGui::GetMainViewport();
    constexpr float kStatusBarHeight = 28.0f;
    const ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoDocking;

    if (!ImGui::BeginViewportSideBar("##StatusBar", viewport, ImGuiDir_Down, kStatusBarHeight, flags)) {
        ImGui::End();
        return;
    }

    const auto& element = sbox::elements::get_element(state.selected_Z);
    const int clamped_l = std::max(state.current_l, 0);
    const int l_index = std::clamp(clamped_l, 0, static_cast<int>(kLLabels.size() - 1));
    const int m_value = std::clamp(state.selected_m, -clamped_l, clamped_l);
    const float fps = ImGui::GetIO().Framerate;

    ImGui::Text("Element: %s (%s)", element.name, element.symbol);
    ImGui::SameLine();
    ImGui::TextUnformatted(" | ");
    ImGui::SameLine();
    ImGui::Text("Orbital: %d%s, m=%d", state.current_n, kLLabels[static_cast<std::size_t>(l_index)], m_value);
    ImGui::SameLine();
    ImGui::TextUnformatted(" | ");
    ImGui::SameLine();
    ImGui::Text("Zeff: %.3f", state.current_Zeff);
    if (show_fps) {
        ImGui::SameLine();
        ImGui::TextUnformatted(" | ");
        ImGui::SameLine();
        ImGui::Text("FPS: %.1f", fps);
    }
    if (state.view_mode == ViewMode::MolecularOrbital && state.molecule_loaded && state.lod_atoms_rendered > 0) {
        ImGui::SameLine();
        ImGui::TextUnformatted(" | ");
        ImGui::SameLine();
        ImGui::Text("Rendered: %d/%d atoms", state.lod_atoms_rendered, state.lod_atoms_rendered + state.lod_atoms_culled);
    }
    if (pending_update != nullptr && pending_update->update_available) {
        ImGui::SameLine();
        ImGui::TextUnformatted(" | ");
        ImGui::SameLine();
        ImGui::TextColored(ImVec4(0.2f, 0.8f, 0.3f, 1.0f), "Update available: %s", pending_update->latest_version.c_str());
        if (update_clicked != nullptr && ImGui::IsItemClicked()) {
            *update_clicked = true;
        }
    }

    ImGui::End();
}

}  // namespace sbox::ui
