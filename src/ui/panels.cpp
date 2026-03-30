#include "ui/panels.h"

#include <cstdint>
#include <string>

namespace {

std::string g_mo_label;

}  // namespace

namespace sbox::ui {

const char* mo_label_for_index(const AppState& state, int mo_index, int homo_index) {
    (void)state;
    if (homo_index < 0) {
        g_mo_label = "MO";
        return g_mo_label.c_str();
    }
    if (mo_index == homo_index) {
        g_mo_label = "HOMO";
    } else if (mo_index == homo_index + 1) {
        g_mo_label = "LUMO";
    } else if (mo_index < homo_index) {
        g_mo_label = "HOMO-" + std::to_string(homo_index - mo_index);
    } else {
        g_mo_label = "LUMO+" + std::to_string(mo_index - homo_index - 1);
    }
    return g_mo_label.c_str();
}

ViewportPanelState draw_viewport(AppState& state, unsigned int texture_id) {
    ViewportPanelState viewport{};

    if (ImGui::Begin("3D Viewport")) {
        int view_mode = static_cast<int>(state.view_mode);
        const char* view_items[] = {"Atomic Orbital", "Molecular Orbital"};
        ImGui::SetNextItemWidth(220.0f);
        if (ImGui::Combo("View", &view_mode, view_items, 2)) {
            state.view_mode = static_cast<ViewMode>(view_mode);
        }

        if (state.view_mode == ViewMode::MolecularOrbital) {
            if (state.num_mo > 0) {
                if (state.selected_mo < 0) {
                    state.selected_mo = state.homo_index >= 0 ? state.homo_index : 0;
                }
                state.selected_mo = std::clamp(state.selected_mo, 0, state.num_mo - 1);
                ImGui::SetNextItemWidth(220.0f);
                ImGui::SliderInt("Molecular Orbital", &state.selected_mo, 0, state.num_mo - 1);
                ImGui::Text("%s", mo_label_for_index(state, state.selected_mo, state.homo_index));
            } else {
                ImGui::TextUnformatted("No molecular orbitals loaded.");
            }
            ImGui::SetNextItemWidth(220.0f);
            const char* mol_items[] = {"Ball-and-Stick", "Space-Filling", "Wireframe", "Stick"};
            ImGui::Combo("Molecule Style", &state.mol_render_mode, mol_items, 4);
        }

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
        ImGui::TextUnformatted("Overlays");
        ImGui::Checkbox("Partial Charges", &state.show_charges);
        ImGui::Checkbox("Bond Orders", &state.show_bond_orders);
        ImGui::Checkbox("Dipole Moment", &state.show_dipole);
        ImGui::Checkbox("ESP Surface", &state.show_esp_surface);
        ImGui::Separator();

        viewport.size = ImGui::GetContentRegionAvail();
        viewport.pos = ImGui::GetCursorScreenPos();

        if (viewport.size.x > 0.0f && viewport.size.y > 0.0f && texture_id != 0U) {
            ImGui::Image((ImTextureID)(intptr_t)texture_id,
                         viewport.size,
                         ImVec2(0.0f, 1.0f),
                         ImVec2(1.0f, 0.0f));
            viewport.hovered = ImGui::IsItemHovered();
        } else {
            viewport.hovered = false;
            ImGui::TextUnformatted("Viewport unavailable.");
        }
    }
    ImGui::End();

    return viewport;
}

}  // namespace sbox::ui
