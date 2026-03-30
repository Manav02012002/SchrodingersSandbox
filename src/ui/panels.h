#pragma once

#include "ui/app_state.h"
#include "ui/energy_diagram.h"
#include "ui/orbital_browser.h"
#include "ui/periodic_table.h"
#include "ui/properties_panel.h"
#include "ui/status_bar.h"

#include <imgui.h>

namespace sbox::ui {

struct ViewportPanelState {
    ImVec2 pos{0.0f, 0.0f};
    ImVec2 size{0.0f, 0.0f};
    bool hovered = false;
};

ViewportPanelState draw_viewport(AppState& state, unsigned int texture_id);
const char* mo_label_for_index(const AppState& state, int mo_index, int homo_index);

}  // namespace sbox::ui
