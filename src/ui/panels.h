#pragma once

#include "ui/app_state.h"
#include "ui/orbital_browser.h"
#include "ui/periodic_table.h"
#include "ui/properties_panel.h"

#include <imgui.h>

namespace sbox::ui {

struct ViewportPanelState {
    ImVec2 size{0.0f, 0.0f};
    bool hovered = false;
};

ViewportPanelState draw_viewport(AppState& state, unsigned int texture_id);

}  // namespace sbox::ui
