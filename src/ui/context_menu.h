#pragma once

#include "editor/command.h"
#include "editor/picking.h"

#include <imgui.h>

namespace sbox::ui {

struct ContextMenuState {
    bool show = false;
    ImVec2 position{0.0f, 0.0f};
    int clicked_atom = -1;
    int clicked_bond = -1;
    bool center_view_requested = false;
    bool fit_view_requested = false;
};

void draw_context_menu(
    ContextMenuState& ctx,
    sbox::chem::MolecularSystem& mol,
    sbox::editor::Selection& selection,
    sbox::editor::CommandStack& commands);

}  // namespace sbox::ui
