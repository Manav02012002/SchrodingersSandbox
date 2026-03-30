#pragma once

#include "backend/backend_manager.h"
#include "ui/app_state.h"

namespace sbox::ui {

void draw_computation_panel(AppState& state, sbox::backend::BackendManager& backend);

}  // namespace sbox::ui
