#pragma once

#include "core/update_checker.h"
#include "ui/app_state.h"

namespace sbox::ui {

void draw_status_bar(const AppState& state, bool show_fps, const sbox::UpdateInfo* pending_update = nullptr, bool* update_clicked = nullptr);

}  // namespace sbox::ui
