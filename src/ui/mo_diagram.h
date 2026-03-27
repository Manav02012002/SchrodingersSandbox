#pragma once

#include "core/basis_set.h"
#include "ui/app_state.h"

namespace sbox::ui {

void draw_mo_diagram(AppState& state, const sbox::basis::MOData& mo_data);

}  // namespace sbox::ui
