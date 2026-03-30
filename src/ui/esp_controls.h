#pragma once

#include "backend/job_types.h"
#include "ui/app_state.h"

namespace sbox::ui {

void draw_esp_controls(AppState& state, const sbox::backend::JobResult& result);

}  // namespace sbox::ui
