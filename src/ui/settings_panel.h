#pragma once

#include "backend/python_env.h"
#include "core/settings.h"

namespace sbox::ui {

void draw_settings_panel(sbox::Settings& settings,
                         sbox::backend::PythonEnvironment& python_env,
                         bool& settings_changed);

bool consume_settings_panel_close_request();

}  // namespace sbox::ui
