#pragma once

#include "backend/python_env.h"

namespace sbox::ui {

void set_about_dialog_context(const sbox::backend::PythonEnvironment* python_env);
void draw_about_dialog(bool& show);
void draw_shortcuts_dialog(bool& show);

}  // namespace sbox::ui
