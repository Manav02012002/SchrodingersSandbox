#pragma once

#include <imgui.h>

namespace sbox::ui {

void setup_dark_plot_style();

void plot_line_styled(const char* label,
                      const float* xs,
                      const float* ys,
                      int count,
                      ImVec4 color,
                      float thickness = 1.5f);

void plot_stems_styled(const char* label,
                       const float* xs,
                       const float* ys,
                       int count,
                       ImVec4 color,
                       float thickness = 1.0f);

void plot_shaded_styled(const char* label,
                        const float* xs,
                        const float* ys,
                        int count,
                        ImVec4 color,
                        float alpha = 0.3f);

ImVec4 orbital_color(int mo_index, int homo_index);

}  // namespace sbox::ui
