#include "ui/plot_utils.h"

#include <implot.h>

namespace sbox::ui {

void setup_dark_plot_style() {
    ImPlotStyle& style = ImPlot::GetStyle();
    style.PlotBorderSize = 1.0f;

    ImVec4* colors = style.Colors;
    colors[ImPlotCol_PlotBg] = ImVec4(0.07f, 0.08f, 0.11f, 0.85f);
    colors[ImPlotCol_FrameBg] = ImVec4(0.0f, 0.0f, 0.0f, 0.0f);
    colors[ImPlotCol_AxisText] = ImVec4(0.82f, 0.84f, 0.88f, 1.0f);
    colors[ImPlotCol_AxisGrid] = ImVec4(0.20f, 0.24f, 0.30f, 0.50f);
    colors[ImPlotCol_AxisTick] = ImVec4(0.72f, 0.74f, 0.78f, 1.0f);
    colors[ImPlotCol_AxisBg] = ImVec4(0.0f, 0.0f, 0.0f, 0.0f);
    colors[ImPlotCol_LegendBg] = ImVec4(0.07f, 0.08f, 0.11f, 0.82f);
    colors[ImPlotCol_LegendBorder] = ImVec4(0.18f, 0.22f, 0.28f, 0.80f);
}

void plot_line_styled(const char* label,
                      const float* xs,
                      const float* ys,
                      int count,
                      ImVec4 color,
                      float thickness) {
    ImPlot::PlotLine(label, xs, ys, count, {ImPlotProp_LineColor, color, ImPlotProp_LineWeight, thickness});
}

void plot_stems_styled(const char* label,
                       const float* xs,
                       const float* ys,
                       int count,
                       ImVec4 color,
                       float thickness) {
    ImPlot::PlotStems(label, xs, ys, count, 0.0, {ImPlotProp_LineColor, color, ImPlotProp_LineWeight, thickness});
}

void plot_shaded_styled(const char* label,
                        const float* xs,
                        const float* ys,
                        int count,
                        ImVec4 color,
                        float alpha) {
    ImVec4 shaded = color;
    shaded.w = alpha;
    ImPlot::PlotShaded(label, xs, ys, count, 0.0, {ImPlotProp_FillColor, shaded, ImPlotProp_FillAlpha, alpha});
}

ImVec4 orbital_color(int mo_index, int homo_index) {
    if (homo_index >= 0) {
        if (mo_index == homo_index) {
            return ImVec4(0.22f, 0.62f, 0.92f, 1.0f);
        }
        if (mo_index == homo_index + 1) {
            return ImVec4(0.95f, 0.58f, 0.18f, 1.0f);
        }
        if (mo_index <= homo_index) {
            return ImVec4(0.18f, 0.46f, 0.82f, 1.0f);
        }
    }
    return ImVec4(0.56f, 0.58f, 0.64f, 1.0f);
}

}  // namespace sbox::ui
