#include "ui/ir_spectrum_panel.h"

#include "ui/plot_utils.h"

#include <imgui.h>
#include <implot.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace sbox::ui {

namespace {

constexpr int kGridPoints = 1000;
constexpr double kCmToMicron = 10000.0;

double lorentzian(double nu, double center, double amplitude, double gamma) {
    const double half = 0.5 * gamma;
    const double dx = nu - center;
    return amplitude * (half * half) / (dx * dx + half * half);
}

const char* assignment_for_frequency(double freq_cm1) {
    const double freq = std::abs(freq_cm1);
    if (freq >= 3200.0 && freq <= 3600.0) {
        return "O-H / N-H stretch";
    }
    if (freq >= 2800.0 && freq <= 3200.0) {
        return "C-H stretch";
    }
    if (freq >= 2000.0 && freq <= 2300.0) {
        return "C≡C / C≡N stretch";
    }
    if (freq >= 1600.0 && freq <= 1800.0) {
        return "C=O stretch";
    }
    if (freq >= 1400.0 && freq <= 1600.0) {
        return "C=C stretch / aromatic";
    }
    if (freq >= 1000.0 && freq <= 1400.0) {
        return "C-O / C-N stretch";
    }
    if (freq >= 600.0 && freq <= 1000.0) {
        return "bending / out-of-plane";
    }
    return "skeletal / torsion";
}

ImVec4 stick_color(bool imaginary, bool selected) {
    if (selected) {
        return ImVec4(0.15f, 0.55f, 0.65f, 1.0f);
    }
    if (imaginary) {
        return ImVec4(0.92f, 0.28f, 0.28f, 1.0f);
    }
    return ImVec4(0.85f, 0.78f, 0.22f, 0.95f);
}

}  // namespace

void draw_ir_spectrum_panel(AppState& state,
                            const sbox::backend::JobResult& result,
                            const sbox::chem::MolecularSystem& mol) {
    (void)mol;
    if (!result.has_frequencies || result.frequencies_cm1.empty()) {
        return;
    }

    if (!ImGui::Begin("IR Spectrum")) {
        ImGui::End();
        return;
    }

    if (ImGui::Button("<- Back to Dashboard")) {
        state.property_view = PropertyView::Dashboard;
    }
    ImGui::Separator();

    static float gamma = 10.0f;
    static float max_wavenumber = 4000.0f;
    static bool show_sticks = true;
    static bool invert_y = false;

    ImGui::SliderFloat("Broadening gamma (cm^-1)", &gamma, 1.0f, 50.0f, "%.1f");
    ImGui::SliderFloat("Wavenumber Range", &max_wavenumber, 1000.0f, 5000.0f, "%.0f");
    ImGui::Checkbox("Show Stick Spectrum", &show_sticks);
    ImGui::SameLine();
    ImGui::Checkbox("Invert Y Axis (Transmittance)", &invert_y);

    const std::size_t mode_count = result.frequencies_cm1.size();
    std::vector<float> grid_x(static_cast<std::size_t>(kGridPoints));
    std::vector<float> spectrum(static_cast<std::size_t>(kGridPoints), 0.0f);
    std::vector<float> spectrum_display(static_cast<std::size_t>(kGridPoints), 0.0f);

    const double nu_min = 0.0;
    const double nu_max = std::max(1000.0, static_cast<double>(max_wavenumber));
    const double dnu = (nu_max - nu_min) / static_cast<double>(kGridPoints - 1);

    double max_intensity = 0.0;
    for (std::size_t i = 0; i < mode_count; ++i) {
        const double intensity = i < result.ir_intensities.size() ? result.ir_intensities[i] : 0.0;
        max_intensity = std::max(max_intensity, intensity);
    }

    for (int i = 0; i < kGridPoints; ++i) {
        const double nu = nu_min + dnu * static_cast<double>(i);
        grid_x[static_cast<std::size_t>(i)] = static_cast<float>(nu);

        double total = 0.0;
        for (std::size_t k = 0; k < mode_count; ++k) {
            const double freq = result.frequencies_cm1[k];
            const double intensity = k < result.ir_intensities.size() ? result.ir_intensities[k] : 1.0;
            total += lorentzian(nu, std::abs(freq), intensity, gamma);
        }
        spectrum[static_cast<std::size_t>(i)] = static_cast<float>(total);
        spectrum_display[static_cast<std::size_t>(i)] =
            invert_y && max_intensity > 1.0e-8
                ? static_cast<float>(std::max(0.0, 1.0 - total / std::max(max_intensity, total)))
                : static_cast<float>(total);
    }

    const float plot_y_max = *std::max_element(spectrum_display.begin(), spectrum_display.end());

    if (ImPlot::BeginPlot("IR Spectrum", ImVec2(-1, 300))) {
        ImPlot::SetupAxes("Wavenumber (cm^-1)", invert_y ? "Transmittance" : "Intensity",
                          ImPlotAxisFlags_Invert, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisLimits(ImAxis_X1, max_wavenumber, 0.0, ImGuiCond_Always);
        if (invert_y) {
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0, 1.05, ImGuiCond_Always);
        }

        plot_shaded_styled("IR Envelope",
                           grid_x.data(),
                           spectrum_display.data(),
                           kGridPoints,
                           ImVec4(0.18f, 0.78f, 0.86f, 0.95f),
                           0.25f);
        plot_line_styled("IR Curve",
                         grid_x.data(),
                         spectrum_display.data(),
                         kGridPoints,
                         ImVec4(0.24f, 0.86f, 0.92f, 1.0f),
                         2.0f);

        if (show_sticks) {
            std::vector<double> real_x;
            std::vector<double> real_y;
            std::vector<double> imag_x;
            std::vector<double> imag_y;
            std::vector<double> sel_x;
            std::vector<double> sel_y;
            for (std::size_t i = 0; i < mode_count; ++i) {
                const double freq = std::abs(result.frequencies_cm1[i]);
                const double intensity = i < result.ir_intensities.size() ? result.ir_intensities[i] : 1.0;
                const double y = invert_y ? std::max(0.0, 1.0 - intensity / std::max(max_intensity, 1.0)) : intensity;
                if (static_cast<int>(i) == state.selected_vibrational_mode) {
                    sel_x.push_back(freq);
                    sel_y.push_back(y);
                } else if (result.frequencies_cm1[i] < 0.0) {
                    imag_x.push_back(freq);
                    imag_y.push_back(y);
                } else {
                    real_x.push_back(freq);
                    real_y.push_back(y);
                }
            }

            if (!real_x.empty()) {
                ImPlot::PlotStems("Modes",
                                  real_x.data(),
                                  real_y.data(),
                                  static_cast<int>(real_x.size()),
                                  0.0,
                                  {ImPlotProp_LineColor, stick_color(false, false), ImPlotProp_LineWeight, 1.2f});
            }
            if (!imag_x.empty()) {
                ImPlot::PlotStems("Imaginary Modes",
                                  imag_x.data(),
                                  imag_y.data(),
                                  static_cast<int>(imag_x.size()),
                                  0.0,
                                  {ImPlotProp_LineColor, stick_color(true, false), ImPlotProp_LineWeight, 1.6f});
            }
            if (!sel_x.empty()) {
                ImPlot::PlotStems("Selected Mode",
                                  sel_x.data(),
                                  sel_y.data(),
                                  static_cast<int>(sel_x.size()),
                                  0.0,
                                  {ImPlotProp_LineColor, stick_color(false, true), ImPlotProp_LineWeight, 2.2f});
            }
        }

        std::vector<std::size_t> top_peaks(mode_count);
        for (std::size_t i = 0; i < mode_count; ++i) {
            top_peaks[i] = i;
        }
        std::sort(top_peaks.begin(), top_peaks.end(), [&](std::size_t a, std::size_t b) {
            const double ia = a < result.ir_intensities.size() ? result.ir_intensities[a] : 0.0;
            const double ib = b < result.ir_intensities.size() ? result.ir_intensities[b] : 0.0;
            return ia > ib;
        });
        const int labels = std::min<int>(static_cast<int>(top_peaks.size()), 8);
        for (int n = 0; n < labels; ++n) {
            const std::size_t idx = top_peaks[static_cast<std::size_t>(n)];
            const double freq = std::abs(result.frequencies_cm1[idx]);
            if (freq > max_wavenumber) {
                continue;
            }
            const double intensity = idx < result.ir_intensities.size() ? result.ir_intensities[idx] : 1.0;
            const double y = invert_y ? std::max(0.0, 1.0 - intensity / std::max(max_intensity, 1.0)) : intensity;
            char label[32];
            std::snprintf(label, sizeof(label), "%.0f", freq);
            ImPlot::Annotation(freq,
                               y + (invert_y ? 0.03 : 0.04 * std::max(1.0f, plot_y_max)),
                               ImVec4(0.92f, 0.93f, 0.96f, 0.86f),
                               ImVec2(0.5f, 0.0f),
                               false,
                               "%s",
                               label);
        }

        if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            const ImPlotPoint mouse = ImPlot::GetPlotMousePos();
            int nearest = -1;
            double best = std::numeric_limits<double>::max();
            for (std::size_t i = 0; i < mode_count; ++i) {
                const double d = std::abs(mouse.x - std::abs(result.frequencies_cm1[i]));
                if (d < best) {
                    best = d;
                    nearest = static_cast<int>(i);
                }
            }
            if (nearest >= 0 && best <= 25.0) {
                state.selected_vibrational_mode = nearest;
            }
        }

        ImPlot::EndPlot();
    }

    ImGui::Separator();
    if (ImGui::BeginTable("IRModes", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_ScrollY, ImVec2(0.0f, 220.0f))) {
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 36.0f);
        ImGui::TableSetupColumn("Frequency (cm^-1)", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableSetupColumn("Intensity", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Assignment", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();

        for (std::size_t i = 0; i < mode_count; ++i) {
            const bool imaginary = result.frequencies_cm1[i] < 0.0;
            const bool selected = static_cast<int>(i) == state.selected_vibrational_mode;
            const double freq_abs = std::abs(result.frequencies_cm1[i]);
            const double intensity = i < result.ir_intensities.size() ? result.ir_intensities[i] : 0.0;
            const char* assignment = assignment_for_frequency(result.frequencies_cm1[i]);

            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            char row_label[32];
            std::snprintf(row_label, sizeof(row_label), "%zu", i + 1);
            if (ImGui::Selectable(row_label, selected, ImGuiSelectableFlags_SpanAllColumns)) {
                state.selected_vibrational_mode = static_cast<int>(i);
            }
            ImGui::TableSetColumnIndex(1);
            if (imaginary) {
                ImGui::TextColored(ImVec4(0.92f, 0.28f, 0.28f, 1.0f), "i%.1f", freq_abs);
            } else {
                ImGui::Text("%.1f", freq_abs);
            }
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%.2f", intensity);
            ImGui::TableSetColumnIndex(3);
            ImGui::TextUnformatted(assignment);
        }
        ImGui::EndTable();
    }

    if (state.selected_vibrational_mode >= 0 && state.selected_vibrational_mode < static_cast<int>(mode_count)) {
        const int idx = state.selected_vibrational_mode;
        const double freq = std::abs(result.frequencies_cm1[static_cast<std::size_t>(idx)]);
        const double wavelength_um = freq > 1.0e-8 ? kCmToMicron / freq : 0.0;
        const double intensity = idx < static_cast<int>(result.ir_intensities.size()) ? result.ir_intensities[static_cast<std::size_t>(idx)] : 0.0;
        const char* assignment = assignment_for_frequency(result.frequencies_cm1[static_cast<std::size_t>(idx)]);

        ImGui::Separator();
        ImGui::Text("Mode #%d: nu = %.1f cm^-1 (%.2f um), Intensity = %.2f", idx + 1, freq, wavelength_um, intensity);
        ImGui::Text("Region: %s", assignment);
        if (ImGui::Button("Animate")) {
            state.animate_vibrational_mode = true;  // TODO: normal mode animation pipeline
        }
        if (state.animate_vibrational_mode) {
            ImGui::SameLine();
            ImGui::TextDisabled("TODO: animation");
        }
    }

    ImGui::End();
}

}  // namespace sbox::ui
