#include "ui/mo_diagram_plot.h"

#include "ui/plot_utils.h"

#include <imgui.h>
#include <implot.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace sbox::ui {

namespace {

constexpr double kHartreeToEv = 27.2114;

struct VisibleOrbital {
    int index = -1;
    double energy_h = 0.0;
    double energy_ev = 0.0;
    double occupation = 0.0;
};

void draw_arrow(ImDrawList* draw_list, const ImVec2& center, float size, ImU32 color, bool up) {
    const float half_w = size * 0.3f;
    const float body_h = size * 0.65f;
    if (up) {
        draw_list->AddLine(ImVec2(center.x, center.y + body_h * 0.5f),
                           ImVec2(center.x, center.y - body_h * 0.2f),
                           color,
                           1.2f);
        draw_list->AddTriangleFilled(ImVec2(center.x, center.y - body_h * 0.55f),
                                     ImVec2(center.x - half_w, center.y - body_h * 0.1f),
                                     ImVec2(center.x + half_w, center.y - body_h * 0.1f),
                                     color);
    } else {
        draw_list->AddLine(ImVec2(center.x, center.y - body_h * 0.5f),
                           ImVec2(center.x, center.y + body_h * 0.2f),
                           color,
                           1.2f);
        draw_list->AddTriangleFilled(ImVec2(center.x, center.y + body_h * 0.55f),
                                     ImVec2(center.x - half_w, center.y + body_h * 0.1f),
                                     ImVec2(center.x + half_w, center.y + body_h * 0.1f),
                                     color);
    }
}

std::string mo_kind_label(int mo_index, int homo, int lumo) {
    if (mo_index == homo) {
        return "HOMO";
    }
    if (mo_index == lumo) {
        return "LUMO";
    }
    if (homo >= 0 && mo_index < homo) {
        return "HOMO-" + std::to_string(homo - mo_index);
    }
    if (lumo >= 0 && mo_index > lumo) {
        return "LUMO+" + std::to_string(mo_index - lumo);
    }
    return "MO";
}

}  // namespace

void draw_mo_diagram_plot(AppState& state, const sbox::backend::JobResult& result) {
    if (!result.has_mo_data) {
        return;
    }

    const auto& mo_data = result.mo_data;
    if (mo_data.energies.size() == 0) {
        return;
    }

    if (!ImGui::Begin("MO Energy Diagram")) {
        ImGui::End();
        return;
    }

    if (ImGui::Button("<- Back to Dashboard")) {
        state.property_view = PropertyView::Dashboard;
    }
    ImGui::Separator();

    // TODO: Split alpha/beta columns for unrestricted calculations when separate spin-channel MO data is available.
    const int homo = result.homo_index();
    const int lumo = result.lumo_index();

    int total_mo = mo_data.energies.size();
    int start = 0;
    int end = total_mo;
    if (total_mo > 12 && homo >= 0) {
        start = std::max(0, homo - 5);
        end = std::min(total_mo, (lumo >= 0 ? lumo + 6 : homo + 7));
        if (end - start < 12) {
            start = std::max(0, end - 12);
            end = std::min(total_mo, start + 12);
        }
    }

    std::vector<VisibleOrbital> orbitals;
    orbitals.reserve(static_cast<std::size_t>(end - start));
    double y_min = 0.0;
    double y_max = 0.0;
    for (int i = start; i < end; ++i) {
        VisibleOrbital orbital;
        orbital.index = i;
        orbital.energy_h = mo_data.energies(i);
        orbital.energy_ev = orbital.energy_h * kHartreeToEv;
        orbital.occupation = i < mo_data.occupations.size() ? mo_data.occupations(i) : 0.0;
        orbitals.push_back(orbital);
        if (i == start || orbital.energy_ev < y_min) {
            y_min = orbital.energy_ev;
        }
        if (i == start || orbital.energy_ev > y_max) {
            y_max = orbital.energy_ev;
        }
    }

    y_min -= 2.0;
    y_max += 2.0;
    if (state.selected_mo < 0 && homo >= 0) {
        state.selected_mo = homo;
    }
    state.selected_mo = std::clamp(state.selected_mo, 0, total_mo - 1);

    if (ImPlot::BeginPlot("##mo_diagram", ImVec2(-1, 400), ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText)) {
        ImPlot::SetupAxes(nullptr, "Energy (eV)", ImPlotAxisFlags_NoTickLabels | ImPlotAxisFlags_NoTickMarks | ImPlotAxisFlags_NoLabel, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxisLimits(ImAxis_X1, 0.0, 1.0, ImPlotCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, y_min, y_max, ImPlotCond_Always);

        ImDrawList* draw_list = ImPlot::GetPlotDrawList();
        const ImPlotPoint gap_a(0.86, homo >= 0 ? mo_data.energies(homo) * kHartreeToEv : 0.0);
        const ImPlotPoint gap_b(0.86, lumo >= 0 ? mo_data.energies(lumo) * kHartreeToEv : 0.0);

        for (const VisibleOrbital& orbital : orbitals) {
            const bool occupied = orbital.occupation > 0.5;
            const bool selected = orbital.index == state.selected_mo;
            ImVec4 color = orbital_color(orbital.index, homo);
            float thickness = 2.0f;
            if (orbital.index == homo) {
                color = ImVec4(0.20f, 0.62f, 0.95f, 1.0f);
                thickness = 3.0f;
            } else if (orbital.index == lumo) {
                color = ImVec4(0.95f, 0.58f, 0.18f, 1.0f);
                thickness = 3.0f;
            }
            if (selected) {
                color = ImVec4(0.15f, 0.55f, 0.65f, 1.0f);
                thickness = 4.0f;
            }

            const ImPlotPoint a(0.22, orbital.energy_ev);
            const ImPlotPoint b(0.72, orbital.energy_ev);
            const ImVec2 pa = ImPlot::PlotToPixels(a);
            const ImVec2 pb = ImPlot::PlotToPixels(b);
            draw_list->AddLine(pa, pb, ImGui::GetColorU32(color), thickness);

            const ImU32 arrow_color = ImGui::GetColorU32(color);
            if (orbital.occupation >= 1.0) {
                draw_arrow(draw_list, ImVec2(pa.x + 14.0f, pa.y - 6.0f), 10.0f, arrow_color, true);
            }
            if (orbital.occupation >= 2.0) {
                draw_arrow(draw_list, ImVec2(pa.x + 24.0f, pa.y + 6.0f), 10.0f, arrow_color, false);
            }

            const std::string label = std::to_string(orbital.index) + "  " + mo_kind_label(orbital.index, homo, lumo);
            draw_list->AddText(ImVec2(pb.x + 10.0f, pb.y - 8.0f), ImGui::GetColorU32(color), label.c_str());
        }

        if (homo >= 0 && lumo >= 0 && lumo < mo_data.energies.size()) {
            const ImVec2 pa = ImPlot::PlotToPixels(gap_a);
            const ImVec2 pb = ImPlot::PlotToPixels(gap_b);
            const ImU32 gap_color = ImGui::GetColorU32(ImVec4(0.92f, 0.92f, 0.96f, 0.9f));
            const int dash_count = 10;
            for (int i = 0; i < dash_count; ++i) {
                const float t0 = static_cast<float>(i) / static_cast<float>(dash_count);
                const float t1 = std::min(1.0f, t0 + 0.05f);
                draw_list->AddLine(ImVec2(pa.x, pa.y + (pb.y - pa.y) * t0),
                                   ImVec2(pa.x, pa.y + (pb.y - pa.y) * t1),
                                   gap_color,
                                   1.5f);
            }
            draw_arrow(draw_list, ImVec2(pa.x, pa.y + 8.0f), 10.0f, gap_color, false);
            draw_arrow(draw_list, ImVec2(pb.x, pb.y - 8.0f), 10.0f, gap_color, true);
            const double gap_ev = result.homo_lumo_gap_eV();
            char gap_label[64];
            std::snprintf(gap_label, sizeof(gap_label), "Gap = %.2f eV", gap_ev);
            draw_list->AddText(ImVec2(pa.x + 10.0f, 0.5f * (pa.y + pb.y) - 8.0f), gap_color, gap_label);
        }

        if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            const ImPlotPoint mouse = ImPlot::GetPlotMousePos();
            int closest = -1;
            double best = 0.3;
            for (const VisibleOrbital& orbital : orbitals) {
                const double d = std::abs(orbital.energy_ev - mouse.y);
                if (d < best) {
                    best = d;
                    closest = orbital.index;
                }
            }
            if (closest >= 0) {
                state.selected_mo = closest;
            }
        }

        ImPlot::EndPlot();
    }

    if (state.selected_mo >= 0 && state.selected_mo < mo_data.energies.size()) {
        const double selected_h = mo_data.energies(state.selected_mo);
        const double selected_ev = selected_h * kHartreeToEv;
        const double occ = state.selected_mo < mo_data.occupations.size() ? mo_data.occupations(state.selected_mo) : 0.0;
        const std::string type = mo_kind_label(state.selected_mo, homo, lumo);

        ImGui::Separator();
        ImGui::Text("Orbital #%d", state.selected_mo);
        ImGui::Text("Energy: %.5f Hartree (%.2f eV)", selected_h, selected_ev);
        ImGui::Text("Occupation: %.1f", occ);
        ImGui::Text("Type: %s", type.c_str());

        if ((state.selected_mo == homo || state.selected_mo == lumo) && result.homo_lumo_gap_eV() > 1.0e-8) {
            const double gap_ev = result.homo_lumo_gap_eV();
            const double wavelength_nm = 1240.0 / gap_ev;
            ImGui::Text("HOMO-LUMO Gap: %.2f eV (%.1f nm)", gap_ev, wavelength_nm);
        }
    }

    ImGui::End();
}

}  // namespace sbox::ui
