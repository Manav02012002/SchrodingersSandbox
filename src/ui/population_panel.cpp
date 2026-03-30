#include "ui/population_panel.h"

#include "core/elements.h"
#include "ui/plot_utils.h"

#include <imgui.h>
#include <implot.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

namespace sbox::ui {

namespace {

ImVec4 charge_color(double charge) {
    if (charge > 0.01) {
        return ImVec4(0.92f, 0.35f, 0.35f, 1.0f);
    }
    if (charge < -0.01) {
        return ImVec4(0.35f, 0.55f, 0.95f, 1.0f);
    }
    return ImVec4(0.82f, 0.82f, 0.84f, 1.0f);
}

void draw_charge_bar(const char* id, double charge, double max_abs) {
    const float width = ImGui::GetContentRegionAvail().x;
    const float height = 14.0f;
    const ImVec2 pos = ImGui::GetCursorScreenPos();
    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    const ImU32 bg = IM_COL32(42, 48, 58, 160);
    draw_list->AddRectFilled(pos, ImVec2(pos.x + width, pos.y + height), bg, 3.0f);

    if (max_abs > 1.0e-8 && std::abs(charge) > 1.0e-8) {
        const float half = width * 0.5f;
        const float mag = static_cast<float>(std::abs(charge) / max_abs);
        const float fill = half * std::min(1.0f, mag);
        const ImU32 color = ImGui::GetColorU32(charge_color(charge));
        if (charge >= 0.0) {
            draw_list->AddRectFilled(ImVec2(pos.x + half, pos.y),
                                     ImVec2(pos.x + half + fill, pos.y + height),
                                     color,
                                     3.0f);
        } else {
            draw_list->AddRectFilled(ImVec2(pos.x + half - fill, pos.y),
                                     ImVec2(pos.x + half, pos.y + height),
                                     color,
                                     3.0f);
        }
        draw_list->AddLine(ImVec2(pos.x + half, pos.y), ImVec2(pos.x + half, pos.y + height), IM_COL32(235, 238, 245, 100), 1.0f);
    }

    ImGui::Dummy(ImVec2(width, height));
    ImGui::PushID(id);
    ImGui::InvisibleButton("##chargebar", ImVec2(width, 0.0f));
    ImGui::PopID();
}

void draw_charge_table(const char* table_id,
                       const char* heading,
                       const std::vector<double>& charges,
                       const sbox::chem::MolecularSystem& mol) {
    if (charges.empty()) {
        return;
    }

    const double max_abs = std::max(1.0e-8, [&]() {
        double value = 0.0;
        for (double q : charges) {
            value = std::max(value, std::abs(q));
        }
        return value;
    }());

    ImGui::TextUnformatted(heading);
    if (ImGui::BeginTable(table_id, 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 36.0f);
        ImGui::TableSetupColumn("Element", ImGuiTableColumnFlags_WidthFixed, 70.0f);
        ImGui::TableSetupColumn("Charge", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Bar", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();

        for (int i = 0; i < mol.num_atoms() && i < static_cast<int>(charges.size()); ++i) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%d", i + 1);
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%s", sbox::elements::get_element(mol.atom(i).Z).symbol);
            ImGui::TableSetColumnIndex(2);
            ImGui::TextColored(charge_color(charges[static_cast<std::size_t>(i)]), "%+.4f", charges[static_cast<std::size_t>(i)]);
            ImGui::TableSetColumnIndex(3);
            draw_charge_bar(table_id, charges[static_cast<std::size_t>(i)], max_abs);
        }
        ImGui::EndTable();
    }

    const double sum = std::accumulate(charges.begin(), charges.end(), 0.0);
    int max_pos_idx = -1;
    int max_neg_idx = -1;
    for (int i = 0; i < static_cast<int>(charges.size()); ++i) {
        if (max_pos_idx < 0 || charges[static_cast<std::size_t>(i)] > charges[static_cast<std::size_t>(max_pos_idx)]) {
            max_pos_idx = i;
        }
        if (max_neg_idx < 0 || charges[static_cast<std::size_t>(i)] < charges[static_cast<std::size_t>(max_neg_idx)]) {
            max_neg_idx = i;
        }
    }

    ImGui::Text("Sum of charges: %.4f", sum);
    if (max_pos_idx >= 0 && max_pos_idx < mol.num_atoms()) {
        ImGui::Text("Max positive: %s%d = %+.4f",
                    sbox::elements::get_element(mol.atom(max_pos_idx).Z).symbol,
                    max_pos_idx + 1,
                    charges[static_cast<std::size_t>(max_pos_idx)]);
    }
    if (max_neg_idx >= 0 && max_neg_idx < mol.num_atoms()) {
        ImGui::Text("Max negative: %s%d = %+.4f",
                    sbox::elements::get_element(mol.atom(max_neg_idx).Z).symbol,
                    max_neg_idx + 1,
                    charges[static_cast<std::size_t>(max_neg_idx)]);
    }
}

void draw_colormap_legend(double max_abs) {
    const float width = ImGui::GetContentRegionAvail().x;
    const float height = 18.0f;
    const ImVec2 pos = ImGui::GetCursorScreenPos();
    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    const int steps = 64;
    for (int i = 0; i < steps; ++i) {
        const float t0 = static_cast<float>(i) / static_cast<float>(steps);
        const float t1 = static_cast<float>(i + 1) / static_cast<float>(steps);
        ImVec4 color;
        const float mid = 0.5f;
        if (t0 < mid) {
            const float u = t0 / mid;
            color = ImVec4(0.35f + 0.65f * u, 0.55f + 0.45f * u, 0.95f + 0.05f * u, 1.0f);
        } else {
            const float u = (t0 - mid) / mid;
            color = ImVec4(1.0f, 1.0f - 0.65f * u, 1.0f - 0.65f * u, 1.0f);
        }
        draw_list->AddRectFilled(ImVec2(pos.x + width * t0, pos.y),
                                 ImVec2(pos.x + width * t1, pos.y + height),
                                 ImGui::GetColorU32(color));
    }
    draw_list->AddRect(pos, ImVec2(pos.x + width, pos.y + height), IM_COL32(210, 214, 222, 120), 3.0f);
    ImGui::Dummy(ImVec2(width, height + 4.0f));
    ImGui::Text("-%.4f          0          +%.4f", max_abs, max_abs);
}

}  // namespace

void draw_population_panel(AppState& state,
                           const sbox::backend::JobResult& result,
                           const sbox::chem::MolecularSystem& mol) {
    if (!result.converged() || result.mulliken_charges.empty()) {
        return;
    }

    if (!ImGui::Begin("Population Analysis")) {
        ImGui::End();
        return;
    }

    if (ImGui::Button("<- Back to Dashboard")) {
        state.property_view = PropertyView::Dashboard;
    }
    ImGui::Separator();

    draw_charge_table("MullikenChargesTable", "Mulliken Charges", result.mulliken_charges, mol);

    if (!result.lowdin_charges.empty()) {
        ImGui::Separator();
        draw_charge_table("LowdinChargesTable", "Lowdin Charges", result.lowdin_charges, mol);
    }

    ImGui::Separator();
    std::vector<float> xs;
    std::vector<float> pos;
    std::vector<float> neg;
    xs.reserve(result.mulliken_charges.size());
    pos.reserve(result.mulliken_charges.size());
    neg.reserve(result.mulliken_charges.size());
    double max_abs = 0.0;
    double min_q = 0.0;
    double max_q = 0.0;
    for (std::size_t i = 0; i < result.mulliken_charges.size(); ++i) {
        const float x = static_cast<float>(i + 1);
        const float q = static_cast<float>(result.mulliken_charges[i]);
        xs.push_back(x);
        pos.push_back(q > 0.0f ? q : 0.0f);
        neg.push_back(q < 0.0f ? q : 0.0f);
        max_abs = std::max(max_abs, std::abs(result.mulliken_charges[i]));
        if (i == 0 || result.mulliken_charges[i] < min_q) {
            min_q = result.mulliken_charges[i];
        }
        if (i == 0 || result.mulliken_charges[i] > max_q) {
            max_q = result.mulliken_charges[i];
        }
    }

    if (ImPlot::BeginPlot("Mulliken Charges", ImVec2(-1, 200))) {
        ImPlot::SetupAxes("Atom Index", "Charge", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        plot_stems_styled("Positive", xs.data(), pos.data(), static_cast<int>(xs.size()), ImVec4(0.92f, 0.35f, 0.35f, 1.0f), 2.0f);
        plot_stems_styled("Negative", xs.data(), neg.data(), static_cast<int>(xs.size()), ImVec4(0.35f, 0.55f, 0.95f, 1.0f), 2.0f);
        const double ref_x[] = {0.0, static_cast<double>(xs.size()) + 1.0};
        const double ref_y[] = {0.0, 0.0};
        ImPlot::PlotLine("Zero", ref_x, ref_y, 2, {ImPlotProp_LineColor, ImVec4(0.85f, 0.85f, 0.88f, 0.60f), ImPlotProp_LineWeight, 1.0f});
        ImPlot::EndPlot();
    }

    ImGui::Separator();
    ImGui::Text("Electronegativity range: %.2f", max_q - min_q);
    ImGui::Text("Dipole moment: %.4f Debye", result.dipole_moment.norm());

    ImGui::Separator();
    ImGui::TextUnformatted("Charge Colormap");
    draw_colormap_legend(std::max(max_abs, 1.0e-6));

    ImGui::End();
}

}  // namespace sbox::ui
