#include "ui/nci_panel.h"

#include <imgui.h>
#include <implot.h>

namespace sbox::ui {

void draw_nci_panel(AppState& state, const sbox::backend::JobResult& result) {
    if (!result.has_density_cube) {
        return;
    }

    if (!ImGui::Begin("Non-Covalent Interactions")) {
        ImGui::End();
        return;
    }

    ImGui::Checkbox("Show NCI Surface", &state.show_nci);
    ImGui::SliderFloat("RDG Isosurface", &state.nci_rdg_iso, 0.1f, 1.0f, "%.2f");
    ImGui::SliderFloat("rho Cutoff", &state.nci_rho_cutoff, 0.01f, 0.1f, "%.3f");
    ImGui::SliderFloat("Color Range", &state.nci_color_range, 0.01f, 0.1f, "%.3f");

    if (ImGui::Button("Compute NCI")) {
        state.nci_compute_requested = true;
    }

    ImGui::Separator();
    ImGui::TextUnformatted("Blue: strong attractive (H-bonds)");
    ImGui::TextUnformatted("Green: weak (van der Waals)");
    ImGui::TextUnformatted("Red: strong repulsive (steric)");

    if (!state.nci_plot_sign_rho.empty() && state.nci_plot_sign_rho.size() == state.nci_plot_rdg.size()) {
        if (ImPlot::BeginPlot("NCI Diagnostic Plot", ImVec2(-1.0f, 240.0f))) {
            ImPlot::SetupAxes("sign(lambda2)*rho", "s(r)");
            ImPlot::SetupAxesLimits(-state.nci_color_range, state.nci_color_range, 0.0, 1.0, ImPlotCond_Always);
            ImPlot::PlotScatter("NCI",
                                state.nci_plot_sign_rho.data(),
                                state.nci_plot_rdg.data(),
                                static_cast<int>(state.nci_plot_rdg.size()));
            ImPlot::EndPlot();
        }
    } else {
        ImGui::TextDisabled("Compute NCI to populate the diagnostic plot.");
    }

    ImGui::End();
}

}  // namespace sbox::ui
