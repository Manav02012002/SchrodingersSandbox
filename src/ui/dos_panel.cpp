#include "ui/dos_panel.h"

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
constexpr int kGridCount = 500;
constexpr double kPi = 3.14159265358979323846;

double gaussian(double x, double mu, double sigma) {
    const double dx = x - mu;
    const double inv = 1.0 / (sigma * std::sqrt(2.0 * kPi));
    return inv * std::exp(-(dx * dx) / (2.0 * sigma * sigma));
}

}  // namespace

void draw_dos_panel(AppState& state, const sbox::backend::JobResult& result) {
    if (!result.has_mo_data || result.mo_data.energies.size() == 0) {
        return;
    }

    if (!ImGui::Begin("Density of States")) {
        ImGui::End();
        return;
    }

    if (ImGui::Button("<- Back to Dashboard")) {
        state.property_view = PropertyView::Dashboard;
    }
    ImGui::Separator();

    static float sigma = 0.3f;
    static bool show_occupied_only = false;
    static bool show_levels = true;

    ImGui::SliderFloat("Broadening sigma (eV)", &sigma, 0.05f, 2.0f, "%.2f");
    ImGui::Checkbox("Show Occupied Only", &show_occupied_only);
    ImGui::SameLine();
    ImGui::Checkbox("Show Individual Levels", &show_levels);

    const auto& mo = result.mo_data;
    std::vector<double> energies_ev;
    std::vector<double> occupations;
    energies_ev.reserve(static_cast<std::size_t>(mo.energies.size()));
    occupations.reserve(static_cast<std::size_t>(mo.energies.size()));

    double e_min = 0.0;
    double e_max = 0.0;
    for (int i = 0; i < mo.energies.size(); ++i) {
        const double e = mo.energies(i) * kHartreeToEv;
        energies_ev.push_back(e);
        occupations.push_back(i < mo.occupations.size() ? mo.occupations(i) : 0.0);
        if (i == 0 || e < e_min) {
            e_min = e;
        }
        if (i == 0 || e > e_max) {
            e_max = e;
        }
    }

    e_min -= 5.0 * sigma;
    e_max += 5.0 * sigma;
    const double step = (e_max - e_min) / static_cast<double>(kGridCount - 1);

    std::vector<double> grid(static_cast<std::size_t>(kGridCount));
    for (int i = 0; i < kGridCount; ++i) {
        grid[static_cast<std::size_t>(i)] = e_min + step * static_cast<double>(i);
    }

    std::vector<double> total_dos(static_cast<std::size_t>(kGridCount), 0.0);
    std::vector<double> occupied_dos(static_cast<std::size_t>(kGridCount), 0.0);
    std::vector<double> virtual_dos(static_cast<std::size_t>(kGridCount), 0.0);

    for (std::size_t i = 0; i < energies_ev.size(); ++i) {
        const double occ = occupations[i];
        for (int g = 0; g < kGridCount; ++g) {
            const double contrib = occ * gaussian(grid[static_cast<std::size_t>(g)], energies_ev[i], sigma);
            total_dos[static_cast<std::size_t>(g)] += contrib;
            if (occ > 0.5) {
                occupied_dos[static_cast<std::size_t>(g)] += contrib;
            } else {
                virtual_dos[static_cast<std::size_t>(g)] += contrib;
            }
        }
    }

    const int homo = result.homo_index();
    const int lumo = result.lumo_index();
    const double homo_ev = (homo >= 0 && homo < static_cast<int>(energies_ev.size())) ? energies_ev[static_cast<std::size_t>(homo)] : 0.0;
    const double lumo_ev = (lumo >= 0 && lumo < static_cast<int>(energies_ev.size())) ? energies_ev[static_cast<std::size_t>(lumo)] : 0.0;
    const double fermi_ev = (homo >= 0 && lumo >= 0) ? 0.5 * (homo_ev + lumo_ev) : 0.0;
    const double selected_ev = (state.selected_mo >= 0 && state.selected_mo < static_cast<int>(energies_ev.size()))
        ? energies_ev[static_cast<std::size_t>(state.selected_mo)]
        : 0.0;

    std::vector<float> grid_f(static_cast<std::size_t>(kGridCount));
    std::vector<float> total_f(static_cast<std::size_t>(kGridCount));
    std::vector<float> occ_f(static_cast<std::size_t>(kGridCount));
    std::vector<float> virt_f(static_cast<std::size_t>(kGridCount));
    for (int i = 0; i < kGridCount; ++i) {
        grid_f[static_cast<std::size_t>(i)] = static_cast<float>(grid[static_cast<std::size_t>(i)]);
        total_f[static_cast<std::size_t>(i)] = static_cast<float>(total_dos[static_cast<std::size_t>(i)]);
        occ_f[static_cast<std::size_t>(i)] = static_cast<float>(occupied_dos[static_cast<std::size_t>(i)]);
        virt_f[static_cast<std::size_t>(i)] = static_cast<float>(virtual_dos[static_cast<std::size_t>(i)]);
    }

    if (ImPlot::BeginPlot("Density of States", ImVec2(-1, 300))) {
        ImPlot::SetupAxes("Energy (eV)", "DOS (states/eV)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);

        if (!show_occupied_only) {
            plot_shaded_styled("Occupied DOS",
                               grid_f.data(),
                               occ_f.data(),
                               kGridCount,
                               ImVec4(0.20f, 0.45f, 0.90f, 0.85f),
                               0.22f);
            plot_shaded_styled("Virtual DOS",
                               grid_f.data(),
                               virt_f.data(),
                               kGridCount,
                               ImVec4(0.55f, 0.58f, 0.66f, 0.85f),
                               0.18f);
        }

        plot_shaded_styled("Total DOS",
                           grid_f.data(),
                           total_f.data(),
                           kGridCount,
                           ImVec4(0.18f, 0.70f, 0.78f, 0.95f),
                           0.25f);
        plot_line_styled("Total DOS Line",
                         grid_f.data(),
                         total_f.data(),
                         kGridCount,
                         ImVec4(0.18f, 0.78f, 0.86f, 1.0f),
                         2.0f);

        if (show_levels) {
            std::vector<double> xs_occ;
            std::vector<double> ys_occ;
            std::vector<double> xs_virt;
            std::vector<double> ys_virt;
            std::vector<double> xs_sel;
            std::vector<double> ys_sel;

            for (std::size_t i = 0; i < energies_ev.size(); ++i) {
                const double occ = occupations[i];
                const double x = energies_ev[i];
                const double y = 0.5;
                if (static_cast<int>(i) == state.selected_mo) {
                    xs_sel.push_back(x);
                    ys_sel.push_back(y);
                } else if (occ > 0.5) {
                    xs_occ.push_back(x);
                    ys_occ.push_back(y);
                } else {
                    xs_virt.push_back(x);
                    ys_virt.push_back(y);
                }
            }

            if (!xs_occ.empty()) {
                ImPlot::PlotStems("Occupied Levels",
                                  xs_occ.data(),
                                  ys_occ.data(),
                                  static_cast<int>(xs_occ.size()),
                                  0.0,
                                  {ImPlotProp_LineColor, ImVec4(0.20f, 0.45f, 0.90f, 0.95f), ImPlotProp_LineWeight, 1.4f});
            }
            if (!xs_virt.empty() && !show_occupied_only) {
                ImPlot::PlotStems("Virtual Levels",
                                  xs_virt.data(),
                                  ys_virt.data(),
                                  static_cast<int>(xs_virt.size()),
                                  0.0,
                                  {ImPlotProp_LineColor, ImVec4(0.58f, 0.60f, 0.66f, 0.95f), ImPlotProp_LineWeight, 1.2f});
            }
            if (!xs_sel.empty()) {
                ImPlot::PlotStems("Selected Level",
                                  xs_sel.data(),
                                  ys_sel.data(),
                                  static_cast<int>(xs_sel.size()),
                                  0.0,
                                  {ImPlotProp_LineColor, ImVec4(0.15f, 0.55f, 0.65f, 1.0f), ImPlotProp_LineWeight, 2.2f});
            }
        }

        if (homo >= 0 && lumo >= 0) {
            const double fx[] = {fermi_ev, fermi_ev};
            const double fy[] = {0.0, *std::max_element(total_dos.begin(), total_dos.end())};
            ImPlot::PlotLine("E_F", fx, fy, 2, {ImPlotProp_LineColor, ImVec4(0.95f, 0.82f, 0.22f, 0.9f), ImPlotProp_LineWeight, 1.5f});
        }
        if (state.selected_mo >= 0 && state.selected_mo < static_cast<int>(energies_ev.size())) {
            const double sx[] = {selected_ev, selected_ev};
            const double sy[] = {0.0, *std::max_element(total_dos.begin(), total_dos.end())};
            ImPlot::PlotLine("Selected MO", sx, sy, 2, {ImPlotProp_LineColor, ImVec4(0.15f, 0.55f, 0.65f, 0.9f), ImPlotProp_LineWeight, 1.5f});
        }

        if (ImPlot::IsPlotHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            const ImPlotPoint mouse = ImPlot::GetPlotMousePos();
            int nearest = -1;
            double best = 1e9;
            for (std::size_t i = 0; i < energies_ev.size(); ++i) {
                const double d = std::abs(mouse.x - energies_ev[i]);
                if (d < best) {
                    best = d;
                    nearest = static_cast<int>(i);
                }
            }
            if (nearest >= 0) {
                state.selected_mo = nearest;
            }
        }

        ImPlot::EndPlot();
    }

    int occupied = 0;
    for (double occ : occupations) {
        if (occ > 0.5) {
            ++occupied;
        }
    }
    const int total = static_cast<int>(energies_ev.size());
    const int virt = total - occupied;

    ImGui::Separator();
    ImGui::Text("Total orbitals: %d (%d occupied, %d virtual)", total, occupied, virt);
    if (homo >= 0) {
        ImGui::Text("HOMO energy: %.2f eV", homo_ev);
    }
    if (lumo >= 0) {
        ImGui::Text("LUMO energy: %.2f eV", lumo_ev);
    }
    if (homo >= 0 && lumo >= 0) {
        ImGui::Text("Band gap: %.2f eV", lumo_ev - homo_ev);
        ImGui::Text("Fermi energy (midgap): %.2f eV", fermi_ev);
    }

    ImGui::End();
}

}  // namespace sbox::ui
