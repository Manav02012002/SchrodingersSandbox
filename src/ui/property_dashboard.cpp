#include "ui/property_dashboard.h"

#include "core/elements.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace sbox::ui {

namespace {

constexpr double kHartreeToEv = 27.2114;

std::string polarity_text(double dipole) {
    if (dipole < 0.5) {
        return "Nonpolar molecule";
    }
    if (dipole < 1.5) {
        return "Weakly polar";
    }
    if (dipole < 3.0) {
        return "Polar molecule";
    }
    return "Highly polar";
}

void draw_card(const char* id,
               const char* title,
               const std::string& value,
               const std::string& subtext,
               float width,
               float height,
               bool enabled,
               AppState& state,
               PropertyView view) {
    ImGui::BeginChild(id, ImVec2(width, height), true);
    ImGui::SetWindowFontScale(1.08f);
    ImGui::TextUnformatted(title);
    ImGui::SetWindowFontScale(1.0f);
    ImGui::Separator();
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.18f, 0.78f, 0.86f, 1.0f));
    ImGui::TextWrapped("%s", value.c_str());
    ImGui::PopStyleColor();
    if (!subtext.empty()) {
        ImGui::TextDisabled("%s", subtext.c_str());
    }
    ImGui::Dummy(ImVec2(0.0f, 6.0f));
    if (!enabled) {
        ImGui::BeginDisabled();
    }
    if (ImGui::Button("Details ->")) {
        if (view == PropertyView::Dashboard) {
            ImGui::SetWindowFocus("Results");
        } else {
            state.property_view = view;
        }
    }
    if (!enabled) {
        ImGui::EndDisabled();
    }
    ImGui::EndChild();
}

}  // namespace

void draw_property_dashboard(AppState& state,
                             const sbox::backend::JobResult& result,
                             const sbox::chem::MolecularSystem& mol) {
    if (!result.converged()) {
        return;
    }

    if (!ImGui::Begin("Properties")) {
        ImGui::End();
        return;
    }

    const float spacing = ImGui::GetStyle().ItemSpacing.x;
    const float panel_width = ImGui::GetContentRegionAvail().x;
    const int columns = panel_width > 760.0f ? 3 : 2;
    const float card_width = (panel_width - spacing * static_cast<float>(columns - 1)) / static_cast<float>(columns);
    const float card_height = 118.0f;

    int card_in_row = 0;
    auto place_card = [&](const char* id,
                          const char* title,
                          const std::string& value,
                          const std::string& subtext,
                          bool enabled,
                          PropertyView view) {
        draw_card(id, title, value, subtext, card_width, card_height, enabled, state, view);
        ++card_in_row;
        if (card_in_row < columns) {
            ImGui::SameLine();
        } else {
            card_in_row = 0;
        }
    };

    char line[256];
    std::snprintf(line, sizeof(line), "%.4f Hartree", result.total_energy);
    char sub[256];
    std::snprintf(sub, sizeof(sub), "(%.2f eV)", result.total_energy * kHartreeToEv);
    place_card("##card_energy", "Energy", line, sub, true, PropertyView::Dashboard);

    const double dipole_mag = result.dipole_moment.norm();
    std::snprintf(line, sizeof(line), "%.2f Debye", dipole_mag);
    place_card("##card_dipole", "Dipole Moment", line, polarity_text(dipole_mag), true, PropertyView::ESPControls);

    const int homo = result.homo_index();
    const int lumo = result.lumo_index();
    if (result.has_mo_data && homo >= 0 && lumo >= 0) {
        const double homo_ev = result.mo_data.energies(homo) * kHartreeToEv;
        const double lumo_ev = result.mo_data.energies(lumo) * kHartreeToEv;
        std::snprintf(line, sizeof(line), "Gap = %.2f eV", lumo_ev - homo_ev);
        std::snprintf(sub, sizeof(sub), "HOMO: %.2f eV, LUMO: %.2f eV", homo_ev, lumo_ev);
        place_card("##card_frontier", "Frontier Orbitals", line, sub, true, PropertyView::MoDiagram);
    } else {
        place_card("##card_frontier", "Frontier Orbitals", "No MO data", "", false, PropertyView::MoDiagram);
    }

    if (!result.mulliken_charges.empty()) {
        const auto max_it = std::max_element(result.mulliken_charges.begin(), result.mulliken_charges.end());
        const auto min_it = std::min_element(result.mulliken_charges.begin(), result.mulliken_charges.end());
        const int max_idx = static_cast<int>(std::distance(result.mulliken_charges.begin(), max_it));
        const int min_idx = static_cast<int>(std::distance(result.mulliken_charges.begin(), min_it));
        const std::string max_sym =
            max_idx >= 0 && max_idx < mol.num_atoms() ? sbox::elements::get_element(mol.atom(max_idx).Z).symbol : "?";
        const std::string min_sym =
            min_idx >= 0 && min_idx < mol.num_atoms() ? sbox::elements::get_element(mol.atom(min_idx).Z).symbol : "?";
        std::snprintf(line, sizeof(line), "Max: %+.2f (%s), Min: %+.2f (%s)", *max_it, max_sym.c_str(), *min_it, min_sym.c_str());
        place_card("##card_charges", "Partial Charges", line, "Mulliken population analysis", true, PropertyView::Population);
    } else {
        place_card("##card_charges", "Partial Charges", "Unavailable", "", false, PropertyView::Population);
    }

    if (result.mayer_bond_orders.size() > 0) {
        double strongest = 0.0;
        int bi = -1;
        int bj = -1;
        int count = 0;
        for (int i = 0; i < result.mayer_bond_orders.rows(); ++i) {
            for (int j = i + 1; j < result.mayer_bond_orders.cols(); ++j) {
                const double bo = result.mayer_bond_orders(i, j);
                if (bo > 0.5) {
                    ++count;
                }
                if (bo > strongest) {
                    strongest = bo;
                    bi = i;
                    bj = j;
                }
            }
        }
        if (bi >= 0 && bj >= 0 && bi < mol.num_atoms() && bj < mol.num_atoms()) {
            std::snprintf(line, sizeof(line), "%d bonds, strongest: %s-%s (%.1f)",
                          count,
                          sbox::elements::get_element(mol.atom(bi).Z).symbol,
                          sbox::elements::get_element(mol.atom(bj).Z).symbol,
                          strongest);
        } else {
            std::snprintf(line, sizeof(line), "%d bonds", count);
        }
        place_card("##card_bond_orders", "Bond Orders", line, "", true, PropertyView::BondOrders);
    } else {
        place_card("##card_bond_orders", "Bond Orders", "Unavailable", "", false, PropertyView::BondOrders);
    }

    if (result.has_frequencies && !result.frequencies_cm1.empty()) {
        int imaginary = 0;
        double strongest_freq = 0.0;
        double strongest_intensity = -1.0;
        for (std::size_t i = 0; i < result.frequencies_cm1.size(); ++i) {
            const double freq = result.frequencies_cm1[i];
            if (freq < 0.0) {
                ++imaginary;
            }
            const double intensity = i < result.ir_intensities.size() ? result.ir_intensities[i] : 0.0;
            if (intensity > strongest_intensity) {
                strongest_intensity = intensity;
                strongest_freq = std::abs(freq);
            }
        }
        std::snprintf(line, sizeof(line), "%zu modes, strongest: %.0f cm^-1", result.frequencies_cm1.size(), strongest_freq);
        if (imaginary == 0) {
            std::snprintf(sub, sizeof(sub), "No imaginary frequencies");
        } else {
            std::snprintf(sub, sizeof(sub), "WARNING: %d imaginary frequencies", imaginary);
        }
        place_card("##card_ir", "Vibrational Modes", line, sub, true, PropertyView::IRSpectrum);
    } else {
        place_card("##card_ir", "Vibrational Modes", "Unavailable", "", false, PropertyView::IRSpectrum);
    }

    if (card_in_row != 0) {
        ImGui::NewLine();
    }

    ImGui::End();
}

}  // namespace sbox::ui
