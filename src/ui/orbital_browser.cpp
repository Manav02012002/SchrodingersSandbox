#include "ui/orbital_browser.h"

#include "core/elements.h"
#include "core/slater.h"

#include <imgui.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>

namespace sbox::ui {
namespace {

constexpr std::array<const char*, 7> kLLabels = {"s", "p", "d", "f", "g", "h", "i"};

std::string superscript_number(int value) {
    static const std::array<const char*, 10> kDigits = {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};

    if (value == 0) {
        return kDigits[0];
    }

    std::string out;
    if (value < 0) {
        out += "⁻";
        value = -value;
    }

    std::string digits = std::to_string(value);
    for (char ch : digits) {
        out += kDigits[static_cast<std::size_t>(ch - '0')];
    }
    return out;
}

}  // namespace

void draw_orbital_browser(AppState& state) {
    ImGui::SetNextWindowSize(ImVec2(400, 300), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Orbital Browser")) {
        ImGui::End();
        return;
    }

    const auto& element = sbox::elements::get_element(state.selected_Z);
    const auto& config = element.config;

    if (config.empty()) {
        ImGui::TextUnformatted("No orbitals available.");
        ImGui::End();
        return;
    }

    int index = state.selected_orbital_index;
    if (index < 0 || index >= static_cast<int>(config.size())) {
        index = static_cast<int>(config.size()) - 1;
    }

    for (std::size_t i = 0; i < config.size(); ++i) {
        const auto& subshell = config[i];
        const double zeff = sbox::slater::compute_zeff(element.Z, config, subshell.n, subshell.l);
        const double energy = -13.6 * std::pow(zeff / static_cast<double>(subshell.n), 2.0);

        char line[128] = {};
        std::snprintf(line,
                      sizeof(line),
                      "%d%s%s   Zeff=%.2f   E=%.1f eV",
                      subshell.n,
                      kLLabels[static_cast<std::size_t>(std::clamp(subshell.l, 0, 6))],
                      superscript_number(subshell.electrons).c_str(),
                      zeff,
                      energy);

        const bool selected = static_cast<int>(i) == index;
        if (ImGui::Selectable(line, selected)) {
            state.selected_orbital_index = static_cast<int>(i);
            state.selected_m = 0;
            state.needs_update = true;
            index = static_cast<int>(i);
        }
    }

    const auto& selected = config[static_cast<std::size_t>(index)];
    const int l = selected.l;

    if (l > 0) {
        ImGui::Separator();
        ImGui::TextUnformatted("m =");

        for (int m = -l; m <= l; ++m) {
            if (m > -l) {
                ImGui::SameLine();
            }

            bool highlighted = (state.selected_m == m);
            if (highlighted) {
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.28f, 0.55f, 0.94f, 0.95f));
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.34f, 0.62f, 1.0f, 1.0f));
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.24f, 0.50f, 0.88f, 1.0f));
            }

            char m_label[16] = {};
            std::snprintf(m_label, sizeof(m_label), "%d", m);
            if (ImGui::Button(m_label)) {
                state.selected_m = m;
                state.needs_update = true;
            }

            if (highlighted) {
                ImGui::PopStyleColor(3);
            }
        }
    }

    const double zeff = sbox::slater::compute_zeff(element.Z, config, selected.n, selected.l);
    const double energy = -13.6 * std::pow(zeff / static_cast<double>(selected.n), 2.0);

    ImGui::Separator();
    ImGui::BeginChild("##orbital_info", ImVec2(0.0f, 48.0f), true);
    ImGui::Text("n = %d   l = %d   m = %d   Zeff = %.2f   E = %.2f eV",
                selected.n,
                selected.l,
                std::clamp(state.selected_m, -selected.l, selected.l),
                zeff,
                energy);
    ImGui::EndChild();

    ImGui::End();
}

}  // namespace sbox::ui
