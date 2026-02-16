#include "ui/properties_panel.h"

#include "core/elements.h"
#include "core/hydrogen.h"

#include <imgui.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

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

std::string config_to_string(const sbox::slater::ElectronConfig& config) {
    std::string out;
    for (std::size_t i = 0; i < config.size(); ++i) {
        const auto& subshell = config[i];
        if (i > 0) {
            out += " ";
        }
        out += std::to_string(subshell.n);
        out += kLLabels[static_cast<std::size_t>(std::clamp(subshell.l, 0, 6))];
        out += superscript_number(subshell.electrons);
    }
    return out;
}

}  // namespace

void draw_properties(AppState& state) {
    ImGui::SetNextWindowSize(ImVec2(400, 300), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Properties")) {
        ImGui::End();
        return;
    }

    const auto& element = sbox::elements::get_element(state.selected_Z);
    const auto& config = element.config;

    if (config.empty()) {
        ImGui::TextUnformatted("No properties available.");
        ImGui::End();
        return;
    }

    int index = state.selected_orbital_index;
    if (index < 0 || index >= static_cast<int>(config.size())) {
        index = static_cast<int>(config.size()) - 1;
    }

    const auto& subshell = config[static_cast<std::size_t>(index)];
    const int m = std::clamp(state.selected_m, -subshell.l, subshell.l);
    const double zeff = static_cast<double>(state.current_Zeff);
    const double energy = -13.6 * std::pow(zeff / static_cast<double>(subshell.n), 2.0);

    ImGui::SetWindowFontScale(1.2f);
    ImGui::Text("%s -- %s", element.symbol, element.name);
    ImGui::SetWindowFontScale(1.0f);

    ImGui::Text("Z = %d     %s", element.Z, element.category);

    const std::string config_text = config_to_string(config);
    ImGui::TextWrapped("Configuration: %s", config_text.c_str());
    ImGui::Text("Mass: %.6g u", element.atomic_mass);

    ImGui::Separator();
    ImGui::Text("Orbital: %d%s%s",
                subshell.n,
                kLLabels[static_cast<std::size_t>(std::clamp(subshell.l, 0, 6))],
                superscript_number(subshell.electrons).c_str());
    ImGui::Text("Quantum numbers: n=%d, l=%d, m=%d", subshell.n, subshell.l, m);
    ImGui::Text("Effective nuclear charge: Zeff = %.4f", zeff);
    ImGui::Text("Orbital energy: %.4f eV", energy);

    ImGui::Separator();
    ImGui::TextUnformatted("Radial probability density P(r) = r²|R(r)|²");

    constexpr int kSamples = 300;
    std::vector<float> values(static_cast<std::size_t>(kSamples), 0.0f);
    const double max_r = 4.0 * static_cast<double>(subshell.n * subshell.n) / std::max(zeff, 1e-6);

    float max_value = 0.0f;
    for (int i = 0; i < kSamples; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(kSamples - 1);
        const double r = max_r * t;
        const double radial = sbox::hydrogen::radial_wavefunction(subshell.n, subshell.l, zeff, r);
        const double probability = r * r * radial * radial;
        values[static_cast<std::size_t>(i)] = static_cast<float>(probability);
        max_value = std::max(max_value, values[static_cast<std::size_t>(i)]);
    }

    if (max_value <= 0.0f) {
        max_value = 1.0f;
    }

    ImGui::PlotLines("##radial_plot",
                     values.data(),
                     kSamples,
                     0,
                     nullptr,
                     0.0f,
                     max_value * 1.05f,
                     ImVec2(0.0f, 120.0f));
    ImGui::TextUnformatted("r (a₀)");

    ImGui::End();
}

}  // namespace sbox::ui
