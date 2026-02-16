#pragma once

#include "core/elements.h"
#include "core/slater.h"

#include <algorithm>

namespace sbox::ui {

struct AppState {
    int selected_Z = 1;
    int selected_orbital_index = -1;  // -1 means "last orbital" (valence)
    int selected_m = 0;
    bool needs_update = true;

    int render_mode = 0;  // 0=volume, 1=isosurface, 2=phase isosurface
    float iso_value = 0.01f;
    float gamma = 0.4f;

    int current_n = 1;
    int current_l = 0;
    float current_Zeff = 1.0f;

    void update() {
        selected_Z = std::clamp(selected_Z, 1, 118);

        const auto& element = sbox::elements::get_element(selected_Z);
        const auto& config = element.config;

        if (config.empty()) {
            selected_orbital_index = -1;
            selected_m = 0;
            current_n = 1;
            current_l = 0;
            current_Zeff = 1.0f;
            needs_update = false;
            return;
        }

        if (selected_orbital_index < 0 || selected_orbital_index >= static_cast<int>(config.size())) {
            selected_orbital_index = static_cast<int>(config.size()) - 1;
        }

        const auto& subshell = config[static_cast<std::size_t>(selected_orbital_index)];
        const int n = subshell.n;
        const int l = subshell.l;
        const double zeff = sbox::slater::compute_zeff(selected_Z, config, n, l);

        selected_m = std::clamp(selected_m, -l, l);

        current_n = n;
        current_l = l;
        current_Zeff = static_cast<float>(zeff);
        needs_update = false;
    }
};

}  // namespace sbox::ui
