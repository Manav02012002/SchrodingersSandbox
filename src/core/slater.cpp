#include "core/slater.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>

namespace sbox::slater {

namespace {

constexpr std::array<const char*, 14> kGroupOrder = {
    "1sp", "2sp", "3sp", "3d", "4sp", "4d", "4f", "5sp", "5d", "5f", "6sp", "6d", "6f", "7sp"};

std::string group_name(int n, int l) {
    if (l <= 1) {
        return std::to_string(n) + "sp";
    }
    if (l == 2) {
        return std::to_string(n) + "d";
    }
    return std::to_string(n) + "f";
}

int group_index(const std::string& name) {
    const auto it = std::find_if(kGroupOrder.begin(), kGroupOrder.end(),
                                 [&name](const char* entry) { return name == entry; });
    if (it == kGroupOrder.end()) {
        return static_cast<int>(kGroupOrder.size());
    }
    return static_cast<int>(std::distance(kGroupOrder.begin(), it));
}

}  // namespace

double compute_zeff(int z, const ElectronConfig& config, int target_n, int target_l) {
    if (z <= 0 || target_n <= 0 || target_l < 0) {
        return 1.0;
    }

    const std::string target_group = group_name(target_n, target_l);
    const int target_group_index = group_index(target_group);

    int target_electrons = 0;
    for (const auto& subshell : config) {
        if (subshell.n == target_n && subshell.l == target_l) {
            target_electrons += subshell.electrons;
        }
    }

    double sigma = 0.0;

    for (const auto& subshell : config) {
        const std::string group = group_name(subshell.n, subshell.l);
        const int electrons = subshell.electrons;
        if (electrons <= 0) {
            continue;
        }

        if (group == target_group) {
            int same_group_count = electrons;
            if (subshell.n == target_n && subshell.l == target_l && target_electrons > 0) {
                same_group_count = std::max(0, electrons - 1);
            }
            const double coeff = (target_group == "1sp") ? 0.30 : 0.35;
            sigma += coeff * static_cast<double>(same_group_count);
            continue;
        }

        if (target_l <= 1) {
            if (subshell.n == target_n - 1) {
                sigma += 0.85 * static_cast<double>(electrons);
            } else if (subshell.n < target_n - 1) {
                sigma += static_cast<double>(electrons);
            }
            continue;
        }

        const int current_group_index = group_index(group);
        if (current_group_index < target_group_index) {
            sigma += static_cast<double>(electrons);
        }
    }

    return std::max(static_cast<double>(z) - sigma, 1.0);
}

}  // namespace sbox::slater
