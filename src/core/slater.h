#pragma once

#include <vector>

namespace sbox::slater {

struct SubshellConfig {
    int n;
    int l;
    int electrons;
};

using ElectronConfig = std::vector<SubshellConfig>;

double compute_zeff(int z, const ElectronConfig& config, int target_n, int target_l);

}  // namespace sbox::slater
