#pragma once

#include "core/slater.h"

#include <array>
#include <string>

namespace sbox::ui {

extern const std::array<const char*, 7> kLLabels;

std::string superscript_number(int value);
std::string config_to_string(const sbox::slater::ElectronConfig& config);

}  // namespace sbox::ui
