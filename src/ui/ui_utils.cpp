#include "ui/ui_utils.h"

#include <algorithm>
#include <array>
#include <string>

namespace sbox::ui {

const std::array<const char*, 7> kLLabels = {"s", "p", "d", "f", "g", "h", "i"};

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

}  // namespace sbox::ui
