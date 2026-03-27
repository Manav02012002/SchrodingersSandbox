#pragma once

#include "core/basis_set.h"

#include <string>

namespace sbox::molden {

struct ParseOptions {
    // If false, parser renormalizes contraction coefficients per shell.
    bool contraction_coefficients_include_shell_normalization = true;
};

sbox::basis::MOData parse_molden_file(const std::string& filepath);
sbox::basis::MOData parse_molden_file(const std::string& filepath, const ParseOptions& options);

}  // namespace sbox::molden
