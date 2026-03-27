#pragma once

#include "core/molecular_system.h"

#include <string>
#include <vector>

namespace sbox::io {

sbox::chem::MolecularSystem read_sdf(const std::string& filepath);
sbox::chem::MolecularSystem read_sdf_string(const std::string& content);
std::vector<sbox::chem::MolecularSystem> read_sdf_multi(const std::string& filepath);
void write_sdf(const std::string& filepath, const sbox::chem::MolecularSystem& mol);

}  // namespace sbox::io
