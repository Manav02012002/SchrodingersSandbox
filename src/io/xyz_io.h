#pragma once

#include "core/molecular_system.h"

#include <string>

namespace sbox::io {

sbox::chem::MolecularSystem read_xyz(const std::string& filepath);
sbox::chem::MolecularSystem read_xyz_string(const std::string& content);
void write_xyz(const std::string& filepath, const sbox::chem::MolecularSystem& mol);

}  // namespace sbox::io
