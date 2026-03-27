#pragma once

#include "core/molecular_system.h"

#include <json.hpp>

#include <string>

namespace sbox::io {

void save_project(const std::string& filepath,
                  const sbox::chem::MolecularSystem& mol,
                  const nlohmann::json& extra_state = {});

sbox::chem::MolecularSystem load_project(const std::string& filepath, nlohmann::json* out_state = nullptr);

}  // namespace sbox::io
