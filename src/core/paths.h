#pragma once

#include <string>

namespace sbox {

std::string get_data_dir();
std::string get_shader_path(const std::string& shader_name);
std::string get_script_path(const std::string& script_name);

}  // namespace sbox
