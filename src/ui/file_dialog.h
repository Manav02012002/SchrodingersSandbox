#pragma once

#include <string>

namespace sbox::ui {

std::string open_file_dialog(const char* title, const char* filters);
std::string save_file_dialog(const char* title, const char* filters, const char* default_name = nullptr);

}  // namespace sbox::ui
