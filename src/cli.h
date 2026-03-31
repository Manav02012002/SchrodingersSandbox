#pragma once

#include <string>

namespace sbox {

struct CLIOptions {
    std::string input_file;
    bool show_help = false;
    bool show_version = false;
    bool headless = false;
    std::string screenshot_path;
    int screenshot_width = 1920;
    int screenshot_height = 1080;
    int orbital = -1;
    std::string render_mode;
    float iso_value = -1.0f;
    std::string method;
    std::string basis;
    bool compute = false;
    bool optimize = false;
    std::string output_file;
    std::string log_level = "info";
    std::string log_file;
};

CLIOptions parse_cli(int argc, char* argv[]);
void print_help();
void print_version();

}  // namespace sbox
