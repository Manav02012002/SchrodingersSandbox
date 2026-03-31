#include "cli.h"

#include "version.h"

#include <cstdlib>
#include <iostream>
#include <string>

namespace sbox {

namespace {

bool consume_value(int& i, int argc, char* argv[], std::string& out) {
    if (i + 1 >= argc) {
        return false;
    }
    out = argv[++i];
    return true;
}

void parse_resolution(const std::string& value, int& width, int& height) {
    const std::size_t x_pos = value.find_first_of("xX");
    if (x_pos == std::string::npos) {
        return;
    }
    try {
        width = std::stoi(value.substr(0, x_pos));
        height = std::stoi(value.substr(x_pos + 1));
    } catch (...) {
    }
}

}  // namespace

CLIOptions parse_cli(int argc, char* argv[]) {
    CLIOptions options;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i] != nullptr ? argv[i] : "";
        if (arg == "-h" || arg == "--help") {
            options.show_help = true;
        } else if (arg == "-v" || arg == "--version") {
            options.show_version = true;
        } else if (arg == "--headless") {
            options.headless = true;
        } else if (arg == "--compute") {
            options.compute = true;
        } else if (arg == "--optimize") {
            options.optimize = true;
        } else if (arg == "--screenshot") {
            std::string value;
            if (consume_value(i, argc, argv, value)) {
                options.screenshot_path = value;
            }
        } else if (arg == "--resolution") {
            std::string value;
            if (consume_value(i, argc, argv, value)) {
                parse_resolution(value, options.screenshot_width, options.screenshot_height);
            }
        } else if (arg == "--orbital") {
            std::string value;
            if (consume_value(i, argc, argv, value)) {
                if (value == "homo") {
                    options.orbital = -1;
                } else if (value == "lumo") {
                    options.orbital = -2;
                } else {
                    try {
                        options.orbital = std::stoi(value);
                    } catch (...) {
                    }
                }
            }
        } else if (arg == "--render-mode") {
            consume_value(i, argc, argv, options.render_mode);
        } else if (arg == "--iso") {
            std::string value;
            if (consume_value(i, argc, argv, value)) {
                try {
                    options.iso_value = std::stof(value);
                } catch (...) {
                }
            }
        } else if (arg == "--method") {
            consume_value(i, argc, argv, options.method);
        } else if (arg == "--basis") {
            consume_value(i, argc, argv, options.basis);
        } else if (arg == "--output") {
            consume_value(i, argc, argv, options.output_file);
        } else if (arg == "--log-level") {
            consume_value(i, argc, argv, options.log_level);
        } else if (arg == "--log-file") {
            consume_value(i, argc, argv, options.log_file);
        } else if (!arg.empty() && arg[0] == '-') {
            continue;
        } else if (options.input_file.empty()) {
            options.input_file = arg;
        }
    }

    return options;
}

void print_help() {
    std::cout
        << "Usage: schrodingers_sandbox [options] [file]\n\n"
        << "Options:\n"
        << "  -h, --help              Show this help message\n"
        << "  -v, --version           Show version and exit\n"
        << "  --headless              Run without GUI (for batch processing)\n"
        << "  --screenshot FILE       Open file, render, save screenshot, and exit\n"
        << "  --resolution WxH        Set screenshot resolution (default: 1920x1080)\n"
        << "  --orbital N             Select orbital N, or homo/lumo\n"
        << "  --render-mode MODE      volume, isosurface, phase\n"
        << "  --iso VALUE             Isosurface threshold (default: 0.01)\n"
        << "  --method METHOD         Computation method for --compute\n"
        << "  --basis BASIS           Basis set for --compute\n"
        << "  --compute               Run a single-point calculation on the input file\n"
        << "  --optimize              Run geometry optimization\n"
        << "  --output FILE           Output file for --compute results\n"
        << "  --log-level LEVEL       trace, debug, info, warn, error\n"
        << "  --log-file FILE         Write log to file\n\n"
        << "Positional:\n"
        << "  file                    Open this file on startup\n";
}

void print_version() {
    std::cout << "Schrodinger's Sandbox " << VERSION << '\n';
}

}  // namespace sbox
