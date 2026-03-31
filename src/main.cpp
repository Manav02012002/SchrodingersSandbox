#include "cli.h"
#include "renderer/screenshot.h"
#include "ui/app.h"
#include "core/crash_handler.h"
#include "core/logging.h"
#include "core/settings.h"
#include "version.h"

#include <exception>
#include <string>

namespace {

sbox::LogLevel parse_log_level(const std::string& level) {
    if (level == "trace") return sbox::LogLevel::Trace;
    if (level == "debug") return sbox::LogLevel::Debug;
    if (level == "warn" || level == "warning") return sbox::LogLevel::Warning;
    if (level == "error") return sbox::LogLevel::Error;
    return sbox::LogLevel::Info;
}

}  // namespace

int main(int argc, char* argv[]) {
    sbox::install_crash_handler();
    const sbox::CLIOptions options = sbox::parse_cli(argc, argv);

    if (options.show_help) {
        sbox::print_help();
        return 0;
    }
    if (options.show_version) {
        sbox::print_version();
        return 0;
    }

    auto& logger = sbox::Logger::instance();
    logger.set_level(parse_log_level(options.log_level));
    logger.set_file(options.log_file.empty() ? sbox::get_app_data_dir() + "schrodingers_sandbox.log" : options.log_file);
    SBOX_LOG_INFO("Schrodinger's Sandbox %s starting", sbox::VERSION);

    if (!options.screenshot_path.empty() && !options.input_file.empty()) {
        try {
            sbox::App app;
            app.load_file_by_extension(options.input_file);

            if (options.orbital >= 0) {
                app.state().selected_mo = options.orbital;
            } else if (options.orbital == -2 && app.state().homo_index >= 0) {
                app.state().selected_mo = app.state().homo_index + 1;
            }

            if (options.iso_value > 0.0f) {
                app.state().iso_value = options.iso_value;
            }

            if (options.render_mode == "volume") {
                app.state().render_mode = 0;
            } else if (options.render_mode == "isosurface") {
                app.state().render_mode = 1;
            } else if (options.render_mode == "phase") {
                app.state().render_mode = 2;
            }

            app.render_single_frame();
            if (!sbox::render::save_screenshot_highres(
                    options.screenshot_path,
                    options.screenshot_width,
                    options.screenshot_height,
                    [&app](unsigned int fbo, int w, int h) { app.render_to_fbo(fbo, w, h); })) {
                SBOX_LOG_FATAL("Screenshot failed: could not write %s", options.screenshot_path.c_str());
                logger.shutdown();
                return 1;
            }

            SBOX_LOG_INFO("Screenshot saved to %s", options.screenshot_path.c_str());
            logger.shutdown();
            return 0;
        } catch (const std::exception& ex) {
            SBOX_LOG_FATAL("Screenshot failed: %s", ex.what());
            logger.shutdown();
            return 1;
        }
    }

    try {
        sbox::App app;
        if (!options.input_file.empty()) {
            app.load_file_by_extension(options.input_file);
        }
        app.run();
    } catch (const std::exception& ex) {
        SBOX_LOG_FATAL("Fatal error: %s", ex.what());
        logger.shutdown();
        return 1;
    }

    SBOX_LOG_INFO("Application exiting normally");
    logger.shutdown();
    return 0;
}
