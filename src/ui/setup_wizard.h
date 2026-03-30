#pragma once

#include "backend/python_env.h"

#include <future>
#include <utility>
#include <string>

namespace sbox::ui {

struct SetupWizardState {
    bool show_wizard = false;
    bool detection_done = false;
    bool detection_running = false;

    sbox::backend::PythonInfo detected_info;

    char python_path_input[512] = "";

    bool installing = false;
    std::string install_log;
    bool install_complete = false;
    bool install_success = false;

    enum class Page {
        Welcome,
        Detecting,
        NotFound,
        Found,
        ManualPath,
        Installing,
        Done,
    };
    Page current_page = Page::Welcome;

    std::future<sbox::backend::PythonInfo> detection_future;
    std::future<std::pair<int, std::string>> install_future;
};

bool draw_setup_wizard(SetupWizardState& wizard_state, sbox::backend::PythonEnvironment& env);

}  // namespace sbox::ui
