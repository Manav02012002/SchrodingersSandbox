#include "core/crash_handler.h"

#include "core/logging.h"
#include "core/settings.h"
#include "version.h"

#include <csignal>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#if defined(_WIN32)
#include <windows.h>
#else
#include <cxxabi.h>
#include <execinfo.h>
#include <sys/utsname.h>
#include <unistd.h>
#endif

namespace sbox {

namespace {

std::string signal_name(int sig) {
    switch (sig) {
    case SIGSEGV: return "SIGSEGV";
    case SIGABRT: return "SIGABRT";
    case SIGFPE: return "SIGFPE";
#ifdef SIGBUS
    case SIGBUS: return "SIGBUS";
#endif
    case SIGILL: return "SIGILL";
    default: return "UNKNOWN";
    }
}

std::string signal_description(int sig) {
    switch (sig) {
    case SIGSEGV: return "Segmentation fault";
    case SIGABRT: return "Abort signal";
    case SIGFPE: return "Floating point exception";
#ifdef SIGBUS
    case SIGBUS: return "Bus error";
#endif
    case SIGILL: return "Illegal instruction";
    default: return "Unknown signal";
    }
}

std::string timestamp_compact() {
    const std::time_t now = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &now);
#else
    localtime_r(&now, &tm);
#endif
    char buffer[32];
    std::strftime(buffer, sizeof(buffer), "%Y%m%d_%H%M%S", &tm);
    return buffer;
}

std::string timestamp_readable() {
    const std::time_t now = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &now);
#else
    localtime_r(&now, &tm);
#endif
    char buffer[64];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S", &tm);
    return buffer;
}

std::string os_description() {
#if defined(_WIN32)
    return "Windows (stack traces not yet implemented)";
#else
    struct utsname info {};
    if (uname(&info) == 0) {
        return std::string(info.sysname) + " " + info.release + " " + info.version + " " + info.machine;
    }
    return "Unknown OS";
#endif
}

std::string get_backtrace(int skip = 2) {
#if defined(_WIN32)
    return "Stack trace unavailable on Windows yet (TODO)\n";
#else
    void* callstack[128];
    const int frames = backtrace(callstack, 128);
    char** symbols = backtrace_symbols(callstack, frames);

    std::string result;
    for (int i = skip; i < frames; ++i) {
        std::string line = symbols != nullptr ? std::string(symbols[i]) : std::string("<unknown>");
        const std::size_t start = line.find('(');
        const std::size_t end = line.find('+', start);
        if (start != std::string::npos && end != std::string::npos && end > start + 1) {
            const std::string mangled = line.substr(start + 1, end - start - 1);
            int status = 0;
            char* demangled = abi::__cxa_demangle(mangled.c_str(), nullptr, nullptr, &status);
            if (status == 0 && demangled != nullptr) {
                line = line.substr(0, start + 1) + demangled + line.substr(end);
                std::free(demangled);
            }
        }
        result += "  " + std::to_string(i - skip) + ": " + line + "\n";
    }
    std::free(symbols);
    return result;
#endif
}

void write_crash_report(int sig) {
    const std::filesystem::path path = std::filesystem::path(get_app_data_dir()) / ("crash_" + timestamp_compact() + ".log");
    std::ofstream out(path);
    if (!out) {
        std::cerr << "Failed to open crash log: " << path << '\n';
        return;
    }

    out << "Schrodinger's Sandbox crash report\n";
    out << "Timestamp: " << timestamp_readable() << "\n";
    out << "Version: " << VERSION << "\n";
    out << "Signal: " << signal_name(sig) << " (" << signal_description(sig) << ")\n";
    out << "OS: " << os_description() << "\n\n";
    out << "Stack trace:\n" << get_backtrace() << "\n";
    out << "Recent log entries:\n" << Logger::instance().get_recent_entries(50) << "\n";
    out.flush();

    std::cerr << "Crash log written to: " << path << '\n';
}

void crash_signal_handler(int sig) {
    write_crash_report(sig);
    std::signal(sig, SIG_DFL);
    std::raise(sig);
}

}  // namespace

void install_crash_handler() {
    std::signal(SIGSEGV, crash_signal_handler);
    std::signal(SIGABRT, crash_signal_handler);
    std::signal(SIGFPE, crash_signal_handler);
#ifdef SIGBUS
    std::signal(SIGBUS, crash_signal_handler);
#endif
    std::signal(SIGILL, crash_signal_handler);
}

}  // namespace sbox
