#include "core/logging.h"

#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace sbox {

Logger& Logger::instance() {
    static Logger logger;
    return logger;
}

Logger::Logger() = default;

Logger::~Logger() {
    shutdown();
}

void Logger::set_level(LogLevel level) {
    std::lock_guard<std::mutex> lock(mutex_);
    level_ = level;
}

void Logger::set_file(const std::string& path) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!path.empty()) {
        std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    }
    if (log_file_.is_open()) {
        log_file_.flush();
        log_file_.close();
    }
    log_file_.open(path, std::ios::out | std::ios::app);
    if (!log_file_) {
        throw std::runtime_error("Failed to open log file: " + path);
    }
}

void Logger::set_console(bool enabled) {
    std::lock_guard<std::mutex> lock(mutex_);
    console_enabled_ = enabled;
}

void Logger::log(LogLevel level, const char* file, int line, const char* fmt, ...) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (static_cast<int>(level) < static_cast<int>(level_)) {
        return;
    }

    va_list args;
    va_start(args, fmt);
    va_list args_copy;
    va_copy(args_copy, args);
    const int required = std::vsnprintf(nullptr, 0, fmt, args_copy);
    va_end(args_copy);

    std::vector<char> buffer(static_cast<std::size_t>(std::max(required, 0)) + 1u, '\0');
    std::vsnprintf(buffer.data(), buffer.size(), fmt, args);
    va_end(args);

    const std::string filename = file != nullptr ? std::filesystem::path(file).filename().string() : std::string("unknown");
    const std::string entry = "[" + timestamp() + "] [" + level_string(level) + "] [" + filename + ":" + std::to_string(line) +
                              "] " + std::string(buffer.data());

    recent_entries_.push_back(entry);
    while (recent_entries_.size() > MAX_RECENT) {
        recent_entries_.pop_front();
    }

    if (console_enabled_) {
        std::cerr << entry << '\n';
    }
    if (log_file_.is_open()) {
        log_file_ << entry << '\n';
        log_file_.flush();
    }
}

std::string Logger::get_recent_entries(int count) const {
    std::lock_guard<std::mutex> lock(mutex_);
    if (count <= 0 || recent_entries_.empty()) {
        return {};
    }
    const std::size_t start = recent_entries_.size() > static_cast<std::size_t>(count)
                                  ? recent_entries_.size() - static_cast<std::size_t>(count)
                                  : 0u;
    std::string joined;
    for (std::size_t i = start; i < recent_entries_.size(); ++i) {
        if (!joined.empty()) {
            joined.push_back('\n');
        }
        joined += recent_entries_[i];
    }
    return joined;
}

void Logger::shutdown() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (log_file_.is_open()) {
        log_file_.flush();
        log_file_.close();
    }
}

const char* Logger::level_string(LogLevel level) const {
    switch (level) {
    case LogLevel::Trace: return "TRACE";
    case LogLevel::Debug: return "DEBUG";
    case LogLevel::Info: return "INFO ";
    case LogLevel::Warning: return "WARN ";
    case LogLevel::Error: return "ERROR";
    case LogLevel::Fatal: return "FATAL";
    }
    return "INFO ";
}

std::string Logger::timestamp() const {
    using namespace std::chrono;
    const auto now = system_clock::now();
    const auto secs = time_point_cast<seconds>(now);
    const auto millis = duration_cast<milliseconds>(now - secs).count();
    const std::time_t t = system_clock::to_time_t(now);

    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif

    char buffer[64];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S", &tm);
    char full[80];
    std::snprintf(full, sizeof(full), "%s.%03lld", buffer, static_cast<long long>(millis));
    return full;
}

}  // namespace sbox
