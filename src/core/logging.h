#pragma once

#include <cstddef>
#include <deque>
#include <fstream>
#include <mutex>
#include <string>

namespace sbox {

enum class LogLevel {
    Trace,
    Debug,
    Info,
    Warning,
    Error,
    Fatal,
};

class Logger {
public:
    static Logger& instance();

    void set_level(LogLevel level);
    void set_file(const std::string& path);
    void set_console(bool enabled);

    void log(LogLevel level, const char* file, int line, const char* fmt, ...);
    std::string get_recent_entries(int count) const;
    void shutdown();

private:
    Logger();
    ~Logger();

    const char* level_string(LogLevel level) const;
    std::string timestamp() const;

    LogLevel level_ = LogLevel::Info;
    bool console_enabled_ = true;
    std::ofstream log_file_;
    mutable std::mutex mutex_;
    std::deque<std::string> recent_entries_;

    static constexpr std::size_t MAX_RECENT = 200;
};

}  // namespace sbox

#define SBOX_LOG_TRACE(...) sbox::Logger::instance().log(sbox::LogLevel::Trace, __FILE__, __LINE__, __VA_ARGS__)
#define SBOX_LOG_DEBUG(...) sbox::Logger::instance().log(sbox::LogLevel::Debug, __FILE__, __LINE__, __VA_ARGS__)
#define SBOX_LOG_INFO(...) sbox::Logger::instance().log(sbox::LogLevel::Info, __FILE__, __LINE__, __VA_ARGS__)
#define SBOX_LOG_WARN(...) sbox::Logger::instance().log(sbox::LogLevel::Warning, __FILE__, __LINE__, __VA_ARGS__)
#define SBOX_LOG_ERROR(...) sbox::Logger::instance().log(sbox::LogLevel::Error, __FILE__, __LINE__, __VA_ARGS__)
#define SBOX_LOG_FATAL(...) sbox::Logger::instance().log(sbox::LogLevel::Fatal, __FILE__, __LINE__, __VA_ARGS__)
