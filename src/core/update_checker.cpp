#include "core/update_checker.h"

#include "backend/python_env.h"
#include "version.h"

#include <json.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <sstream>
#include <vector>

namespace sbox {

namespace {

using json = nlohmann::json;

std::string trim_version(std::string value) {
    if (!value.empty() && (value.front() == 'v' || value.front() == 'V')) {
        value.erase(value.begin());
    }
    return value;
}

std::vector<int> parse_semver(const std::string& version) {
    std::vector<int> parts;
    std::string cleaned = trim_version(version);
    const std::size_t dash = cleaned.find('-');
    if (dash != std::string::npos) {
        cleaned = cleaned.substr(0, dash);
    }
    std::stringstream ss(cleaned);
    std::string piece;
    while (std::getline(ss, piece, '.')) {
        std::string digits;
        for (char ch : piece) {
            if (ch >= '0' && ch <= '9') {
                digits.push_back(ch);
            } else {
                break;
            }
        }
        if (digits.empty()) {
            parts.push_back(0);
        } else {
            parts.push_back(std::stoi(digits));
        }
    }
    while (parts.size() < 3) {
        parts.push_back(0);
    }
    return parts;
}

std::string current_platform_asset_url(const json& assets) {
    if (!assets.is_array()) {
        return {};
    }
    const auto matches = [](const std::string& name, const std::vector<std::string>& needles) {
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        for (const std::string& needle : needles) {
            if (lower.find(needle) != std::string::npos) {
                return true;
            }
        }
        return false;
    };

#if defined(__APPLE__)
    const std::vector<std::string> preferred = {".dmg", "macos"};
#elif defined(_WIN32)
    const std::vector<std::string> preferred = {".msi", ".exe", "windows"};
#else
    const std::vector<std::string> preferred = {".appimage", "linux"};
#endif

    for (const auto& asset : assets) {
        const std::string name = asset.value("name", std::string{});
        if (matches(name, preferred)) {
            return asset.value("browser_download_url", std::string{});
        }
    }
    if (!assets.empty()) {
        return assets.front().value("browser_download_url", std::string{});
    }
    return {};
}

}  // namespace

bool is_newer(const std::string& latest, const std::string& current) {
    const std::vector<int> latest_parts = parse_semver(latest);
    const std::vector<int> current_parts = parse_semver(current);
    for (std::size_t i = 0; i < std::min(latest_parts.size(), current_parts.size()); ++i) {
        if (latest_parts[i] > current_parts[i]) {
            return true;
        }
        if (latest_parts[i] < current_parts[i]) {
            return false;
        }
    }
    return false;
}

UpdateInfo parse_github_release_response(const std::string& json_str) {
    UpdateInfo info;
    info.current_version = trim_version(VERSION);

    const json j = json::parse(json_str);
    info.latest_version = trim_version(j.value("tag_name", std::string{}));
    info.release_url = j.value("html_url", std::string{});
    info.published_date = j.value("published_at", std::string{});
    info.download_url = current_platform_asset_url(j.value("assets", json::array()));
    info.changelog = j.value("body", std::string{});
    if (info.changelog.size() > 500) {
        info.changelog = info.changelog.substr(0, 500) + "...";
    }
    info.update_available = !info.latest_version.empty() && is_newer(info.latest_version, info.current_version);
    return info;
}

UpdateChecker::UpdateChecker(const std::string& github_repo) : github_repo_(github_repo) {}

void UpdateChecker::check_async() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (checking_) {
        return;
    }
    cached_result_.reset();
    checking_ = true;
    future_ = std::async(std::launch::async, [this]() { return do_check(); });
}

std::optional<UpdateInfo> UpdateChecker::get_result() const {
    std::lock_guard<std::mutex> lock(mutex_);
    if (cached_result_.has_value()) {
        return cached_result_;
    }
    if (!future_.valid()) {
        return std::nullopt;
    }
    if (future_.wait_for(std::chrono::seconds(0)) != std::future_status::ready) {
        return std::nullopt;
    }
    cached_result_ = future_.get();
    checking_ = false;
    return cached_result_;
}

UpdateInfo UpdateChecker::check_sync(int timeout_seconds) {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (cached_result_.has_value()) {
            return *cached_result_;
        }
        if (!future_.valid()) {
            checking_ = true;
            future_ = std::async(std::launch::async, [this]() { return do_check(); });
        }
    }

    std::lock_guard<std::mutex> lock(mutex_);
    if (future_.wait_for(std::chrono::seconds(timeout_seconds)) != std::future_status::ready) {
        UpdateInfo info;
        info.current_version = trim_version(VERSION);
        info.check_failed = true;
        info.error_message = "Update check timed out";
        checking_ = false;
        cached_result_ = info;
        return info;
    }
    cached_result_ = future_.get();
    checking_ = false;
    return *cached_result_;
}

bool UpdateChecker::is_checking() const {
    return checking_;
}

UpdateInfo UpdateChecker::do_check() {
    UpdateInfo info;
    info.current_version = trim_version(VERSION);
    try {
        const std::string url = "https://api.github.com/repos/" + github_repo_ + "/releases/latest";
        std::string response;
        const int status = sbox::backend::PythonEnvironment::run_capture("curl -sL --max-time 10 " + url, response);
        if (status != 0 || response.empty()) {
            info.check_failed = true;
            info.error_message = "Failed to fetch release metadata";
            return info;
        }
        return parse_github_response(response);
    } catch (const std::exception& ex) {
        info.check_failed = true;
        info.error_message = ex.what();
        return info;
    }
}

UpdateInfo UpdateChecker::parse_github_response(const std::string& json_str) {
    try {
        return parse_github_release_response(json_str);
    } catch (const std::exception& ex) {
        UpdateInfo info;
        info.current_version = trim_version(VERSION);
        info.check_failed = true;
        info.error_message = ex.what();
        return info;
    }
}

}  // namespace sbox
