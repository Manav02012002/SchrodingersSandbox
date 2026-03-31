#pragma once

#include <atomic>
#include <future>
#include <mutex>
#include <optional>
#include <string>

namespace sbox {

struct UpdateInfo {
    bool update_available = false;
    std::string current_version;
    std::string latest_version;
    std::string release_url;
    std::string download_url;
    std::string changelog;
    std::string published_date;
    bool check_failed = false;
    std::string error_message;
};

bool is_newer(const std::string& latest, const std::string& current);
UpdateInfo parse_github_release_response(const std::string& json_str);

class UpdateChecker {
public:
    explicit UpdateChecker(const std::string& github_repo);

    void check_async();
    std::optional<UpdateInfo> get_result() const;
    UpdateInfo check_sync(int timeout_seconds = 10);
    bool is_checking() const;

private:
    std::string github_repo_;
    mutable std::future<UpdateInfo> future_;
    mutable std::mutex mutex_;
    mutable std::optional<UpdateInfo> cached_result_;
    mutable std::atomic<bool> checking_{false};

    UpdateInfo do_check();
    UpdateInfo parse_github_response(const std::string& json_str);
};

}  // namespace sbox
