#include "core/update_checker.h"

#include <gtest/gtest.h>

namespace {

TEST(UpdateChecker, VersionComparison) {
    EXPECT_TRUE(sbox::is_newer("1.2.0", "1.1.0"));
    EXPECT_FALSE(sbox::is_newer("1.1.0", "1.2.0"));
    EXPECT_FALSE(sbox::is_newer("1.1.0", "1.1.0"));
    EXPECT_TRUE(sbox::is_newer("2.0.0", "1.9.9"));
    EXPECT_TRUE(sbox::is_newer("1.0.1", "1.0.0"));
}

TEST(UpdateChecker, ParseGithubResponse) {
    const std::string json = R"({
        "tag_name": "v1.2.3",
        "html_url": "https://github.com/example/repo/releases/tag/v1.2.3",
        "body": "Release notes here",
        "published_at": "2026-03-30T12:00:00Z",
        "assets": [
            {
                "name": "SchrodingersSandbox-1.2.3-x86_64.AppImage",
                "browser_download_url": "https://example.com/appimage"
            }
        ]
    })";

    const sbox::UpdateInfo info = sbox::parse_github_release_response(json);
    EXPECT_EQ(info.latest_version, "1.2.3");
    EXPECT_EQ(info.release_url, "https://github.com/example/repo/releases/tag/v1.2.3");
    EXPECT_EQ(info.published_date, "2026-03-30T12:00:00Z");
    EXPECT_EQ(info.changelog, "Release notes here");
#if defined(__linux__)
    EXPECT_EQ(info.download_url, "https://example.com/appimage");
#endif
}

}  // namespace
