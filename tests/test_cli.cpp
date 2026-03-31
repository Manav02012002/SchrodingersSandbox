#include "cli.h"

#include <gtest/gtest.h>

#include <string>
#include <vector>

namespace {

std::vector<char*> make_argv(std::vector<std::string>& storage) {
    std::vector<char*> argv;
    argv.reserve(storage.size());
    for (std::string& item : storage) {
        argv.push_back(item.data());
    }
    return argv;
}

TEST(CLI, DefaultsWithNoArguments) {
    std::vector<std::string> args = {"schrodingers_sandbox"};
    std::vector<char*> argv = make_argv(args);
    const sbox::CLIOptions options = sbox::parse_cli(static_cast<int>(argv.size()), argv.data());
    EXPECT_FALSE(options.show_help);
    EXPECT_FALSE(options.show_version);
    EXPECT_TRUE(options.input_file.empty());
    EXPECT_EQ(options.screenshot_width, 1920);
    EXPECT_EQ(options.screenshot_height, 1080);
    EXPECT_EQ(options.orbital, -1);
}

TEST(CLI, VersionFlag) {
    std::vector<std::string> args = {"schrodingers_sandbox", "--version"};
    std::vector<char*> argv = make_argv(args);
    const sbox::CLIOptions options = sbox::parse_cli(static_cast<int>(argv.size()), argv.data());
    EXPECT_TRUE(options.show_version);
}

TEST(CLI, PositionalInputFile) {
    std::vector<std::string> args = {"schrodingers_sandbox", "molecule.xyz"};
    std::vector<char*> argv = make_argv(args);
    const sbox::CLIOptions options = sbox::parse_cli(static_cast<int>(argv.size()), argv.data());
    EXPECT_EQ(options.input_file, "molecule.xyz");
}

TEST(CLI, ScreenshotAndResolution) {
    std::vector<std::string> args = {"schrodingers_sandbox", "--screenshot", "out.png", "molecule.xyz", "--resolution", "1920x1080"};
    std::vector<char*> argv = make_argv(args);
    const sbox::CLIOptions options = sbox::parse_cli(static_cast<int>(argv.size()), argv.data());
    EXPECT_EQ(options.screenshot_path, "out.png");
    EXPECT_EQ(options.input_file, "molecule.xyz");
    EXPECT_EQ(options.screenshot_width, 1920);
    EXPECT_EQ(options.screenshot_height, 1080);
}

TEST(CLI, OrbitalHomo) {
    std::vector<std::string> args = {"schrodingers_sandbox", "--orbital", "homo"};
    std::vector<char*> argv = make_argv(args);
    const sbox::CLIOptions options = sbox::parse_cli(static_cast<int>(argv.size()), argv.data());
    EXPECT_EQ(options.orbital, -1);
}

TEST(CLI, OrbitalInteger) {
    std::vector<std::string> args = {"schrodingers_sandbox", "--orbital", "5"};
    std::vector<char*> argv = make_argv(args);
    const sbox::CLIOptions options = sbox::parse_cli(static_cast<int>(argv.size()), argv.data());
    EXPECT_EQ(options.orbital, 5);
}

TEST(CLI, LogLevelDebug) {
    std::vector<std::string> args = {"schrodingers_sandbox", "--log-level", "debug"};
    std::vector<char*> argv = make_argv(args);
    const sbox::CLIOptions options = sbox::parse_cli(static_cast<int>(argv.size()), argv.data());
    EXPECT_EQ(options.log_level, "debug");
}

TEST(CLI, UnknownOptionIgnored) {
    std::vector<std::string> args = {"schrodingers_sandbox", "--foo", "molecule.xyz"};
    std::vector<char*> argv = make_argv(args);
    EXPECT_NO_THROW({
        const sbox::CLIOptions options = sbox::parse_cli(static_cast<int>(argv.size()), argv.data());
        EXPECT_EQ(options.input_file, "molecule.xyz");
    });
}

}  // namespace
