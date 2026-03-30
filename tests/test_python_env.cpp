#define private public
#include "backend/python_env.h"
#undef private

#include <gtest/gtest.h>

TEST(PythonEnvTest, RunCaptureSucceedsForEcho) {
    std::string out;
    const int status = sbox::backend::PythonEnvironment::run_capture("echo hello", out);
    EXPECT_EQ(status, 0);
    EXPECT_NE(out.find("hello"), std::string::npos);
}

TEST(PythonEnvTest, RunCaptureReportsFailure) {
    std::string out;
    const int status = sbox::backend::PythonEnvironment::run_capture("false", out);
    EXPECT_NE(status, 0);
}

TEST(PythonEnvTest, DetectFindsWorkingPython) {
    sbox::backend::PythonEnvironment env;
    env.detect();
    EXPECT_TRUE(env.info().valid);
    EXPECT_FALSE(env.info().version.empty());
    EXPECT_FALSE(env.info().python_path.empty());
}

TEST(PythonEnvTest, CheckPackagesDoesNotCrash) {
    sbox::backend::PythonEnvironment env;
    env.detect();
    ASSERT_TRUE(env.info().valid);
    env.check_packages();
    EXPECT_TRUE(env.info().valid);
}
