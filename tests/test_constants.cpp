#include "core/constants.h"

#include <gtest/gtest.h>

TEST(ConstantsTest, BohrRadiusAngstromMatchesExpectedValue) {
    EXPECT_DOUBLE_EQ(sbox::constants::kBohrRadiusAngstrom, 0.529177);
}

TEST(ConstantsTest, HartreeEnergyEvMatchesExpectedValue) {
    EXPECT_DOUBLE_EQ(sbox::constants::kHartreeEnergyEv, 27.2114);
}
