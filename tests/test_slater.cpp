#include "core/elements.h"
#include "core/slater.h"

#include <gtest/gtest.h>

namespace {

int electrons_in_subshell(const sbox::slater::ElectronConfig& config, int n, int l) {
    for (const auto& subshell : config) {
        if (subshell.n == n && subshell.l == l) {
            return subshell.electrons;
        }
    }
    return 0;
}

}  // namespace

TEST(SlaterTest, EffectiveNuclearChargeExamples) {
    EXPECT_NEAR(sbox::slater::compute_zeff(1, sbox::elements::get_element(1).config, 1, 0), 1.00, 1e-2);
    EXPECT_NEAR(sbox::slater::compute_zeff(2, sbox::elements::get_element(2).config, 1, 0), 1.70, 1e-2);
    EXPECT_NEAR(sbox::slater::compute_zeff(3, sbox::elements::get_element(3).config, 2, 0), 1.30, 1e-2);
    EXPECT_NEAR(sbox::slater::compute_zeff(6, sbox::elements::get_element(6).config, 2, 1), 3.25, 1e-2);
    EXPECT_NEAR(sbox::slater::compute_zeff(11, sbox::elements::get_element(11).config, 3, 0), 2.20, 1e-2);
    EXPECT_NEAR(sbox::slater::compute_zeff(26, sbox::elements::get_element(26).config, 3, 2), 6.25, 1e-2);
    EXPECT_NEAR(sbox::slater::compute_zeff(26, sbox::elements::get_element(26).config, 4, 0), 3.75, 1e-2);
}

TEST(SlaterTest, ElementAndLayoutLookups) {
    EXPECT_STREQ(sbox::elements::get_element(1).symbol, "H");
    EXPECT_STREQ(sbox::elements::get_element(118).symbol, "Og");
    EXPECT_EQ(sbox::elements::get_element("Fe").Z, 26);

    const auto& cu = sbox::elements::get_element("Cu");
    EXPECT_EQ(electrons_in_subshell(cu.config, 3, 2), 10);

    const auto h_pos = sbox::elements::PT_LAYOUT[0];
    const auto he_pos = sbox::elements::PT_LAYOUT[1];
    const auto la_pos = sbox::elements::PT_LAYOUT[56];
    const auto ac_pos = sbox::elements::PT_LAYOUT[88];

    EXPECT_EQ(h_pos.row, 0);
    EXPECT_EQ(h_pos.col, 0);
    EXPECT_EQ(he_pos.row, 0);
    EXPECT_EQ(he_pos.col, 17);
    EXPECT_EQ(la_pos.row, 8);
    EXPECT_EQ(la_pos.col, 2);
    EXPECT_EQ(ac_pos.row, 9);
    EXPECT_EQ(ac_pos.col, 2);
}
