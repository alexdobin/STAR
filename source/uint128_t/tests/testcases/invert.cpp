#include <gtest/gtest.h>

#include "uint128_t.h"

TEST(BitWise, invert){
    EXPECT_EQ(~uint128_t(0x0000000000000000ULL, 0x0000000000000000ULL), uint128_t(0xffffffffffffffffULL, 0xffffffffffffffffULL));
    EXPECT_EQ(~uint128_t(0x0000000000000000ULL, 0xffffffffffffffffULL), uint128_t(0xffffffffffffffffULL, 0x0000000000000000ULL));
    EXPECT_EQ(~uint128_t(0xffffffffffffffffULL, 0xffffffffffffffffULL), uint128_t(0x0000000000000000ULL, 0x0000000000000000ULL));
}