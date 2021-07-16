#include <gtest/gtest.h>

#include "uint128_t.h"

TEST(Comparison, equals){
    EXPECT_EQ( (uint128_t(0xdeadbeefULL) == uint128_t(0xdeadbeefULL)), true);
    EXPECT_EQ(!(uint128_t(0xdeadbeefULL) == uint128_t(0xfee1baadULL)), true);
}

TEST(External, equals){
    const bool     t   = true;
    const bool     f   = false;
    const uint8_t  u8  = 0xaaULL;
    const uint16_t u16 = 0xaaaaULL;
    const uint32_t u32 = 0xaaaaaaaaULL;
    const uint64_t u64 = 0xaaaaaaaaaaaaaaaaULL;

    EXPECT_EQ(t,   uint128_t(t));
    EXPECT_EQ(f,   uint128_t(f));
    EXPECT_EQ(u8,  uint128_t(u8));
    EXPECT_EQ(u16, uint128_t(u16));
    EXPECT_EQ(u32, uint128_t(u32));
    EXPECT_EQ(u64, uint128_t(u64));
}
