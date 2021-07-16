#include <gtest/gtest.h>

#include "uint128_t.h"

TEST(Arithmetic, modulo){
    // has remainder
    const uint128_t val    (0xffffffffffffffffULL, 0xffffffffffffffffULL);
    const uint128_t val_mod(0xfedcba9876543210ULL);

    EXPECT_EQ(val % val_mod, uint128_t(0x7f598f328cc265bfULL));

    // no remainder
    const uint128_t val_0  (0xfedcba9876543210, 0);
    EXPECT_EQ(val_0 % val_mod, 0);

    // mod 0
    EXPECT_THROW(uint128_t(1) % uint128_t(0), std::domain_error);
}

TEST(External, modulo){
    bool     t   = true;
    bool     f   = false;
    uint8_t  u8  = 0xaaULL;
    uint16_t u16 = 0xaaaaULL;
    uint32_t u32 = 0xaaaaaaaaULL;
    uint64_t u64 = 0xaaaaaaaaaaaaaaaaULL;

    const uint128_t val(0xd03ULL); // prime

    EXPECT_EQ(t   %  val, true);
    EXPECT_EQ(f   %  val, false);
    EXPECT_EQ(u8  %  val, uint128_t(0xaaULL));
    EXPECT_EQ(u16 %  val, uint128_t(0x183ULL));
    EXPECT_EQ(u32 %  val, uint128_t(0x249ULL));
    EXPECT_EQ(u64 %  val, uint128_t(0xc7fULL));

    EXPECT_EQ(t   %= val, true);
    EXPECT_EQ(f   %= val, false);
    EXPECT_EQ(u8  %= val, (uint8_t)  0xaaULL);
    EXPECT_EQ(u16 %= val, (uint16_t) 0x183ULL);
    EXPECT_EQ(u32 %= val, (uint32_t) 0x249ULL);
    EXPECT_EQ(u64 %= val, (uint64_t) 0xc7fULL);
}