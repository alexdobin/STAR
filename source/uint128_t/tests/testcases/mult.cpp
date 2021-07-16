#include <gtest/gtest.h>

#include "uint128_t.h"

TEST(Arithmetic, multiply){
    uint128_t val(0xfedbca9876543210ULL);

    EXPECT_EQ(val * val, uint128_t(0xfdb8e2bacbfe7cefULL, 0x010e6cd7a44a4100ULL));

    const uint128_t zero = 0;
    EXPECT_EQ(val  * zero, zero);
    EXPECT_EQ(zero * val,  zero);

    const uint128_t one = 1;
    EXPECT_EQ(val * one, val);
    EXPECT_EQ(one * val, val);
}

TEST(External, multiply){
    bool     t   = true;
    bool     f   = false;
    uint8_t  u8  = 0xaaULL;
    uint16_t u16 = 0xaaaaULL;
    uint32_t u32 = 0xaaaaaaaaULL;
    uint64_t u64 = 0xaaaaaaaaaaaaaaaaULL;

    const uint128_t val(0xf0f0f0f0f0f0f0f0, 0xf0f0f0f0f0f0f0f0ULL);

    EXPECT_EQ(t   *  val, val);
    EXPECT_EQ(f   *  val, 0);
    EXPECT_EQ(u8  *  val, uint128_t(0xffffffffffffffff, 0xffffffffffffff60ULL));
    EXPECT_EQ(u16 *  val, uint128_t(0xffffffffffffffff, 0xffffffffffff5f60ULL));
    EXPECT_EQ(u32 *  val, uint128_t(0xffffffffffffffff, 0xffffffff5f5f5f60ULL));
    EXPECT_EQ(u64 *  val, uint128_t(0xffffffffffffffff, 0x5f5f5f5f5f5f5f60ULL));

    EXPECT_EQ(t   *= val, true);
    EXPECT_EQ(f   *= val, false);
    EXPECT_EQ(u8  *= val, (uint8_t)                0x60ULL);
    EXPECT_EQ(u16 *= val, (uint16_t)             0x5f60ULL);
    EXPECT_EQ(u32 *= val, (uint32_t)         0x5f5f5f60ULL);
    EXPECT_EQ(u64 *= val, (uint64_t) 0x5f5f5f5f5f5f5f60ULL);
}
