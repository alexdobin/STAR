#include <gtest/gtest.h>

#include "uint128_t.h"

TEST(BitShift, left){
    // operator<<
    uint128_t val(0x1);
    uint64_t exp_val = 1;
    for(uint8_t i = 0; i < 64; i++){
        EXPECT_EQ(val << i, exp_val << i);
    }

    uint128_t zero(0);
    for(uint8_t i = 0; i < 64; i++){
        EXPECT_EQ(zero << i, 0);
    }

    // operator<<=
    for(uint8_t i = 0; i < 63; i++){ // 1 is already a bit
        EXPECT_EQ(val  <<= 1, exp_val <<= 1);
    }

    for(uint8_t i = 0; i < 63; i++){
        EXPECT_EQ(zero <<= 1, 0);
    }
}

TEST(External, shift_left){
    bool     t   = true;
    bool     f   = false;
    uint8_t  u8  = 0xffULL;
    uint16_t u16 = 0xffffULL;
    uint32_t u32 = 0xffffffffULL;
    uint64_t u64 = 0xffffffffffffffffULL;

    const uint128_t zero(0);
    const uint128_t one(1);

    EXPECT_EQ(t   << zero, t);
    EXPECT_EQ(f   << zero, f);
    EXPECT_EQ(u8  << zero, u8);
    EXPECT_EQ(u16 << zero, u16);
    EXPECT_EQ(u32 << zero, u32);
    EXPECT_EQ(u64 << zero, u64);

    EXPECT_EQ(t   <<= zero, t);
    EXPECT_EQ(f   <<= zero, f);
    EXPECT_EQ(u8  <<= zero, u8);
    EXPECT_EQ(u16 <<= zero, u16);
    EXPECT_EQ(u32 <<= zero, u32);
    EXPECT_EQ(u64 <<= zero, u64);

    EXPECT_EQ(t   << one, uint128_t(t)   << 1);
    EXPECT_EQ(f   << one, uint128_t(f)   << 1);
    EXPECT_EQ(u8  << one, uint128_t(u8)  << 1);
    EXPECT_EQ(u16 << one, uint128_t(u16) << 1);
    EXPECT_EQ(u32 << one, uint128_t(u32) << 1);
    EXPECT_EQ(u64 << one, uint128_t(u64) << 1);

    EXPECT_EQ(t   <<= one, true);
    EXPECT_EQ(f   <<= one, false);
    EXPECT_EQ(u8  <<= one, (uint8_t)  0xfeULL);
    EXPECT_EQ(u16 <<= one, (uint16_t) 0xfffeULL);
    EXPECT_EQ(u32 <<= one, (uint32_t) 0xfffffffeULL);
    EXPECT_EQ(u64 <<= one, (uint64_t) 0xfffffffffffffffeULL);

    EXPECT_EQ(u8  << uint128_t(7),  uint128_t(0x7f00ULL));
    EXPECT_EQ(u16 << uint128_t(15), uint128_t(0x7fff0000ULL));
    EXPECT_EQ(u32 << uint128_t(31), uint128_t(0x7fffffff00000000ULL));
    EXPECT_EQ(u64 << uint128_t(63), uint128_t(0x7fffffffffffffff, 0x0000000000000000ULL));

    EXPECT_EQ(u8  <<= uint128_t(7),  (uint8_t)  0);
    EXPECT_EQ(u16 <<= uint128_t(15), (uint16_t) 0);
    EXPECT_EQ(u32 <<= uint128_t(31), (uint32_t) 0);
    EXPECT_EQ(u64 <<= uint128_t(63), (uint64_t) 0);
}
