#include <gtest/gtest.h>

#include "uint128_t.h"

TEST(Typecast, all){
    const uint128_t val(0xaaaaaaaaaaaaaaaaULL, 0xaaaaaaaaaaaaaaaaULL);

    EXPECT_EQ(static_cast <bool>     (uint128_t(true)),          true);
    EXPECT_EQ(static_cast <bool>     (uint128_t(false)),         false);
    EXPECT_EQ(static_cast <uint8_t>  (val),           (uint8_t)  0xaaULL);
    EXPECT_EQ(static_cast <uint16_t> (val),           (uint16_t) 0xaaaaULL);
    EXPECT_EQ(static_cast <uint32_t> (val),           (uint32_t) 0xaaaaaaaaULL);
    EXPECT_EQ(static_cast <uint64_t> (val),           (uint64_t) 0xaaaaaaaaaaaaaaaaULL);
}