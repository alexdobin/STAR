#include <gtest/gtest.h>

#include "uint128_t.h"

TEST(Accessor, bits){
    uint128_t value = 1;
    for(uint32_t i = 0; i < 127; i++){
        EXPECT_EQ(value.bits(), i + 1);                     // before shift
        value <<= 1;
    }

    EXPECT_EQ(uint128_t(0).bits(), 0);
}

TEST(Accessor, data){
    const uint128_t value(0xfedcba9876543210ULL, 0x0123456789abcdefULL);
    EXPECT_EQ(value.upper(), 0xfedcba9876543210ULL);
    EXPECT_EQ(value.lower(), 0x0123456789abcdefULL);
}
