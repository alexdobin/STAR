#include <gtest/gtest.h>

#include "uint128_t.h"

TEST(Comparison, greater_than_or_equals){
    const uint128_t big  (0xffffffffffffffffULL, 0xffffffffffffffffULL);
    const uint128_t small(0x0000000000000000ULL, 0x0000000000000000ULL);

    EXPECT_EQ(small >= small,  true);
    EXPECT_EQ(small >= big,   false);

    EXPECT_EQ(big >= small,    true);
    EXPECT_EQ(big >= big,      true);
}

#define unsigned_compare_gte(Z)                                         \
do                                                                      \
{                                                                       \
    static_assert(std::is_unsigned <Z>::value, "Type must be signed");  \
                                                                        \
    const Z small = std::numeric_limits <Z>::min();                     \
    const Z big   = std::numeric_limits <Z>::max();                     \
                                                                        \
    const uint128_t int_small(small);                                   \
    const uint128_t int_big(big);                                       \
                                                                        \
    EXPECT_EQ(small >= int_small,  true);                               \
    EXPECT_EQ(small >= int_big,   false);                               \
                                                                        \
    EXPECT_EQ(big >= int_small,    true);                               \
    EXPECT_EQ(big >= int_big,      true);                               \
}                                                                       \
while (0)

#define signed_compare_gte(Z)                                           \
do                                                                      \
{                                                                       \
    static_assert(std::is_signed <Z>::value, "Type must be signed");    \
                                                                        \
    const Z small = 1;                                                  \
    const Z big = std::numeric_limits <Z>::max();                       \
                                                                        \
    const uint128_t int_small(small);                                   \
    const uint128_t int_big(big);                                       \
                                                                        \
    EXPECT_EQ(small >= int_small,  true);                               \
    EXPECT_EQ(small >= int_big,   false);                               \
                                                                        \
    EXPECT_EQ(big >= int_small,    true);                               \
    EXPECT_EQ(big >= int_big,      true);                               \
}                                                                       \
while (0)

TEST(External, greater_than_or_equals){
    unsigned_compare_gte(bool);
    unsigned_compare_gte(uint8_t);
    unsigned_compare_gte(uint16_t);
    unsigned_compare_gte(uint32_t);
    unsigned_compare_gte(uint64_t);
    signed_compare_gte(int8_t);
    signed_compare_gte(int16_t);
    signed_compare_gte(int32_t);
    signed_compare_gte(int64_t);
}
