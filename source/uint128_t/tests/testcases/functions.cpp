#include <map>

#include <gtest/gtest.h>

#include "uint128_t.h"

static const std::map <uint32_t, std::string> tests = {
    std::make_pair(2,  "10000100000101011000010101101100"),
    std::make_pair(3,  "12201102210121112101"),
    std::make_pair(4,  "2010011120111230"),
    std::make_pair(5,  "14014244043144"),
    std::make_pair(6,  "1003520344444"),
    std::make_pair(7,  "105625466632"),
    std::make_pair(8,  "20405302554"),
    std::make_pair(9,  "5642717471"),
    std::make_pair(10, "2216002924"),
    std::make_pair(11, "a3796a883"),
    std::make_pair(12, "51a175124"),
    std::make_pair(13, "294145645"),
    std::make_pair(14, "170445352"),
    std::make_pair(15, "ce82d6d4"),
    std::make_pair(16, "8415856c"),
    // std::make_pair(256, "uint128_t"),
};

TEST(Function, str){
    // number of leading 0s
    const std::string::size_type leading = 5;

    // make sure all of the test strings create the ASCII version of the string
    const uint128_t original(2216002924);
    for(std::pair <uint32_t const, std::string>  t : tests){
        EXPECT_EQ(original.str(t.first), t.second);
    }

    // add leading zeros
    for(uint32_t base = 2; base <= 16; base++){
        EXPECT_EQ(original.str(base, tests.at(base).size() + leading), std::string(leading, '0') + tests.at(base));
    }
}

TEST(Function, export_bits){
    const uint64_t u64 = 0x0123456789abcdefULL;
    const uint128_t value = u64;

    EXPECT_EQ(value, u64);

    const std::vector<uint8_t> full = {
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef
    };

    std::vector<uint8_t> bits;
    value.export_bits(bits);

    EXPECT_EQ(bits, full);
}

TEST(External, ostream){
    const uint128_t value(0xfedcba9876543210ULL);

    // write out octal uint128_t
    std::stringstream oct; oct << std::oct << value;
    EXPECT_EQ(oct.str(), "1773345651416625031020");

    // write out decimal uint128_t
    std::stringstream dec; dec << std::dec << value;
    EXPECT_EQ(dec.str(), "18364758544493064720");

    // write out hexadecimal uint128_t
    std::stringstream hex; hex << std::hex << value;
    EXPECT_EQ(hex.str(), "fedcba9876543210");

    // zero
    std::stringstream zero; zero << uint128_t();
    EXPECT_EQ(zero.str(), "0");
}
