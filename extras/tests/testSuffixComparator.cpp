#include "funCompareUintAndSuffixesMemcmp.h"

#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

namespace {

struct Entry {
    uint64_t insertionPoint;
    uint64_t suffix;
};

bool expectSign(const Entry &left, const Entry &right, int expectedSign, const char *label)
{
    int result=funCompareUintAndSuffixesMemcmp(&left, &right);
    int resultSign=(result>0)-(result<0);
    if (resultSign!=expectedSign)
    {
        std::cerr << label << ": expected sign " << expectedSign << ", observed " << resultSign << '\n';
        return false;
    };
    return true;
}

}

int main()
{
    std::vector<char> sequence={0,1,3,5, 0,1,2,5, 0,1,3,5, 2,2,2,5};
    g_funCompareUintAndSuffixesMemcmp_G=sequence.data();
    g_funCompareUintAndSuffixesMemcmp_N=sequence.size()-1;
    g_funCompareUintAndSuffixesMemcmp_L=std::numeric_limits<uint64_t>::max();

    Entry a={7,0};
    Entry b={7,4};
    Entry sameSuffix={7,8};
    Entry earlierBucket={6,12};

    bool ok=true;
    ok &= expectSign(a, b, 1, "suffix order");
    ok &= expectSign(b, a, -1, "reverse suffix order");
    ok &= expectSign(a, sameSuffix, -1, "identical suffix tie break");
    ok &= expectSign(a, a, 0, "self comparison");
    ok &= expectSign(earlierBucket, a, -1, "insertion point order");

    g_funCompareUintAndSuffixesMemcmp_L=2;
    ok &= expectSign(a, b, -1, "finite prefix tie break");

    return ok ? 0 : 1;
}
