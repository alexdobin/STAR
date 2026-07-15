#include "PackedArray.h"

#include <iostream>
#include <vector>

int main()
{
    const uint bits=13;
    const uint count=257;
    const uint mask=(1LLU<<bits)-1;

    PackedArray packed;
    packed.defineBits(bits, count);

    std::vector<char> storage(packed.lengthByte+1, 0);
    packed.pointArray(storage.data()+1);

    for (uint ii=0; ii<count; ii++)
    {
        packed.writePacked(ii, (ii*7919+17)&mask);
    };

    for (uint ii=0; ii<count; ii++)
    {
        uint expected=(ii*7919+17)&mask;
        if (packed[ii]!=expected)
        {
            std::cerr << "packed value mismatch at " << ii << ": expected "
                      << expected << ", observed " << packed[ii] << '\n';
            return 1;
        };
    };

    return 0;
}
