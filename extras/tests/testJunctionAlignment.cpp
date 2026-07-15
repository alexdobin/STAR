#include "OutSJ.h"

#include <cstdint>
#include <type_traits>
#include <utility>

using JunctionStorage=decltype(std::declval<OutSJ>().dataVec);

static_assert(std::is_same<JunctionStorage::value_type, uint64>::value,
              "junction storage must provide uint64 alignment");
static_assert(Junction::startP % alignof(uint) == 0, "unaligned junction start");
static_assert(Junction::gapP % alignof(uint32) == 0, "unaligned junction gap");
static_assert(Junction::countUniqueP % alignof(uint32) == 0, "unaligned unique count");
static_assert(Junction::countMultipleP % alignof(uint32) == 0, "unaligned multiple count");
static_assert(Junction::overhangLeftP % alignof(uint16) == 0, "unaligned left overhang");
static_assert(Junction::overhangRightP % alignof(uint16) == 0, "unaligned right overhang");
static_assert(Junction::dataSize % alignof(uint64) == 0, "unaligned junction records");

int main()
{
    return 0;
}
