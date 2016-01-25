#ifndef PACKEDARRAY_DEF
#define PACKEDARRAY_DEF

#include "IncludeDefine.h"

class PackedArray {
    private:
        uint bitRecMask, wordCompLength;
        bool arrayAllocated; //true if charArray was allocated
    public:
        uint wordLength, length, lengthByte;
        uint operator [] (uint ii);
        char* charArray;

    PackedArray();
    void defineBits (uint Nbits, uint lengthIn);
    void writePacked(uint jj, uint x);
    void allocateArray();
    void deallocateArray();
    void pointArray(char* pointerCharIn);
//     PackedArray(uint N);
};

inline uint PackedArray::operator [] (uint ii) {
   uint b=ii*wordLength;
   uint B=b/8;
   uint S=b%8;

   uint a1 = *((uint*) (charArray+B));
   a1 = ((a1>>S)<<wordCompLength)>>wordCompLength;
   return a1;
};

#endif
