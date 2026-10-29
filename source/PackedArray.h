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

void adviseHugePages(char* array, uint nBytes);//ask for transparent huge pages on the 2 MB blocks inside a large array, before it is first touched

inline uint PackedArray::operator [] (uint ii) {
   uint b=ii*wordLength;
   uint B=b/8;
   uint S=b%8;

   uint a1 = *((uint*) (charArray+B));
   a1 = ((a1>>S)<<wordCompLength)>>wordCompLength;
   return a1;
};

inline void PackedArray::writePacked( uint jj, uint x) {
   uint b=jj*wordLength;
   uint B=b/8LLU;
   uint S=b%8LLU;

   x = x << S;
   uint* a1 = (uint*) (charArray+B);
   *a1 = ( (*a1) & ~(bitRecMask<<S) ) | x;
};

#endif
