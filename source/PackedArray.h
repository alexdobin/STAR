#ifndef PACKEDARRAY_DEF
#define PACKEDARRAY_DEF

#include "IncludeDefine.h"

class PackedArray {
    private:
        uint bitRecMask, wordCompLength;    
    public:
        uint wordLength, length, lengthByte;
        uint operator [] (uint ii);
        char* charArray;
        
    void defineBits (uint Nbits, uint lengthIn);
    void writePacked(uint jj, uint x);
    void allocateArray();
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
