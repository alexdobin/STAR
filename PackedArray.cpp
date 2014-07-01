# include "PackedArray.h"

void PackedArray::defineBits(uint Nbits, uint lengthIn){
    wordLength=Nbits;
    wordCompLength=sizeof(uint)*8LLU-wordLength;
    bitRecMask=(~0LLU)>>wordCompLength;
    length=lengthIn;
    lengthByte=length*Nbits/8LLU+1LLU;
};

void PackedArray::writePacked( uint jj, uint x) {
   uint b=jj*wordLength;
   uint B=b/8LLU;
   uint S=b%8LLU;

   x = x << S;
   uint* a1 = (uint*) (charArray+B);
   *a1 = ( (*a1) & ~(bitRecMask<<S) ) | x;
};

void PackedArray::pointArray(char* pointerCharIn) {
    charArray=pointerCharIn;
};

void PackedArray::allocateArray() {
    charArray=new char[lengthByte];
};
