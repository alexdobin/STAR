# include "PackedArray.h"

PackedArray::PackedArray() {
    charArray=NULL;
    arrayAllocated=false;
};

void PackedArray::defineBits(uint Nbits, uint lengthIn){
    wordLength=Nbits;
    wordCompLength=sizeof(uint)*8LLU-wordLength;
    bitRecMask=(~0LLU)>>wordCompLength;
    length=lengthIn;
    lengthByte=(length-1)*wordLength/8LLU+sizeof(uint);
//     lengthByte=((lengthByte+sizeof(uint)-1LLU)/sizeof(uint))*sizeof(uint);
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
    memset(charArray+lengthByte-sizeof(uint),0,sizeof(uint));//set the last 8 bytes to zero, since some of them may lnever be written
    arrayAllocated=true;
};

void PackedArray::deallocateArray() {
    if (arrayAllocated) {
        delete[] charArray;
        arrayAllocated=false;
    };
    charArray=NULL;
};
