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

void PackedArray::pointArray(char* pointerCharIn) {
    charArray=pointerCharIn;
};

void adviseHugePages(char* array, uint nBytes) {
#ifdef MADV_HUGEPAGE //Linux only: elsewhere the array keeps the default pages
    const uintptr_t hugeBytes=2LLU<<20;
    uintptr_t start=((uintptr_t) array + hugeBytes - 1) & ~(hugeBytes-1);
    uintptr_t end=((uintptr_t) array + nBytes) & ~(hugeBytes-1);
    if (end>start)
        madvise((void*) start, end-start, MADV_HUGEPAGE);//a hint: on failure the pages stay 4 KB, the contents are unaffected
#endif
};

void PackedArray::allocateArray() {
    charArray=new char[lengthByte];
    adviseHugePages(charArray, lengthByte);
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
