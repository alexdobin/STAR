#ifndef H_SoloCommon
#define H_SoloCommon

typedef struct{
    uint64 cb; 
    uint32 umi;
} readInfoStruct;

typedef struct{
    uint32 tr; 
    uint8  type;
} trTypeStruct;

typedef uint32 uintUMI;
typedef uint64 uintCB;
typedef uint32 uintRead;

#define uintUMIbits 32
#define velocytoTypeGeneBits 4
#define velocytoTypeGeneBitShift 28
#define geneMultMark (((uint32)1)<<31)

#endif