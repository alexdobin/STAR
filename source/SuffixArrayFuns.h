#ifndef CODE_SuffixArrayFuns
#define CODE_SuffixArrayFuns

#include "IncludeDefine.h"
#include "Parameters.h"
#include "PackedArray.h"
#include "Genome.h"

#if defined(__GNUC__)
//always_inline: GCC 13 drops the calls to a copy left out of line, which it judges pure
__attribute__((always_inline)) inline void prefetchSuffixStart(Genome &mapGen, uint iSA, uint L)
{//start loading the genome byte where compareSeqToGenome(..., L, iSA, ...) begins its comparison; no effect on results
    uint SAstr=mapGen.SA[iSA];
    uint SApos=SAstr & mapGen.GstrandMask;
    uint gFwd=SApos+L, gRev=mapGen.nGenome-1-SApos-L;
    uint revMask=0-(SAstr>>mapGen.GstrandBit); //no branch: the strand is a coin flip, and a mispredict would discard the compare that follows
    __builtin_prefetch(mapGen.G + (gFwd ^ ((gFwd^gRev) & revMask)));
};

__attribute__((always_inline)) inline void prefetchPackedWords(PackedArray &A, uint i1, uint i2)
{//start loading the bytes that A[i1] to A[i2] read (i1<=i2); no effect on results
    __builtin_prefetch(A.charArray + i1*A.wordLength/8);
    __builtin_prefetch(A.charArray + i2*A.wordLength/8 + 7);
};
#else
inline void prefetchSuffixStart(Genome &, uint, uint) {};
inline void prefetchPackedWords(PackedArray &, uint, uint) {};
#endif

uint medianUint2(uint, uint);
uint compareSeqToGenome(Genome &mapGen, char** s2, uint S, uint N, uint L, uint iSA, bool dirR, bool& comparRes);
uint findMultRange(Genome &mapGen, uint i3, uint L3, uint i1, uint L1, uint i1a, uint L1a, uint i1b, uint L1b, char** s, bool dirR, uint S);
uint maxMappableLength(Genome &mapGen, char** s, uint S, uint N, uint i1, uint i2, bool dirR, uint& L, uint* indStartEnd);
void writePacked(Genome &mapGen, char* a, uint jj, uint x);
uint readPacked(Genome &mapGen, char* a, uint jj);
uint suffixArraySearch1(Genome &mapGen, char** s2, uint S, uint N, uint64 gInsert, bool dirR, uint i1, uint i2, uint L);
uint compareSeqToGenome1(Genome &mapGen, char** s2, uint S, uint N, uint L, uint iSA, bool dirR, uint64 gInsert, int & compRes);
int64 funCalcSAi(char *G, uint iL);
uint funCalcSAiFromSA(char* gSeq, PackedArray& gSA, Genome &mapGen, uint iSA, int L, int & iL4);
#endif
