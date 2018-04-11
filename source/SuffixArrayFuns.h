#ifndef CODE_SuffixArrayFuns
#define CODE_SuffixArrayFuns

#include "IncludeDefine.h"
#include "Parameters.h"
#include "PackedArray.h"
#include "Genome.h"

uint medianUint2(uint, uint);
uint compareSeqToGenome(Genome &mapGen, char** s2, uint S, uint N, uint L, uint iSA, bool dirR, bool& comparRes);
uint findMultRange(Genome &mapGen, uint i3, uint L3, uint i1, uint L1, uint i1a, uint L1a, uint i1b, uint L1b, char** s, bool dirR, uint S);
uint maxMappableLength(Genome &mapGen, char** s, uint S, uint N, uint i1, uint i2, bool dirR, uint& L, uint* indStartEnd);
void writePacked(Genome &mapGen, char* a, uint jj, uint x);
uint readPacked(Genome &mapGen, char* a, uint jj);
uint suffixArraySearch1(Genome &mapGen, char** s2, uint S, uint N, uint64 gInsert, bool dirR, uint i1, uint i2, uint L);
int64 funCalcSAi(char *G, uint iL);
uint funCalcSAiFromSA(char* gSeq, PackedArray& gSA, Genome &mapGen, uint iSA, int L, int & iL4);
#endif
