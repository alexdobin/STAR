#ifndef SUFFIXARRAYSFUNS_DEF
#define SUFFIXARRAYSFUNS_DEF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "PackedArray.h"

uint medianUint2(uint, uint);
uint compareSeqToGenome(char** s2, uint S, uint N, uint L, char* g, PackedArray& SA, uint iSA, bool dirR, bool& comparRes, Parameters* P); 
uint findMultRange(uint i3, uint L3, uint i1, uint L1, uint i1a, uint L1a, uint i1b, uint L1b, char** s, char* g, PackedArray& SA, bool dirR, uint S, Parameters* P);
uint maxMappableLength(char** s, uint S, uint N, char* g, PackedArray& SA, uint i1, uint i2, bool dirR, uint& L, uint* indStartEnd, Parameters* P);
void writePacked( char* a, uint jj, uint x);
uint readPacked(char* a, uint jj);
uint suffixArraySearch(char** s, uint S, uint N, char* g, PackedArray& SA, bool dirR, uint i1, uint i2, uint L, Parameters* P);
uint suffixArraySearch1(char** s2, uint S, uint N, char* G, uint64 gInsert, PackedArray& SA, bool dirR, uint i1, uint i2, uint L, Parameters* P);
int64 funCalcSAi(char* G, uint iL);
uint funCalcSAiFromSA(char* G, PackedArray& SA, uint iSA, int L, Parameters* P, int & iL4);
#endif
