#ifndef CODE_insertSeqSA
#define CODE_insertSeqSA

#include "IncludeDefine.h"
#include "PackedArray.h"
#include "Parameters.h"
#include "Genome.h"

uint insertSeqSA(PackedArray & SA, PackedArray & SA1, PackedArray & SAi, char * G, char * G1, uint64 nG, uint64 nG1, uint64 nG2, Parameters & P, Genome &mapGen);

#endif
