#ifndef CODE_genomeSAindex
#define CODE_genomeSAindex
#include "PackedArray.h"
#include "Parameters.h"
#include "Genome.h"

void genomeSAindex(char * G, PackedArray & SA, Parameters & P, PackedArray & SAip, Genome &mapGen);
void genomeSAindexChunk(char * G, PackedArray & SA, Parameters & P, PackedArray & SAi, uint iSA1, uint iSA2, Genome &mapGen);
void funSAiFindNextIndex(char *G, PackedArray &SA, uint isaStep, uint & isa, uint & indFull, int & iL4, Genome &mapGen);

#endif
