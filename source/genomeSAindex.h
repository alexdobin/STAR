#ifndef CODE_genomeSAindex
#define CODE_genomeSAindex
#include "PackedArray.h"
#include "Parameters.h"

void genomeSAindex(char * G, PackedArray & SA, Parameters * P, PackedArray & SAip);
void genomeSAindexChunk(char * G, PackedArray & SA, Parameters * P, PackedArray & SAi, uint iSA1, uint iSA2);
void funSAiFindNextIndex(Parameters * P, char * G, PackedArray & SA, uint isaStep, uint & isa, uint & indFull, int & iL4);

#endif
