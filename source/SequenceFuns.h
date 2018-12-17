/**  general basic functions to operate on sequences, no classes
 */

#ifndef SEQUENCEFUNS_DEF
#define SEQUENCEFUNS_DEF

#include "IncludeDefine.h"

void complementSeqNumbers(char*, char*, uint);
void convertNucleotidesToNumbers(const char* R0, char* R1, const uint Lread);
uint convertNucleotidesToNumbersRemoveControls(const char* R0, char* R1, const uint Lread);
void revComplementNucleotides(char* ReadsIn, char* ReadsOut, uint Lread); //complement the numeric sequences
void revComplementNucleotides(string &seq);
char nuclToNumBAM(char cc);
void nuclPackBAM(char* ReadsIn, char* ReadsOut, uint Lread);
char convertNt01234(const char R0);//transform sequence  from ACGT into 0-1-2-3-4 code    

uint chrFind(uint, uint, uint*); // find chromosome from global locus
uint localSearch(const char*, uint, const char*, uint, double); //local search to clip adapter
uint localSearchNisMM(const char *x, uint nx, const char *y, uint ny, double pMM);

uint qualitySplit(char*, char*, uint, char, uint, uint, uint**);

int32 convertNuclStrToInt32(const string S, uint32 &intOut);
string convertNuclInt32toString(const uint32 nuclNum, const uint32 L);

#endif
