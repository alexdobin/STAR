/**  general basic functions to operate on sequences, no classes
 *
 *
 * A longer description.
 *  
 * @see something
 */

#ifndef SEQUENCEFUNS_DEF
#define SEQUENCEFUNS_DEF

#include "IncludeDefine.h"

void complementSeqNumbers(char*, char*, uint);
void convertNucleotidesToNumbers(const char*, char*, uint);
void revComplementNucleotides(char* ReadsIn, char* ReadsOut, uint Lread); //complement the numeric sequences
char nuclToNumBAM(char cc);
void nuclPackBAM(char* ReadsIn, char* ReadsOut, uint Lread);
uint chrFind(uint, uint, uint*); // find chromosome from global locus
uint localSearch(const char*, uint, const char*, uint, double); //local search to clip adapter
uint qualitySplit(char*, char*, uint, char, uint, uint, uint**);
#endif
