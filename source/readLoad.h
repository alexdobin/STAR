#ifndef READLOAD_DEF
#define READLOAD_DEF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "SequenceFuns.h"

int readLoad(istream& readInStream, Parameters &P, uint iMate, uint& Lread, uint& readLengthPairOriginal, char* readName, char* Seq, char* SeqNum, char* Qual, char* QualNum, uint &clip3pNtotal, uint &clip5pNtotal, uint &clip3pAdapterN, uint &iReadAll, uint &readFilesIndex, char &readFilter, string &readNameExtra);

#endif
