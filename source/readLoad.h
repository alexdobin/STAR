#ifndef READLOAD_DEF
#define READLOAD_DEF

#include "IncludeDefine.h"
#include "Parameters.h"
#include "SequenceFuns.h"

int readLoad(istream& readInStream, Parameters& P, uint32 iMate, uint& Lread, uint& LreadOriginal, \
		     char* readName, char* Seq, char* SeqNum, char* Qual, char* QualNum, \
			 uint32 &clip3pNtotal, uint32 &clip5pNtotal, \
			 uint &iReadAll, uint32 &readFilesIndex, char &readFilter, string &readNameExtra);

#endif
