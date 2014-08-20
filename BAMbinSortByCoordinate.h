#ifndef CODE_BAMbinSortByCoordinate
#define CODE_BAMbinSortByCoordinate
#include "IncludeDefine.h"
#include "Parameters.h"
#include SAMTOOLS_BGZF_H

void BAMbinSortByCoordinate(uint32 iBin, uint binN, uint binS, uint nThreads, string dirBAMsort, BGZF *bgzfBAM, Parameters *P);
       
#endif