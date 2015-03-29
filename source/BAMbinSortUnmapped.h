#ifndef CODE_BAMbinSortUnmapped
#define CODE_BAMbinSortUnmapped
#include "IncludeDefine.h"
#include "Parameters.h"
#include SAMTOOLS_BGZF_H

void BAMbinSortUnmapped(uint32 iBin, uint nThreads, string dirBAMsort, BGZF *bgzfBAM, Parameters *P);
       
#endif
