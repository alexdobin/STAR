#ifndef CODE_BAMbinSortUnmapped
#define CODE_BAMbinSortUnmapped
#include "IncludeDefine.h"
#include "Parameters.h"

#include <htslib/bgzf.h>

void BAMbinSortUnmapped(uint32 iBin, uint nThreads, string dirBAMsort, BGZF *bgzfBAM, Parameters *P);
       
#endif
