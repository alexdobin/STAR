#ifndef CODE_bam_cat
#define CODE_bam_cat

#include <htslib/sam.h>

int bam_cat(int nfn, char * const *fn, const bam_hdr_t *h, const char* outbam);

#endif
