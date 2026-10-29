#ifndef CODE_bam_cat
#define CODE_bam_cat

#include "htslib/htslib/sam.h"

int bam_cat(int nfn, char * const *fn, const bam_hdr_t *h, const char* outbam);

// bam_cat in steps, so that the output can be written while later inputs are
// still being made: bam_cat(n, fn, h, out) is bam_cat_open(h, out), then
// bam_cat_append(fn[i]) for i = 0..n-1 in order, then bam_cat_close.
typedef struct bam_cat_t bam_cat_t;
bam_cat_t *bam_cat_open(const bam_hdr_t *h, const char *outbam);
int bam_cat_append(bam_cat_t *bc, const char *fn);
int bam_cat_close(bam_cat_t *bc);

#endif
