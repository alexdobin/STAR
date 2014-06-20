#ifndef BAM_SAMPLE_H
#define BAM_SAMPLE_H

#include "kstring.h"

typedef struct {
	int n, m;
	char **smpl;
	void *rg2smid, *sm2id;
} bam_sample_t;

bam_sample_t *bam_smpl_init(void);
int bam_smpl_add(bam_sample_t *sm, const char *abs, const char *txt);
int bam_smpl_rg2smid(const bam_sample_t *sm, const char *fn, const char *rg, kstring_t *str);
void bam_smpl_destroy(bam_sample_t *sm);

#endif
