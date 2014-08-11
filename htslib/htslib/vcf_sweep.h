#ifndef __VCF_SWEEP_H__
#define __VCF_SWEEP_H__

#include "hts.h"
#include "vcf.h"

typedef struct _bcf_sweep_t bcf_sweep_t;

bcf_sweep_t *bcf_sweep_init(const char *fname);
void bcf_sweep_destroy(bcf_sweep_t *sw);
bcf_hdr_t *bcf_sweep_hdr(bcf_sweep_t *sw);
bcf1_t *bcf_sweep_fwd(bcf_sweep_t *sw);
bcf1_t *bcf_sweep_bwd(bcf_sweep_t *sw);

#endif
