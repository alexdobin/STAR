/*  vcf_sweep.h -- forward/reverse sweep API.

    Copyright (C) 2013 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef HTSLIB_VCF_SWEEP_H
#define HTSLIB_VCF_SWEEP_H

#include "hts.h"
#include "vcf.h"

typedef struct _bcf_sweep_t bcf_sweep_t;

bcf_sweep_t *bcf_sweep_init(const char *fname);
void bcf_sweep_destroy(bcf_sweep_t *sw);
bcf_hdr_t *bcf_sweep_hdr(bcf_sweep_t *sw);
bcf1_t *bcf_sweep_fwd(bcf_sweep_t *sw);
bcf1_t *bcf_sweep_bwd(bcf_sweep_t *sw);

#endif
