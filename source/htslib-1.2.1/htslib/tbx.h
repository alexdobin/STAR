/*  tbx.h -- tabix API functions.

    Copyright (C) 2009, 2012-2014 Genome Research Ltd.
    Copyright (C) 2010, 2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#ifndef HTSLIB_TBX_H
#define HTSLIB_TBX_H

#include "hts.h"

#define TBX_MAX_SHIFT 31

#define TBX_GENERIC 0
#define TBX_SAM     1
#define TBX_VCF     2
#define TBX_UCSC    0x10000

typedef struct {
    int32_t preset;
    int32_t sc, bc, ec; // seq col., beg col. and end col.
    int32_t meta_char, line_skip;
} tbx_conf_t;

typedef struct {
    tbx_conf_t conf;
    hts_idx_t *idx;
    void *dict;
} tbx_t;

extern tbx_conf_t tbx_conf_gff, tbx_conf_bed, tbx_conf_psltbl, tbx_conf_sam, tbx_conf_vcf;

#ifdef __cplusplus
extern "C" {
#endif

    #define tbx_itr_destroy(iter) hts_itr_destroy(iter)
    #define tbx_itr_queryi(tbx, tid, beg, end) hts_itr_query((tbx)->idx, (tid), (beg), (end), tbx_readrec)
    #define tbx_itr_querys(tbx, s) hts_itr_querys((tbx)->idx, (s), (hts_name2id_f)(tbx_name2id), (tbx), hts_itr_query, tbx_readrec)
    #define tbx_itr_next(htsfp, tbx, itr, r) hts_itr_next(hts_get_bgzfp(htsfp), (itr), (r), (tbx))
    #define tbx_bgzf_itr_next(bgzfp, tbx, itr, r) hts_itr_next((bgzfp), (itr), (r), (tbx))

    int tbx_name2id(tbx_t *tbx, const char *ss);

    /* Internal helper function used by tbx_itr_next() */
    BGZF *hts_get_bgzfp(htsFile *fp);
    int tbx_readrec(BGZF *fp, void *tbxv, void *sv, int *tid, int *beg, int *end);

    int tbx_index_build(const char *fn, int min_shift, const tbx_conf_t *conf);
    tbx_t *tbx_index_load(const char *fn);
    const char **tbx_seqnames(tbx_t *tbx, int *n);  // free the array but not the values
    void tbx_destroy(tbx_t *tbx);

#ifdef __cplusplus
}
#endif

#endif
