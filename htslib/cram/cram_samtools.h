/*
Copyright (c) 2010-2013 Genome Research Ltd.
Author: James Bonfield <jkb@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef _CRAM_SAMTOOLS_H_
#define _CRAM_SAMTOOLS_H_

/* Samtools compatible API */
#define bam_blk_size(b)  ((b)->l_data)
#define bam_set_blk_size(b,v) ((b)->data_len = (v))

#define bam_ref(b)       (b)->core.tid
#define bam_pos(b)       (b)->core.pos
#define bam_mate_pos(b)  (b)->core.mpos
#define bam_mate_ref(b)  (b)->core.mtid
#define bam_ins_size(b)  (b)->core.isize
#define bam_seq_len(b)   (b)->core.l_qseq
#define bam_cigar_len(b) (b)->core.n_cigar
#define bam_flag(b)      (b)->core.flag
#define bam_bin(b)       (b)->core.bin
#define bam_map_qual(b)  (b)->core.qual
#define bam_name_len(b)  (b)->core.l_qname
#define bam_name(b)      bam_get_qname((b))
#define bam_qual(b)      bam_get_qual((b))
#define bam_seq(b)       bam_get_seq((b))
#define bam_cigar(b)     bam_get_cigar((b))
#define bam_aux(b)       bam_get_aux((b))

#define bam_dup(b)       bam_copy1(bam_init1(), (b))

#define bam_free(b)      bam_destroy1((b))

#define bam_reg2bin(beg,end) hts_reg2bin((beg),(end),14,5)

#include "htslib/sam.h"

enum cigar_op {
    BAM_CMATCH_=BAM_CMATCH,
    BAM_CINS_=BAM_CINS,
    BAM_CDEL_=BAM_CDEL,
    BAM_CREF_SKIP_=BAM_CREF_SKIP,
    BAM_CSOFT_CLIP_=BAM_CSOFT_CLIP,
    BAM_CHARD_CLIP_=BAM_CHARD_CLIP,
    BAM_CPAD_=BAM_CPAD,
    BAM_CBASE_MATCH=BAM_CEQUAL,
    BAM_CBASE_MISMATCH=BAM_CDIFF
};

typedef bam1_t bam_seq_t;

#include "cram/sam_header.h"

bam_hdr_t *cram_header_to_bam(SAM_hdr *h);
SAM_hdr *bam_header_to_cram(bam_hdr_t *h);

int bam_construct_seq(bam_seq_t **bp, size_t extra_len,
		      const char *qname, size_t qname_len,
		      int flag,
		      int rname,      // Ref ID
		      int pos,
		      int end,        // aligned start/end coords
		      int mapq,
		      uint32_t ncigar, const uint32_t *cigar,
		      int mrnm,       // Mate Ref ID
		      int mpos,
		      int isize,
		      int len,
		      const char *seq,
		      const char *qual);

#endif /* _CRAM_SAMTOOLS_H_ */
