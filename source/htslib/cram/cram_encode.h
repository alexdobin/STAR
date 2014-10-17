/*
Copyright (c) 2012-2013 Genome Research Ltd.
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

/*! \file
 * Include cram.h instead.
 *
 * This is an internal part of the CRAM system and is automatically included
 * when you #include cram.h.
 *
 * Implements the encoding portion of CRAM I/O. Also see
 * cram_codecs.[ch] for the actual encoding functions themselves.
 */

#ifndef _CRAM_WRITE_H_
#define _CRAM_WRITE_H_

#ifdef __cplusplus
extern "C" {
#endif

/* ----------------------------------------------------------------------
 * CRAM sequence iterators.
 */

/*! Write iterator: put BAM format sequences into a CRAM file.
 *
 * We buffer up a containers worth of data at a time.
 *
 * FIXME: break this into smaller pieces.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_put_bam_seq(cram_fd *fd, bam_seq_t *b);


/* ----------------------------------------------------------------------
 * Internal functions
 */

/*! INTERNAL:
 * Encodes a compression header block into a generic cram_block structure.
 *
 * @return
 * Returns cram_block ptr on success;
 *         NULL on failure
 */
cram_block *cram_encode_compression_header(cram_fd *fd, cram_container *c,
					   cram_block_compression_hdr *h);

/*! INTERNAL:
 * Encodes a slice compression header. 
 *
 * @return
 * Returns cram_block on success;
 *         NULL on failure
 */
cram_block *cram_encode_slice_header(cram_fd *fd, cram_slice *s);

/*! INTERNAL:
 * Encodes all slices in a container into blocks.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 *
 * FIXME: separate into encode_container and write_container. Ideally
 * we should be able to do read_container / write_container or
 * decode_container / encode_container.
 */
int cram_encode_container(cram_fd *fd, cram_container *c);

#ifdef __cplusplus
}
#endif

#endif
