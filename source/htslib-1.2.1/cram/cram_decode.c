/*
Copyright (c) 2012-2014 Genome Research Ltd.
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

/*
 * - In-memory decoding of CRAM data structures.
 * - Iterator for reading CRAM record by record.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <ctype.h>

#include "cram/cram.h"
#include "cram/os.h"
#include "cram/md5.h"

//Whether CIGAR has just M or uses = and X to indicate match and mismatch
//#define USE_X

/* ----------------------------------------------------------------------
 * CRAM compression headers
 */

/*
 * Decodes the Tag Dictionary record in the preservation map
 * Updates the cram compression header.
 * 
 * Returns number of bytes decoded on success
 *        -1 on failure
 */
int cram_decode_TD(char *cp, cram_block_compression_hdr *h) {
    char *op = cp;
    unsigned char *dat;
    cram_block *b;
    int32_t blk_size;
    int nTL, i, sz;

	if (!(b = cram_new_block((cram_content_type)0, 0)))
	return -1;
    h->TD_blk = b;

    /* Decode */
    cp += itf8_get(cp, &blk_size);
    if (!blk_size) {
	h->nTL = 0;
	h->TL = NULL;
	cram_free_block(b);
        return cp - op;
    }

    BLOCK_APPEND(b, cp, blk_size);
    cp += blk_size;
    sz = cp - op;

    // Force nul termination if missing
    if (BLOCK_DATA(b)[BLOCK_SIZE(b)-1])
	BLOCK_APPEND_CHAR(b, '\0');

    /* Set up TL lookup table */
    dat = BLOCK_DATA(b);

    // Count
    for (nTL = i = 0; i < BLOCK_SIZE(b); i++) {
	nTL++;
	while (dat[i])
	    i++;
    }

    // Copy
    h->nTL = nTL;
	if (!(h->TL = (unsigned char**)calloc(h->nTL, sizeof(unsigned char *))))
	return -1;
    for (nTL = i = 0; i < BLOCK_SIZE(b); i++) {
	h->TL[nTL++] = &dat[i];
	while (dat[i])
	    i++;
    }
    
    return sz;
}

/*
 * Decodes a CRAM block compression header.
 * Returns header ptr on success
 *         NULL on failure
 */
cram_block_compression_hdr *cram_decode_compression_header(cram_fd *fd,
							   cram_block *b) {
    char *cp, *cp_copy;
	cram_block_compression_hdr *hdr = (cram_block_compression_hdr*)calloc(1, sizeof(*hdr));
    int i;
    int32_t map_size, map_count;

    if (!hdr)
	return NULL;

    if (b->method != RAW) {
	if (cram_uncompress_block(b)) {
	    free(hdr);
	    return NULL;
	}
    }

    cp = (char *)b->data;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
	cp += itf8_get(cp, &hdr->ref_seq_id);
	cp += itf8_get(cp, &hdr->ref_seq_start);
	cp += itf8_get(cp, &hdr->ref_seq_span);
	cp += itf8_get(cp, &hdr->num_records);
	cp += itf8_get(cp, &hdr->num_landmarks);
	if (!(hdr->landmark = (int32_t*)malloc(hdr->num_landmarks * sizeof(int32_t)))) {
	    free(hdr);
	    return NULL;
	}
	for (i = 0; i < hdr->num_landmarks; i++) {
	    cp += itf8_get(cp, &hdr->landmark[i]);
	}
    }

    hdr->preservation_map = kh_init(map);

    memset(hdr->rec_encoding_map, 0,
	   CRAM_MAP_HASH * sizeof(hdr->rec_encoding_map[0]));
    memset(hdr->tag_encoding_map, 0,
	   CRAM_MAP_HASH * sizeof(hdr->tag_encoding_map[0]));

    if (!hdr->preservation_map) {
	cram_free_compression_header(hdr);
	return NULL;
    }

    /* Initialise defaults for preservation map */
    hdr->mapped_qs_included = 0;
    hdr->unmapped_qs_included = 0;
    hdr->unmapped_placed = 0;
    hdr->qs_included = 0;
    hdr->read_names_included = 0;
    hdr->AP_delta = 1;
    memcpy(hdr->substitution_matrix, "CGTNAGTNACTNACGNACGT", 20);

    /* Preservation map */
    cp += itf8_get(cp, &map_size); cp_copy = cp;
    cp += itf8_get(cp, &map_count);
    for (i = 0; i < map_count; i++) {
	pmap_t hd;
	khint_t k;
	int r;

	cp += 2;
	switch(CRAM_KEY(cp[-2],cp[-1])) {
	case CRAM_KEY('M','I'):
	    hd.i = *cp++;
	    k = kh_put(map, hdr->preservation_map, "MI", &r);
	    if (-1 == r) {
		cram_free_compression_header(hdr);
                return NULL;
            }

	    kh_val(hdr->preservation_map, k) = hd;
	    hdr->mapped_qs_included = hd.i;
	    break;

	case CRAM_KEY('U','I'):
	    hd.i = *cp++;
	    k = kh_put(map, hdr->preservation_map, "UI", &r);
	    if (-1 == r) {
		cram_free_compression_header(hdr);
                return NULL;
            }

	    kh_val(hdr->preservation_map, k) = hd;
	    hdr->unmapped_qs_included = hd.i;
	    break;

	case CRAM_KEY('P','I'):
	    hd.i = *cp++;
	    k = kh_put(map, hdr->preservation_map, "PI", &r);
	    if (-1 == r) {
		cram_free_compression_header(hdr);
                return NULL;
            }

	    kh_val(hdr->preservation_map, k) = hd;
	    hdr->unmapped_placed = hd.i;
	    break;

	case CRAM_KEY('R','N'):
	    hd.i = *cp++;
	    k = kh_put(map, hdr->preservation_map, "RN", &r);
	    if (-1 == r) {
		cram_free_compression_header(hdr);
                return NULL;
            }

	    kh_val(hdr->preservation_map, k) = hd;
	    hdr->read_names_included = hd.i;
	    break;

	case CRAM_KEY('A','P'):
	    hd.i = *cp++;
	    k = kh_put(map, hdr->preservation_map, "AP", &r);
	    if (-1 == r) {
		cram_free_compression_header(hdr);
                return NULL;
            }

	    kh_val(hdr->preservation_map, k) = hd;
	    hdr->AP_delta = hd.i;
	    break;

	case CRAM_KEY('R','R'):
	    hd.i = *cp++;
	    k = kh_put(map, hdr->preservation_map, "RR", &r);
	    if (-1 == r) {
		cram_free_compression_header(hdr);
                return NULL;
            }

	    kh_val(hdr->preservation_map, k) = hd;
	    fd->no_ref = !hd.i;
	    break;

	case CRAM_KEY('S','M'):
	    hdr->substitution_matrix[0][(cp[0]>>6)&3] = 'C';
	    hdr->substitution_matrix[0][(cp[0]>>4)&3] = 'G';
	    hdr->substitution_matrix[0][(cp[0]>>2)&3] = 'T';
	    hdr->substitution_matrix[0][(cp[0]>>0)&3] = 'N';

	    hdr->substitution_matrix[1][(cp[1]>>6)&3] = 'A';
	    hdr->substitution_matrix[1][(cp[1]>>4)&3] = 'G';
	    hdr->substitution_matrix[1][(cp[1]>>2)&3] = 'T';
	    hdr->substitution_matrix[1][(cp[1]>>0)&3] = 'N';

	    hdr->substitution_matrix[2][(cp[2]>>6)&3] = 'A';
	    hdr->substitution_matrix[2][(cp[2]>>4)&3] = 'C';
	    hdr->substitution_matrix[2][(cp[2]>>2)&3] = 'T';
	    hdr->substitution_matrix[2][(cp[2]>>0)&3] = 'N';

	    hdr->substitution_matrix[3][(cp[3]>>6)&3] = 'A';
	    hdr->substitution_matrix[3][(cp[3]>>4)&3] = 'C';
	    hdr->substitution_matrix[3][(cp[3]>>2)&3] = 'G';
	    hdr->substitution_matrix[3][(cp[3]>>0)&3] = 'N';

	    hdr->substitution_matrix[4][(cp[4]>>6)&3] = 'A';
	    hdr->substitution_matrix[4][(cp[4]>>4)&3] = 'C';
	    hdr->substitution_matrix[4][(cp[4]>>2)&3] = 'G';
	    hdr->substitution_matrix[4][(cp[4]>>0)&3] = 'T';

	    hd.p = cp;
	    cp += 5;

	    k = kh_put(map, hdr->preservation_map, "SM", &r);
	    if (-1 == r) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	    kh_val(hdr->preservation_map, k) = hd;
	    break;

	case CRAM_KEY('T','D'): {
	    int sz = cram_decode_TD(cp, hdr); // tag dictionary
	    if (sz < 0) {
		cram_free_compression_header(hdr);
		return NULL;
	    }

	    hd.p = cp;
	    cp += sz;

	    k = kh_put(map, hdr->preservation_map, "TD", &r);
	    if (-1 == r) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	    kh_val(hdr->preservation_map, k) = hd;
	    break;
	}

	default:
	    fprintf(stderr, "Unrecognised preservation map key %c%c\n",
		    cp[-2], cp[-1]);
	    // guess byte;
	    cp++;
	    break;
	}
    }
    if (cp - cp_copy != map_size) {
	cram_free_compression_header(hdr);
	return NULL;
    }

    /* Record encoding map */
    cp += itf8_get(cp, &map_size); cp_copy = cp;
    cp += itf8_get(cp, &map_count);
    for (i = 0; i < map_count; i++) {
	char *key = cp;
	int32_t encoding;
	int32_t size;
	cram_map *m = (cram_map*) malloc(sizeof(*m)); // FIXME: use pooled_alloc

	if (!m) {
	    cram_free_compression_header(hdr);
	    return NULL;
	}

	cp += 2;
	cp += itf8_get(cp, &encoding);
	cp += itf8_get(cp, &size);

	// Fill out cram_map purely for cram_dump to dump out.
	m->key = (key[0]<<8)|key[1];
	m->encoding = (cram_encoding)encoding;
	m->size     = size;
	m->offset   = cp - (char *)b->data;
	m->codec = NULL;

	if (m->encoding == E_NULL)
	    continue;

	//printf("%s codes for %.2s\n", cram_encoding2str(encoding), key);

	/*
	 * For CRAM1.0 CF and BF are Byte and not Int.
	 * Practically speaking it makes no difference unless we have a
	 * 1.0 format file that stores these in EXTERNAL as only then
	 * does Byte vs Int matter.
	 *
	 * Neither this C code nor Java reference implementations did this,
	 * so we gloss over it and treat them as int.
	 */

	if (key[0] == 'B' && key[1] == 'F') {
	    if (!(hdr->codecs[DS_BF] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'C' && key[1] == 'F') {
		if (!(hdr->codecs[DS_CF] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'I') {
		if (!(hdr->codecs[DS_RI] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'L') {
		if (!(hdr->codecs[DS_RL] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'A' && key[1] == 'P') {
		if (!(hdr->codecs[DS_AP] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'G') {
		if (!(hdr->codecs[DS_RG] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'M' && key[1] == 'F') {
		if (!(hdr->codecs[DS_MF] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'N' && key[1] == 'S') {
		if (!(hdr->codecs[DS_NS] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'N' && key[1] == 'P') {
		if (!(hdr->codecs[DS_NP] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'S') {
		if (!(hdr->codecs[DS_TS] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'N' && key[1] == 'F') {
		if (!(hdr->codecs[DS_NF] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'C') {
		if (!(hdr->codecs[DS_TC] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_BYTE,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'N') {
		if (!(hdr->codecs[DS_TN] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'F' && key[1] == 'N') {
		if (!(hdr->codecs[DS_FN] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'F' && key[1] == 'C') {
		if (!(hdr->codecs[DS_FC] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_BYTE,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'F' && key[1] == 'P') {
		if (!(hdr->codecs[DS_FP] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'B' && key[1] == 'S') {
		if (!(hdr->codecs[DS_BS] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_BYTE,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'I' && key[1] == 'N') {
		if (!(hdr->codecs[DS_IN] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_BYTE_ARRAY,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'S' && key[1] == 'C') {
		if (!(hdr->codecs[DS_SC] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_BYTE_ARRAY,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'D' && key[1] == 'L') {
		if (!(hdr->codecs[DS_DL] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'B' && key[1] == 'A') {
		if (!(hdr->codecs[DS_BA] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_BYTE,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'B' && key[1] == 'B') {
		if (!(hdr->codecs[DS_BB] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_BYTE_ARRAY,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'S') {
		if (!(hdr->codecs[DS_RS] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'P' && key[1] == 'D') {
		if (!(hdr->codecs[DS_PD] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'H' && key[1] == 'C') {
		if (!(hdr->codecs[DS_HC] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'M' && key[1] == 'Q') {
		if (!(hdr->codecs[DS_MQ] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'N') {
		if (!(hdr->codecs[DS_RN] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_BYTE_ARRAY_BLOCK,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'Q' && key[1] == 'S') {
		if (!(hdr->codecs[DS_QS] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_BYTE,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'Q' && key[1] == 'Q') {
		if (!(hdr->codecs[DS_QQ] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_BYTE_ARRAY,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'L') {
		if (!(hdr->codecs[DS_TL] = cram_decoder_init((cram_encoding)encoding, cp, size,
							 E_INT,
							 fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'M') {
	} else if (key[0] == 'T' && key[1] == 'V') {
	} else
	    fprintf(stderr, "Unrecognised key: %.2s\n", key);

	cp += size;

	m->next = hdr->rec_encoding_map[CRAM_MAP(key[0], key[1])];
	hdr->rec_encoding_map[CRAM_MAP(key[0], key[1])] = m;
    }
    if (cp - cp_copy != map_size) {
	cram_free_compression_header(hdr);
	return NULL;
    }

    /* Tag encoding map */
    cp += itf8_get(cp, &map_size); cp_copy = cp;
    cp += itf8_get(cp, &map_count);
    for (i = 0; i < map_count; i++) {
	int32_t encoding;
	int32_t size;
	cram_map *m = (cram_map*)malloc(sizeof(*m)); // FIXME: use pooled_alloc
	char *key = cp+1;

	if (!m) {
	    cram_free_compression_header(hdr);
	    return NULL;
	}

	m->key = (key[0]<<16)|(key[1]<<8)|key[2];

	cp += 4; // Strictly ITF8, but this suffices
	cp += itf8_get(cp, &encoding);
	cp += itf8_get(cp, &size);

	m->encoding = (cram_encoding)encoding;
	m->size     = size;
	m->offset   = cp - (char *)b->data;
	if (!(m->codec = cram_decoder_init((cram_encoding)encoding, cp, size,
					   E_BYTE_ARRAY_BLOCK, fd->version))) {
	    cram_free_compression_header(hdr);
	    free(m);
	    return NULL;
	}
	
	cp += size;

	m->next = hdr->tag_encoding_map[CRAM_MAP(key[0],key[1])];
	hdr->tag_encoding_map[CRAM_MAP(key[0],key[1])] = m;
    }
    if (cp - cp_copy != map_size) {
	cram_free_compression_header(hdr);
	return NULL;
    }

    return hdr;
}

/*
 * Note we also need to scan through the record encoding map to
 * see which data series share the same block, either external or
 * CORE. For example if we need the BF data series but MQ and CF
 * are also encoded in the same block then we need to add those in
 * as a dependency in order to correctly decode BF.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_dependent_data_series(cram_fd *fd,
			       cram_block_compression_hdr *hdr,
			       cram_slice *s) {
    int *block_used;
    int core_used = 0;
    int i;
    static int i_to_id[] = {
	DS_BF, DS_AP, DS_FP, DS_RL, DS_DL, DS_NF, DS_BA, DS_QS,
	DS_FC, DS_FN, DS_BS, DS_IN, DS_RG, DS_MQ, DS_TL, DS_RN,
	DS_NS, DS_NP, DS_TS, DS_MF, DS_CF, DS_RI, DS_RS, DS_PD,
	DS_HC, DS_SC, DS_BB, DS_QQ,
    };
    uint32_t orig_ds;

    /*
     * Set the data_series bit field based on fd->required_fields
     * contents.
     */
    if (fd->required_fields && fd->required_fields != INT_MAX) {
	hdr->data_series = 0;

	if (fd->required_fields & SAM_QNAME)
	    hdr->data_series |= CRAM_RN;

	if (fd->required_fields & SAM_FLAG)
	    hdr->data_series |= CRAM_BF;

	if (fd->required_fields & SAM_RNAME)
	    hdr->data_series |= CRAM_RI | CRAM_BF;

	if (fd->required_fields & SAM_POS)
	    hdr->data_series |= CRAM_AP | CRAM_BF;

	if (fd->required_fields & SAM_MAPQ)
	    hdr->data_series |= CRAM_MQ;

	if (fd->required_fields & SAM_CIGAR)
	    hdr->data_series |= CRAM_CIGAR;

	if (fd->required_fields & SAM_RNEXT)
	    hdr->data_series |= CRAM_CF | CRAM_NF | CRAM_RI | CRAM_NS |CRAM_BF;

	if (fd->required_fields & SAM_PNEXT)
	    hdr->data_series |= CRAM_CF | CRAM_NF | CRAM_AP | CRAM_NP | CRAM_BF;

	if (fd->required_fields & SAM_TLEN)
	    hdr->data_series |= CRAM_CF | CRAM_NF | CRAM_AP | CRAM_TS |
		CRAM_BF | CRAM_MF | CRAM_RI | CRAM_CIGAR;

	if (fd->required_fields & SAM_SEQ)
	    hdr->data_series |= CRAM_SEQ;

	if (!(fd->required_fields & SAM_AUX))
	    // No easy way to get MD/NM without other tags at present
	    fd->decode_md = 0;

	if (fd->required_fields & SAM_QUAL)
	    hdr->data_series |= CRAM_SEQ;

	if (fd->required_fields & SAM_AUX)
	    hdr->data_series |= CRAM_RG | CRAM_TL | CRAM_aux;

	if (fd->required_fields & SAM_RGAUX)
	    hdr->data_series |= CRAM_RG | CRAM_BF;

	// Always uncompress CORE block
	if (cram_uncompress_block(s->block[0]))
	    return -1;
    } else {
	hdr->data_series = CRAM_ALL;

	for (i = 0; i < s->hdr->num_blocks; i++) {
	    if (cram_uncompress_block(s->block[i]))
		return -1;
	}

	return 0;
    }

    block_used = (int*)calloc(s->hdr->num_blocks+1, sizeof(int));
    if (!block_used)
	return -1;

    do {
	/*
	 * Also set data_series based on code prerequisites. Eg if we need
	 * CRAM_QS then we also need to know CRAM_RL so we know how long it
	 * is, or if we need FC/FP then we also need FN (number of features).
	 *
	 * It's not reciprocal though. We may be needing to decode FN
	 * but have no need to decode FC, FP and cigar ops.
	 */
	if (hdr->data_series & CRAM_RS)    hdr->data_series |= CRAM_FC|CRAM_FP;
	if (hdr->data_series & CRAM_PD)    hdr->data_series |= CRAM_FC|CRAM_FP;
	if (hdr->data_series & CRAM_HC)    hdr->data_series |= CRAM_FC|CRAM_FP;
	if (hdr->data_series & CRAM_QS)    hdr->data_series |= CRAM_FC|CRAM_FP;
	if (hdr->data_series & CRAM_IN)    hdr->data_series |= CRAM_FC|CRAM_FP;
	if (hdr->data_series & CRAM_SC)    hdr->data_series |= CRAM_FC|CRAM_FP;
	if (hdr->data_series & CRAM_BS)    hdr->data_series |= CRAM_FC|CRAM_FP;
	if (hdr->data_series & CRAM_DL)    hdr->data_series |= CRAM_FC|CRAM_FP;
	if (hdr->data_series & CRAM_BA)    hdr->data_series |= CRAM_FC|CRAM_FP;
	if (hdr->data_series & CRAM_BB)    hdr->data_series |= CRAM_FC|CRAM_FP;
	if (hdr->data_series & CRAM_QQ)    hdr->data_series |= CRAM_FC|CRAM_FP;

	// cram_decode_seq() needs seq[] array
	if (hdr->data_series & (CRAM_SEQ|CRAM_CIGAR)) hdr->data_series |= CRAM_RL;

	if (hdr->data_series & CRAM_FP)    hdr->data_series |= CRAM_FC;
	if (hdr->data_series & CRAM_FC)    hdr->data_series |= CRAM_FN;
	if (hdr->data_series & CRAM_aux)   hdr->data_series |= CRAM_TL;
	if (hdr->data_series & CRAM_MF)    hdr->data_series |= CRAM_CF;
	if (hdr->data_series & CRAM_MQ)    hdr->data_series |= CRAM_BF;
	if (hdr->data_series & CRAM_BS)    hdr->data_series |= CRAM_RI;
	if (hdr->data_series & (CRAM_MF |CRAM_NS |CRAM_NP |CRAM_TS |CRAM_NF))
	    hdr->data_series |= CRAM_CF;
	if (!hdr->read_names_included && hdr->data_series & CRAM_RN)
	    hdr->data_series |= CRAM_CF | CRAM_NF;
	if (hdr->data_series & (CRAM_BA | CRAM_QS | CRAM_BB | CRAM_QQ))
	    hdr->data_series |= CRAM_BF | CRAM_CF | CRAM_RL;

	orig_ds = hdr->data_series;

	// Find which blocks are in use.
	for (i = 0; i < sizeof(i_to_id)/sizeof(*i_to_id); i++) {
	    int bnum1, bnum2, j;
	    cram_codec *c = hdr->codecs[i_to_id[i]];

	    if (!(hdr->data_series & (1<<i)))
		continue;

	    if (!c)
		continue;

	    bnum1 = cram_codec_to_id(c, &bnum2);

	    for (;;) {
		switch (bnum1) {
		case -2:
		    break;

		case -1:
		    core_used = 1;
		    break;

		default:
		    for (j = 0; j < s->hdr->num_blocks; j++) {
			if (s->block[j]->content_type == EXTERNAL &&
			    s->block[j]->content_id == bnum1) {
			    block_used[j] = 1;
			    if (cram_uncompress_block(s->block[j])) {
				free(block_used);
				return -1;
			    }
			}
		    }
		    break;
		}

		if (bnum2 == -2 || bnum1 == bnum2)
		    break;

		bnum1 = bnum2; // 2nd pass
	    }
	}

	// Tags too
	if ((fd->required_fields & SAM_AUX) ||
	    (hdr->data_series & CRAM_aux)) {
	    for (i = 0; i < CRAM_MAP_HASH; i++) {
		int bnum1, bnum2, j;
		cram_map *m = hdr->tag_encoding_map[i];

		while (m) {
		    cram_codec *c = m->codec;
		    if (!c)
			continue;

		    bnum1 = cram_codec_to_id(c, &bnum2);

		    for (;;) {
			switch (bnum1) {
			case -2:
			    break;

			case -1:
			    core_used = 1;
			    break;

			default:
			    for (j = 0; j < s->hdr->num_blocks; j++) {
				if (s->block[j]->content_type == EXTERNAL &&
				    s->block[j]->content_id == bnum1) {
				    block_used[j] = 1;
				    if (cram_uncompress_block(s->block[j])) {
					free(block_used);
					return -1;
				    }
				}
			    }
			    break;
			}

			if (bnum2 == -2 || bnum1 == bnum2)
			    break;

			bnum1 = bnum2; // 2nd pass
		    }

		    m = m->next;
		}
	    }
	}

	// We now know which blocks are in used, so repeat and find
	// which other data series need to be added.
	for (i = 0; i < sizeof(i_to_id)/sizeof(*i_to_id); i++) {
	    int bnum1, bnum2, j;
	    cram_codec *c = hdr->codecs[i_to_id[i]];

	    if (!c)
		continue;

	    bnum1 = cram_codec_to_id(c, &bnum2);

	    for (;;) {
		switch (bnum1) {
		case -2:
		    break;

		case -1:
		    if (core_used) {
			//printf(" + data series %08x:\n", 1<<i);
			hdr->data_series |= 1<<i;
		    }
		    break;

		default:
		    for (j = 0; j < s->hdr->num_blocks; j++) {
			if (s->block[j]->content_type == EXTERNAL &&
			    s->block[j]->content_id == bnum1) {
			    if (block_used[j]) {
				//printf(" + data series %08x:\n", 1<<i);
				hdr->data_series |= 1<<i;
			    }
			}
		    }
		    break;
		}

		if (bnum2 == -2 || bnum1 == bnum2)
		    break;

		bnum1 = bnum2; // 2nd pass
	    }
	}

	// Tags too
	for (i = 0; i < CRAM_MAP_HASH; i++) {
	    int bnum1, bnum2, j;
	    cram_map *m = hdr->tag_encoding_map[i];

	    while (m) {
		cram_codec *c = m->codec;
		if (!c)
		    continue;

		bnum1 = cram_codec_to_id(c, &bnum2);
		
		for (;;) {
		    switch (bnum1) {
		    case -2:
			break;

		    case -1:
			//printf(" + data series %08x:\n", CRAM_aux);
			hdr->data_series |= CRAM_aux;
			break;

		    default:
			for (j = 0; j < s->hdr->num_blocks; j++) {
			    if (s->block[j]->content_type &&
				s->block[j]->content_id == bnum1) {
				if (block_used[j]) {
				    //printf(" + data series %08x:\n",
				    //       CRAM_aux);
				    hdr->data_series |= CRAM_aux;
				}
			    }
			}
			break;
		    }

		    if (bnum2 == -2 || bnum1 == bnum2)
			break;

		    bnum1 = bnum2; // 2nd pass
		}

		m = m->next;
	    }
	}
    } while (orig_ds != hdr->data_series);

    free(block_used);
    return 0;
}

/* ----------------------------------------------------------------------
 * CRAM slices
 */

/*
 * Decodes a CRAM (un)mapped slice header block.
 * Returns slice header ptr on success
 *         NULL on failure
 */
cram_block_slice_hdr *cram_decode_slice_header(cram_fd *fd, cram_block *b) {
    cram_block_slice_hdr *hdr;
    char *cp = (char *)b->data;
    int i;

    if (b->content_type != MAPPED_SLICE &&
	b->content_type != UNMAPPED_SLICE)
	return NULL;

	if (!(hdr = (cram_block_slice_hdr*)calloc(1, sizeof(*hdr))))
	return NULL;

    hdr->content_type = b->content_type;

    if (b->content_type == MAPPED_SLICE) {
	cp += itf8_get(cp, &hdr->ref_seq_id);
	cp += itf8_get(cp, &hdr->ref_seq_start);
	cp += itf8_get(cp, &hdr->ref_seq_span);
    }
    cp += itf8_get(cp, &hdr->num_records);
    hdr->record_counter = 0;
    if (CRAM_MAJOR_VERS(fd->version) == 2) {
	int32_t i32;
	cp += itf8_get(cp, &i32);
	hdr->record_counter = i32;
    } else if (CRAM_MAJOR_VERS(fd->version) >= 3) {
	cp += ltf8_get(cp, &hdr->record_counter);
    }

    cp += itf8_get(cp, &hdr->num_blocks);

    cp += itf8_get(cp, &hdr->num_content_ids);
    hdr->block_content_ids = (int32_t*)malloc(hdr->num_content_ids * sizeof(int32_t));
    if (!hdr->block_content_ids) {
	free(hdr);
	return NULL;
    }

    for (i = 0; i < hdr->num_content_ids; i++) {
	cp += itf8_get(cp, &hdr->block_content_ids[i]);
    }

    if (b->content_type == MAPPED_SLICE) {
	cp += itf8_get(cp, &hdr->ref_base_id);
    }

    if (CRAM_MAJOR_VERS(fd->version) != 1) {
	memcpy(hdr->md5, cp, 16);
    } else {
	memset(hdr->md5, 0, 16);
    }

    return hdr;
}


#if 0
/* Returns the number of bits set in val; it the highest bit used */
static int nbits(int v) {
    static const int MultiplyDeBruijnBitPosition[32] = {
	1, 10, 2, 11, 14, 22, 3, 30, 12, 15, 17, 19, 23, 26, 4, 31,
	9, 13, 21, 29, 16, 18, 25, 8, 20, 28, 24, 7, 27, 6, 5, 32
    };

    v |= v >> 1; // first up to set all bits 1 after the first 1 */
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;

    // DeBruijn magic to find top bit
    return MultiplyDeBruijnBitPosition[(uint32_t)(v * 0x07C4ACDDU) >> 27];
}
#endif

#if 0
static int sort_freqs(const void *vp1, const void *vp2) {
    const int i1 = *(const int *)vp1;
    const int i2 = *(const int *)vp2;
    return i1-i2;
}
#endif

/* ----------------------------------------------------------------------
 * Primary CRAM sequence decoder
 */

/*
 * Internal part of cram_decode_slice().
 * Generates the sequence, quality and cigar components.
 */
static int cram_decode_seq(cram_fd *fd, cram_container *c, cram_slice *s,
			   cram_block *blk, cram_record *cr, SAM_hdr *bfd,
			   int cf, char *seq, char *qual) {
    int prev_pos = 0, f, r = 0, out_sz = 1;
    int seq_pos = 1;
    int cig_len = 0, ref_pos = cr->apos;
    int32_t fn, i32;
	enum cigar_op cig_op = (cigar_op)BAM_CMATCH;
    uint32_t *cigar = s->cigar;
    uint32_t ncigar = s->ncigar;
    uint32_t cigar_alloc = s->cigar_alloc;
    uint32_t nm = 0, md_dist = 0;
    int orig_aux = 0;
    int decode_md = fd->decode_md && s->ref;
    uint32_t ds = c->comp_hdr->data_series;

    if ((ds & CRAM_QS) && !(cf & CRAM_FLAG_PRESERVE_QUAL_SCORES)) {
	memset(qual, 30, cr->len);
    }

    if (decode_md) {
	orig_aux = BLOCK_SIZE(s->aux_blk);
	BLOCK_APPEND(s->aux_blk, "MDZ", 3);
    }
    
    if (ds & CRAM_FN) {
	if (!c->comp_hdr->codecs[DS_FN]) return -1;
	r |= c->comp_hdr->codecs[DS_FN]->decode(s,c->comp_hdr->codecs[DS_FN],
						blk, (char *)&fn, &out_sz);
    } else {
	fn = 0;
    }

    ref_pos--; // count from 0
    cr->cigar = ncigar;

    if (!(ds & (CRAM_FC | CRAM_FP)))
	goto skip_cigar;

    for (f = 0; f < fn; f++) {
	int32_t pos = 0;
	char op;

	if (ncigar+2 >= cigar_alloc) {
	    cigar_alloc = cigar_alloc ? cigar_alloc*2 : 1024;
	    s->cigar = cigar;
	    if (!(cigar = (uint32_t*)realloc(cigar, cigar_alloc * sizeof(*cigar))))
		return -1;
	}

	if (ds & CRAM_FC) {
	    if (!c->comp_hdr->codecs[DS_FC]) return -1;
	    r |= c->comp_hdr->codecs[DS_FC]->decode(s,
						    c->comp_hdr->codecs[DS_FC],
						    blk,
						    &op,  &out_sz);
	}

	if (!(ds & CRAM_FP))
	    continue;

	if (!c->comp_hdr->codecs[DS_FP]) return -1;
	r |= c->comp_hdr->codecs[DS_FP]->decode(s,
						c->comp_hdr->codecs[DS_FP],
						blk,
						(char *)&pos, &out_sz);
	pos += prev_pos;

	if (pos > seq_pos) {
	    if (pos > cr->len+1)
		return -1;

	    if (s->ref && cr->ref_id >= 0) {
		if (ref_pos + pos - seq_pos > bfd->ref[cr->ref_id].len) {
		    static int whinged = 0;
		    if (!whinged)
			fprintf(stderr, "Ref pos outside of ref "
				"sequence boundary\n");
		    whinged = 1;
		} else {
		    memcpy(&seq[seq_pos-1], &s->ref[ref_pos - s->ref_start +1],
			   pos - seq_pos);
		}
	    }
#ifdef USE_X
	    if (cig_len && cig_op != BAM_CBASE_MATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    cig_op = BAM_CBASE_MATCH;
#else
	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
		cig_op = (cigar_op)BAM_CMATCH;
#endif
	    cig_len += pos - seq_pos;
	    ref_pos += pos - seq_pos;
	    md_dist += pos - seq_pos;
	    seq_pos = pos;
	}

	prev_pos = pos;

	if (!(ds & CRAM_FC))
	    goto skip_cigar;

	if (!(ds & CRAM_FC))
	    continue;

	switch(op) {
	case 'S': { // soft clip: IN
	    int32_t out_sz2 = 1;

	    if (cig_len) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_IN) {
		switch (CRAM_MAJOR_VERS(fd->version)) {
		case 1:
		    r |= c->comp_hdr->codecs[DS_IN]
			? c->comp_hdr->codecs[DS_IN]
			             ->decode(s, c->comp_hdr->codecs[DS_IN],
					      blk, &seq[pos-1], &out_sz2)
			: (seq[pos-1] = 'N', out_sz2 = 1, 0);
		    break;

		case 2:
		default:
		    r |= c->comp_hdr->codecs[DS_SC]
			? c->comp_hdr->codecs[DS_SC]
			             ->decode(s, c->comp_hdr->codecs[DS_SC],
					      blk, &seq[pos-1], &out_sz2)
			: (seq[pos-1] = 'N', out_sz2 = 1, 0);
		    break;

//		default:
//		    r |= c->comp_hdr->codecs[DS_BB]
//			? c->comp_hdr->codecs[DS_BB]
//			             ->decode(s, c->comp_hdr->codecs[DS_BB],
//					      blk, &seq[pos-1], &out_sz2)
//			: (seq[pos-1] = 'N', out_sz2 = 1, 0);
		}
		cigar[ncigar++] = (out_sz2<<4) + BAM_CSOFT_CLIP;
		cig_op = (cigar_op)BAM_CSOFT_CLIP;
		seq_pos += out_sz2;
	    }
	    break;
	}

	case 'X': { // Substitution; BS
	    unsigned char base;
#ifdef USE_X
	    if (cig_len && cig_op != BAM_CBASE_MISMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_BS) {
		if (!c->comp_hdr->codecs[DS_BS]) return -1;
		r |= c->comp_hdr->codecs[DS_BS]
		                ->decode(s, c->comp_hdr->codecs[DS_BS], blk,
					 (char *)&base, &out_sz);
		seq[pos-1] = 'N'; // FIXME look up BS=base value
	    }
	    cig_op = BAM_CBASE_MISMATCH;
#else
	    int ref_base;
	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_BS) {
		if (!c->comp_hdr->codecs[DS_BS]) return -1;
		r |= c->comp_hdr->codecs[DS_BS]
		                ->decode(s, c->comp_hdr->codecs[DS_BS], blk,
					 (char *)&base, &out_sz);
		if (ref_pos >= bfd->ref[cr->ref_id].len || !s->ref) {
		    seq[pos-1] = 'N';
		} else {
		    ref_base = fd->L1[(uc)s->ref[ref_pos - s->ref_start +1]];
		    seq[pos-1] = c->comp_hdr->
			substitution_matrix[ref_base][base];
		    if (decode_md) {
			BLOCK_APPEND_UINT(s->aux_blk, md_dist);
			BLOCK_APPEND_CHAR(s->aux_blk,
					  s->ref[ref_pos-s->ref_start +1]);
			md_dist = 0;
		    }
		}
	    }
		cig_op = (cigar_op)BAM_CMATCH;
#endif
	    nm++;
	    cig_len++;
	    seq_pos++;
	    ref_pos++;
	    break;
	}

	case 'D': { // Deletion; DL
	    if (cig_len && cig_op != BAM_CDEL) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_DL) {
		if (!c->comp_hdr->codecs[DS_DL]) return -1;
		r |= c->comp_hdr->codecs[DS_DL]
		                ->decode(s, c->comp_hdr->codecs[DS_DL], blk,
					 (char *)&i32, &out_sz);
		if (decode_md) {
		    BLOCK_APPEND_UINT(s->aux_blk, md_dist);
		    BLOCK_APPEND_CHAR(s->aux_blk, '^');
		    BLOCK_APPEND(s->aux_blk,
				 &s->ref[ref_pos - s->ref_start +1],
				 i32);
		    md_dist = 0;
		}
		cig_op = (cigar_op)BAM_CDEL;
		cig_len += i32;
		ref_pos += i32;
		nm      += i32;
		//printf("  %d: DL = %d (ret %d)\n", f, i32, r);
	    }
	    break;
	}

	case 'I': { // Insertion (several bases); IN
	    int32_t out_sz2 = 1;

	    if (cig_len && cig_op != BAM_CINS) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }

	    if (ds & CRAM_IN) {
		if (!c->comp_hdr->codecs[DS_IN]) return -1;
		r |= c->comp_hdr->codecs[DS_IN]
		                ->decode(s, c->comp_hdr->codecs[DS_IN], blk,
					 &seq[pos-1], &out_sz2);
		cig_op = (cigar_op)BAM_CINS;
		cig_len += out_sz2;
		seq_pos += out_sz2;
		nm      += out_sz2;
		//printf("  %d: IN(I) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, dat, r, out_sz2);
	    }
	    break;
	}

	case 'i': { // Insertion (single base); BA
	    if (cig_len && cig_op != BAM_CINS) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_BA) {
		if (!c->comp_hdr->codecs[DS_BA]) return -1;
		r |= c->comp_hdr->codecs[DS_BA]
		                ->decode(s, c->comp_hdr->codecs[DS_BA], blk,
					 (char *)&seq[pos-1], &out_sz);
	    }
		cig_op = (cigar_op)BAM_CINS;
	    cig_len++;
	    seq_pos++;
	    nm++;
	    break;
	}

	case 'b': { // Several bases
	    int32_t len = 1;

	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }

	    if (ds & CRAM_BB) {
		if (!c->comp_hdr->codecs[DS_BB]) return -1;
		r |= c->comp_hdr->codecs[DS_BB]
		    ->decode(s, c->comp_hdr->codecs[DS_BB], blk,
			     (char *)&seq[pos-1], &len);
	    }

		cig_op = (cigar_op)BAM_CMATCH;

	    cig_len+=len;
	    seq_pos+=len;
	    ref_pos+=len;
	    //prev_pos+=len;
	    break;
	}

	case 'q': { // Several quality values
	    int32_t len = 1;

	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }

	    if (ds & CRAM_QQ) {
		if (!c->comp_hdr->codecs[DS_QQ]) return -1;
		r |= c->comp_hdr->codecs[DS_QQ]
		    ->decode(s, c->comp_hdr->codecs[DS_QQ], blk,
			     (char *)&qual[pos-1], &len);
	    }

		cig_op = (cigar_op)BAM_CMATCH;

	    cig_len+=len;
	    seq_pos+=len;
	    ref_pos+=len;
	    //prev_pos+=len;
	    break;
	}

	case 'B': { // Read base; BA, QS
#ifdef USE_X
	    if (cig_len && cig_op != BAM_CBASE_MISMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
#else
	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
#endif
	    if (ds & CRAM_BA) {
		if (!c->comp_hdr->codecs[DS_BA]) return -1;
		r |= c->comp_hdr->codecs[DS_BA]
		                ->decode(s, c->comp_hdr->codecs[DS_BA], blk,
					 (char *)&seq[pos-1], &out_sz);
	    }
	    if (ds & CRAM_QS) {
		if (!c->comp_hdr->codecs[DS_QS]) return -1;
		r |= c->comp_hdr->codecs[DS_QS]
		                ->decode(s, c->comp_hdr->codecs[DS_QS], blk,
					 (char *)&qual[pos-1], &out_sz);
	    }
#ifdef USE_X
	    cig_op = BAM_CBASE_MISMATCH;
#else
		cig_op = (cigar_op)BAM_CMATCH;
#endif
	    cig_len++;
	    seq_pos++;
	    ref_pos++;
	    //printf("  %d: BA/QS(B) = %c/%d (ret %d)\n", f, i32, qc, r);
	    break;
	}

	case 'Q': { // Quality score; QS
	    if (ds & CRAM_QS) {
		if (!c->comp_hdr->codecs[DS_QS]) return -1;
		r |= c->comp_hdr->codecs[DS_QS]
		                ->decode(s, c->comp_hdr->codecs[DS_QS], blk,
					 (char *)&qual[pos-1], &out_sz);
		//printf("  %d: QS = %d (ret %d)\n", f, qc, r);
	    }
	    break;
	}

	case 'H': { // hard clip; HC
	    if (cig_len && cig_op != BAM_CHARD_CLIP) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_HC) {
		if (!c->comp_hdr->codecs[DS_HC]) return -1;
		r |= c->comp_hdr->codecs[DS_HC]
		                ->decode(s, c->comp_hdr->codecs[DS_HC], blk,
					 (char *)&i32, &out_sz);
		cig_op = (cigar_op)BAM_CHARD_CLIP;
		cig_len += i32;
		nm      += i32;
	    }
	    break;
	}

	case 'P': { // padding; PD
	    if (cig_len && cig_op != BAM_CPAD) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_PD) {
		if (!c->comp_hdr->codecs[DS_PD]) return -1;
		r |= c->comp_hdr->codecs[DS_PD]
		                ->decode(s, c->comp_hdr->codecs[DS_PD], blk,
					 (char *)&i32, &out_sz);
		cig_op = (cigar_op)BAM_CPAD;
		cig_len += i32;
		nm      += i32;
	    }
	    break;
	}

	case 'N': { // Ref skip; RS
	    if (cig_len && cig_op != BAM_CREF_SKIP) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (ds & CRAM_RS) {
		if (!c->comp_hdr->codecs[DS_RS]) return -1;
		r |= c->comp_hdr->codecs[DS_RS]
		                ->decode(s, c->comp_hdr->codecs[DS_RS], blk,
					 (char *)&i32, &out_sz);
		cig_op = (cigar_op)BAM_CREF_SKIP;
		cig_len += i32;
		ref_pos += i32;
		nm      += i32;
	    }
	    break;
	}

	default:
	    abort();
	}
    }

    if (!(ds & CRAM_FC))
	goto skip_cigar;

    /* An implement match op for any unaccounted for bases */
    if ((ds & CRAM_FN) && cr->len >= seq_pos) {
	if (s->ref) {
	    if (ref_pos + cr->len - seq_pos + 1 > bfd->ref[cr->ref_id].len) {
		static int whinged = 0;
		if (!whinged)
		    fprintf(stderr, "Ref pos outside of ref sequence boundary\n");
		whinged = 1;
	    } else {
		memcpy(&seq[seq_pos-1], &s->ref[ref_pos - s->ref_start +1],
		       cr->len - seq_pos + 1);
		ref_pos += cr->len - seq_pos + 1;
		md_dist += cr->len - seq_pos + 1;
	    }
	}

	if (ncigar+1 >= cigar_alloc) {
	    cigar_alloc = cigar_alloc ? cigar_alloc*2 : 1024;
	    s->cigar = cigar;
		if (!(cigar = (uint32_t*)realloc(cigar, cigar_alloc * sizeof(*cigar))))
		return -1;
	}
#ifdef USE_X
	if (cig_len && cig_op != BAM_CBASE_MATCH) {
	    cigar[ncigar++] = (cig_len<<4) + cig_op;
	    cig_len = 0;
	}
	cig_op = BAM_CBASE_MATCH;
#else
	if (cig_len && cig_op != BAM_CMATCH) {
	    cigar[ncigar++] = (cig_len<<4) + cig_op;
	    cig_len = 0;
	}
	cig_op = (cigar_op)BAM_CMATCH;
#endif
	cig_len += cr->len - seq_pos+1;
    }

 skip_cigar:

    if ((ds & CRAM_FN) && decode_md) {
	BLOCK_APPEND_UINT(s->aux_blk, md_dist);
    }

    if (cig_len) {
	if (ncigar >= cigar_alloc) {
	    cigar_alloc = cigar_alloc ? cigar_alloc*2 : 1024;
	    s->cigar = cigar;
		if (!(cigar = (uint32_t*)realloc(cigar, cigar_alloc * sizeof(*cigar))))
		return -1;
	}

	cigar[ncigar++] = (cig_len<<4) + cig_op;
    }

    cr->ncigar = ncigar - cr->cigar;
    cr->aend = ref_pos;

    //printf("2: %.*s %d .. %d\n", cr->name_len, DSTRING_STR(name_ds) + cr->name, cr->apos, ref_pos);

    if (ds & CRAM_MQ) {
	if (!c->comp_hdr->codecs[DS_MQ]) return -1;
	r |= c->comp_hdr->codecs[DS_MQ]
	                ->decode(s, c->comp_hdr->codecs[DS_MQ], blk,
				 (char *)&cr->mqual, &out_sz);
    } else {
	cr->mqual = 40;
    }

    if ((ds & CRAM_QS) && (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES)) {
	int32_t out_sz2 = cr->len;

	if (ds & CRAM_QS) {
	    if (!c->comp_hdr->codecs[DS_QS]) return -1;
	    r |= c->comp_hdr->codecs[DS_QS]
		            ->decode(s, c->comp_hdr->codecs[DS_QS], blk,
				     qual, &out_sz2);
	}
    }

    s->cigar = cigar;
    s->cigar_alloc = cigar_alloc;
    s->ncigar = ncigar;

    if (decode_md) {
	char buf[7];
	BLOCK_APPEND_CHAR(s->aux_blk, '\0'); // null terminate MD:Z:
	cr->aux_size += BLOCK_SIZE(s->aux_blk) - orig_aux;
	buf[0] = 'N'; buf[1] = 'M'; buf[2] = 'I';
	buf[3] = (nm>> 0) & 0xff;
	buf[4] = (nm>> 8) & 0xff;
	buf[5] = (nm>>16) & 0xff;
	buf[6] = (nm>>24) & 0xff;
	BLOCK_APPEND(s->aux_blk, buf, 7);
	cr->aux_size += 7;
    }

    return r;
}

/*
 * Quick and simple hash lookup for cram_map arrays
 */
static cram_map *map_find(cram_map **map, unsigned char *key, int id) {
    cram_map *m;

    m = map[CRAM_MAP(key[0],key[1])];
    while (m && m->key != id)
	m= m->next;

    return m;
}

//#define map_find(M,K,I) M[CRAM_MAP(K[0],K[1])];while (m && m->key != I);m= m->next


static int cram_decode_aux_1_0(cram_container *c, cram_slice *s,
			       cram_block *blk, cram_record *cr) {
    int i, r = 0, out_sz = 1;
    unsigned char ntags;
	    
    if (!c->comp_hdr->codecs[DS_TC]) return -1;
    r |= c->comp_hdr->codecs[DS_TC]->decode(s, c->comp_hdr->codecs[DS_TC], blk,
					    (char *)&ntags, &out_sz);
    cr->ntags = ntags;

    //printf("TC=%d\n", cr->ntags);
    cr->aux_size = 0;
    cr->aux = BLOCK_SIZE(s->aux_blk);

    for (i = 0; i < cr->ntags; i++) {
	int32_t id, out_sz = 1;
	unsigned char tag_data[3];
	cram_map *m;

	//printf("Tag %d/%d\n", i+1, cr->ntags);
	if (!c->comp_hdr->codecs[DS_TN]) return -1;
	r |= c->comp_hdr->codecs[DS_TN]->decode(s, c->comp_hdr->codecs[DS_TN],
						blk, (char *)&id, &out_sz);
	if (out_sz == 3) {
	    tag_data[0] = ((char *)&id)[0];
	    tag_data[1] = ((char *)&id)[1];
	    tag_data[2] = ((char *)&id)[2];
	} else {
	    tag_data[0] = (id>>16) & 0xff;
	    tag_data[1] = (id>>8)  & 0xff;
	    tag_data[2] = id       & 0xff;
	} 

	m = map_find(c->comp_hdr->tag_encoding_map, tag_data, id);
	if (!m)
	    return -1;
	BLOCK_APPEND(s->aux_blk, (char *)tag_data, 3);

	if (!m->codec) return -1;
	r |= m->codec->decode(s, m->codec, blk, (char *)s->aux_blk, &out_sz);

	cr->aux_size += out_sz + 3;
    }
    
    return r;
}

static int cram_decode_aux(cram_container *c, cram_slice *s,
			   cram_block *blk, cram_record *cr) {
    int i, r = 0, out_sz = 1;
    int32_t TL = 0;
    unsigned char *TN;
    uint32_t ds = c->comp_hdr->data_series;
	    
    if (!(ds & (CRAM_TL|CRAM_aux))) {
	cr->aux = 0;
	cr->aux_size = 0;
	return 0;
    }

    if (!c->comp_hdr->codecs[DS_TL]) return -1;
    r |= c->comp_hdr->codecs[DS_TL]->decode(s, c->comp_hdr->codecs[DS_TL], blk,
					    (char *)&TL, &out_sz);
    if (r || TL < 0 || TL >= c->comp_hdr->nTL)
	return -1;

    TN = c->comp_hdr->TL[TL];
    cr->ntags = strlen((char *)TN)/3; // optimise to remove strlen

    //printf("TC=%d\n", cr->ntags);
    cr->aux_size = 0;
    cr->aux = BLOCK_SIZE(s->aux_blk);

    if (!(ds & CRAM_aux))
	return 0;

    for (i = 0; i < cr->ntags; i++) {
	int32_t id, out_sz = 1;
	unsigned char tag_data[3];
	cram_map *m;

	//printf("Tag %d/%d\n", i+1, cr->ntags);
	tag_data[0] = *TN++;
	tag_data[1] = *TN++;
	tag_data[2] = *TN++;
	id = (tag_data[0]<<16) | (tag_data[1]<<8) | tag_data[2];

	m = map_find(c->comp_hdr->tag_encoding_map, tag_data, id);
	if (!m)
	    return -1;
	BLOCK_APPEND(s->aux_blk, (char *)tag_data, 3);

	if (!m->codec) return -1;
	r |= m->codec->decode(s, m->codec, blk, (char *)s->aux_blk, &out_sz);
	cr->aux_size += out_sz + 3;
    }
    
    return r;
}

/* Resolve mate pair cross-references between recs within this slice */
static void cram_decode_slice_xref(cram_slice *s, int required_fields) {
    int rec;

    if (!(required_fields & (SAM_RNEXT | SAM_PNEXT | SAM_TLEN))) {
	for (rec = 0; rec < s->hdr->num_records; rec++) {
	    cram_record *cr = &s->crecs[rec];

	    cr->tlen = 0;
	    cr->mate_pos = 0;
	    cr->mate_ref_id = -1;
	}

	return;
    }

    for (rec = 0; rec < s->hdr->num_records; rec++) {
	cram_record *cr = &s->crecs[rec];

	if (cr->mate_line >= 0) {
	    if (cr->mate_line < s->hdr->num_records) {
		/*
		 * On the first read, loop through computing lengths.
		 * It's not perfect as we have one slice per reference so we
		 * cannot detect when TLEN should be zero due to seqs that
		 * map to multiple references.
		 *
		 * We also cannot set tlen correct when it spans a slice for
		 * other reasons. This may make tlen too small. Should we
		 * fix this by forcing TLEN to be stored verbatim in such cases?
		 *
		 * Or do we just admit defeat and output 0 for tlen? It's the
		 * safe option...
		 */
		if (cr->tlen == INT_MIN) {
		    int id1 = rec, id2 = rec;
		    int aleft = cr->apos, aright = cr->aend;
		    int tlen;
		    int ref = cr->ref_id;

		    // number of segments starting at the same point.
		    int left_cnt = 0;

		    do {
			if (aleft > s->crecs[id2].apos)
			    aleft = s->crecs[id2].apos, left_cnt = 1;
			else if (aleft == s->crecs[id2].apos)
			    left_cnt++;
			if (aright < s->crecs[id2].aend)
			    aright = s->crecs[id2].aend;
			if (s->crecs[id2].mate_line == -1) {
			    s->crecs[id2].mate_line = rec;
			    break;
			}
			assert(s->crecs[id2].mate_line > id2);
			id2 = s->crecs[id2].mate_line;

			if (s->crecs[id2].ref_id != ref)
			    ref = -1;
		    } while (id2 != id1);

		    if (ref != -1) {
			tlen = aright - aleft + 1;
			id1 = id2 = rec;

			/*
			 * When we have two seqs with identical start and
			 * end coordinates, set +/- tlen based on 1st/last
			 * bit flags instead, as a tie breaker.
			 */
			if (s->crecs[id2].apos == aleft) {
			    if (left_cnt == 1 || 
				(s->crecs[id2].flags & BAM_FREAD1))
				s->crecs[id2].tlen = tlen;
			    else
				s->crecs[id2].tlen = -tlen;
			} else {
			    s->crecs[id2].tlen = -tlen;
			}

			id2 = s->crecs[id2].mate_line;
			while (id2 != id1) {
			    if (s->crecs[id2].apos == aleft) {
				if (left_cnt == 1 || 
				    (s->crecs[id2].flags & BAM_FREAD1))
				    s->crecs[id2].tlen = tlen;
				else
				    s->crecs[id2].tlen = -tlen;
			    } else {
				s->crecs[id2].tlen = -tlen;
			    }
			    id2 = s->crecs[id2].mate_line;
			}
		    } else {
			id1 = id2 = rec;

			s->crecs[id2].tlen = 0;
			id2 = s->crecs[id2].mate_line;
			while (id2 != id1) {
			    s->crecs[id2].tlen = 0;
			    id2 = s->crecs[id2].mate_line;
			}
		    }
		}

		cr->mate_pos = s->crecs[cr->mate_line].apos;
		cr->mate_ref_id = s->crecs[cr->mate_line].ref_id;

		// paired
		cr->flags |= BAM_FPAIRED;

		// set mate unmapped if needed
		if (s->crecs[cr->mate_line].flags & BAM_FUNMAP) {
		    cr->flags |= BAM_FMUNMAP;
		    cr->tlen = 0;
		}
		if (cr->flags & BAM_FUNMAP) {
		    cr->tlen = 0;
		}

		// set mate reversed if needed
		if (s->crecs[cr->mate_line].flags & BAM_FREVERSE)
		    cr->flags |= BAM_FMREVERSE;
	    } else {
		fprintf(stderr, "Mate line out of bounds: %d vs [0, %d]\n",
			cr->mate_line, s->hdr->num_records-1);
	    }

	    /* FIXME: construct read names here too if needed */
	} else {
	    if (cr->mate_flags & CRAM_M_REVERSE) {
		cr->flags |= BAM_FPAIRED | BAM_FMREVERSE;
	    }
	    if (cr->mate_flags & CRAM_M_UNMAP) {
		cr->flags |= BAM_FMUNMAP;
		//cr->mate_ref_id = -1;
	    }
	    if (!(cr->flags & BAM_FPAIRED))
		cr->mate_ref_id = -1;
	}

	if (cr->tlen == INT_MIN)
	    cr->tlen = 0; // Just incase
    }    
}

static char *md5_print(unsigned char *md5, char *out) {
    int i;
    for (i = 0; i < 16; i++) {
	out[i*2+0] = "0123456789abcdef"[md5[i]>>4];
	out[i*2+1] = "0123456789abcdef"[md5[i]&15];
    }
    out[32] = 0;

    return out;
}

typedef struct {
	cram_fd *fd;
	cram_container *c;
	cram_slice *s;
	SAM_hdr *h;
	int exit_code;
} cram_decode_job;

void *cram_decode_slice_thread(void *arg) {
	cram_decode_job *j = (cram_decode_job *)arg;

	j->exit_code = cram_decode_slice(j->fd, j->c, j->s, j->h);

	return j;
}

/*
* Spawn a multi-threaded version of cram_decode_slice().
*/
int cram_decode_slice_mt(cram_fd *fd, cram_container *c, cram_slice *s,
	SAM_hdr *bfd) {
	cram_decode_job *j;
	int nonblock;

	if (!fd->pool)
		return cram_decode_slice(fd, c, s, bfd);

	if (!(j = (cram_decode_job*)malloc(sizeof(*j))))
		return -1;

	j->fd = fd;
	j->c = c;
	j->s = s;
	j->h = bfd;

	nonblock = t_pool_results_queue_sz(fd->rqueue) ? 1 : 0;

	if (-1 == t_pool_dispatch2(fd->pool, fd->rqueue, cram_decode_slice_thread,
		j, nonblock)) {
		/* Would block */
		fd->job_pending = j;
	}
	else {
		fd->job_pending = NULL;
	}

	// flush too
	return 0;
}

/* ----------------------------------------------------------------------
* CRAM sequence iterators.
*/

/*
* Converts a cram in-memory record into a bam in-memory record. We
* pass a pointer to a bam_seq_t pointer along with the a pointer to
* the allocated size. These can initially be pointers to NULL and zero.
*
* This function will reallocate the bam buffer as required and update
* (*bam)->alloc accordingly, allowing it to be used within a loop
* efficiently without needing to allocate new bam objects over and
* over again.
*
* Returns the used size of the bam record on success
*         -1 on failure.
*/
static int cram_to_bam(SAM_hdr *bfd, cram_fd *fd, cram_slice *s,
	cram_record *cr, int rec, bam_seq_t **bam) {
	int bam_idx, rg_len;
	char name_a[1024], *name;
	int name_len;
	char *aux, *aux_orig;
	char *seq, *qual;

	/* Assign names if not explicitly set */
	if (fd->required_fields & SAM_QNAME) {
		if (cr->name_len) {
			name = (char *)BLOCK_DATA(s->name_blk) + cr->name;
			name_len = cr->name_len;
		}
		else {
			name = name_a;
			name_len = strlen(fd->prefix);
			memcpy(name, fd->prefix, name_len);
			name += name_len;
			*name++ = ':';
			if (cr->mate_line >= 0 && cr->mate_line < rec)
				name = (char *)append_uint64((unsigned char *)name,
				s->hdr->record_counter +
				cr->mate_line + 1);
			else
				name = (char *)append_uint64((unsigned char *)name,
				s->hdr->record_counter +
				rec + 1);
			name_len = name - name_a;
			name = name_a;
		}
	}
	else {
		name = "?";
		name_len = 1;
	}

	/* Generate BAM record */
	if (cr->rg < -1 || cr->rg >= bfd->nrg)
		return -1;
	rg_len = (cr->rg != -1) ? bfd->rg[cr->rg].name_len + 4 : 0;

	if (fd->required_fields & (SAM_SEQ | SAM_QUAL)) {
		if (!BLOCK_DATA(s->seqs_blk))
			return -1;
		seq = (char *)BLOCK_DATA(s->seqs_blk) + cr->seq;
	}
	else {
		seq = "*";
		cr->len = 1;
	}


	if (fd->required_fields & SAM_QUAL) {
		if (!BLOCK_DATA(s->qual_blk))
			return -1;
		qual = (char *)BLOCK_DATA(s->qual_blk) + cr->qual;
	}
	else {
		qual = NULL;
	}

	bam_idx = bam_construct_seq(bam, cr->aux_size + rg_len,
		name, name_len,
		cr->flags,
		cr->ref_id,
		cr->apos,
		cr->aend,
		cr->mqual,
		cr->ncigar, &s->cigar[cr->cigar],
		cr->mate_ref_id,
		cr->mate_pos,
		cr->tlen,
		cr->len,
		seq,
		qual);
	if (bam_idx == -1)
		return -1;

	aux = aux_orig = (char *)bam_aux(*bam);

	/* Auxiliary strings */
	if (cr->aux_size != 0) {
		memcpy(aux, BLOCK_DATA(s->aux_blk) + cr->aux, cr->aux_size);
		aux += cr->aux_size;
	}

	/* RG:Z: */
	if (cr->rg != -1) {
		int len = bfd->rg[cr->rg].name_len;
		*aux++ = 'R'; *aux++ = 'G'; *aux++ = 'Z';
		memcpy(aux, bfd->rg[cr->rg].name, len);
		aux += len;
		*aux++ = 0;
	}

	return bam_idx + (aux - aux_orig);
}

/*
* Read the next cram record and convert it to a bam_seq_t struct.
*
* Returns 0 on success
*        -1 on EOF or failure (check fd->err)
*/
int cram_get_bam_seq(cram_fd *fd, bam_seq_t **bam) {
	cram_record *cr;
	cram_container *c;
	cram_slice *s;

	if (!(cr = cram_get_seq(fd)))
		return -1;

	c = fd->ctr;
	s = c->slice;

	return cram_to_bam(fd->header, fd, s, cr, c->curr_rec - 1, bam);
}


#if !defined(_WIN32) && defined(USE_PTHREAD)

/*
* Decode an entire slice from container blocks. Fills out s->crecs[] array.
* Returns 0 on success
*        -1 on failure
*/
int cram_decode_slice(cram_fd *fd, cram_container *c, cram_slice *s,
	SAM_hdr *bfd) {
	cram_block *blk = s->block[0];
	int32_t bf, ref_id;
	unsigned char cf;
	int out_sz, r = 0;
	int rec;
	char *seq = NULL, *qual = NULL;
	int unknown_rg = -1;
	int embed_ref;
	char **refs = NULL;
	uint32_t ds;

	if (cram_dependent_data_series(fd, c->comp_hdr, s) != 0)
		return -1;

	ds = c->comp_hdr->data_series;

	blk->bit = 7; // MSB first

	/* Look for unknown RG, added as last by Java CRAM? */
	if (bfd->nrg > 0 &&
		!strcmp(bfd->rg[bfd->nrg - 1].name, "UNKNOWN"))
		unknown_rg = bfd->nrg - 1;

	if (blk->content_type != CORE)
		return -1;

	if (s->crecs)
		free(s->crecs);
	if (!(s->crecs = (cram_record*)malloc(s->hdr->num_records * sizeof(*s->crecs))))
		return -1;

	ref_id = s->hdr->ref_seq_id;
	embed_ref = s->hdr->ref_base_id >= 0 ? 1 : 0;

	if (ref_id >= 0) {
		if (embed_ref) {
			cram_block *b;
			if (s->hdr->ref_base_id < 0) {
				fprintf(stderr, "No reference specified and "
					"no embedded reference is available.\n");
				return -1;
			}
			if (!s->block_by_id ||
				!(b = s->block_by_id[s->hdr->ref_base_id]))
				return -1;
			cram_uncompress_block(b);
			s->ref = (char *)BLOCK_DATA(b);
			s->ref_start = s->hdr->ref_seq_start;
			s->ref_end = s->hdr->ref_seq_start + s->hdr->ref_seq_span - 1;
		}
		else if (!fd->no_ref) {
			//// Avoid Java cramtools bug by loading entire reference seq 
			//s->ref = cram_get_ref(fd, s->hdr->ref_seq_id, 1, 0);
			//s->ref_start = 1;

			if (fd->required_fields & SAM_SEQ)
				s->ref =
				cram_get_ref(fd, s->hdr->ref_seq_id,
				s->hdr->ref_seq_start,
				s->hdr->ref_seq_start + s->hdr->ref_seq_span - 1);
			s->ref_start = s->hdr->ref_seq_start;
			s->ref_end = s->hdr->ref_seq_start + s->hdr->ref_seq_span - 1;

			/* Sanity check */
			if (s->ref_start < 0) {
				fprintf(stderr, "Slice starts before base 1.\n");
				s->ref_start = 0;
			}

			pthread_mutex_lock(&fd->ref_lock);
			pthread_mutex_lock(&fd->refs->lock);
			if ((fd->required_fields & SAM_SEQ) &&
				s->ref_end > fd->refs->ref_id[ref_id]->length) {
				fprintf(stderr, "Slice ends beyond reference end.\n");
				s->ref_end = fd->refs->ref_id[ref_id]->length;
			}
			pthread_mutex_unlock(&fd->refs->lock);
			pthread_mutex_unlock(&fd->ref_lock);
		}
	}

	if ((fd->required_fields & SAM_SEQ) &&
		s->ref == NULL && s->hdr->ref_seq_id >= 0 && !fd->no_ref) {
		fprintf(stderr, "Unable to fetch reference #%d %d..%d\n",
			s->hdr->ref_seq_id, s->hdr->ref_seq_start,
			s->hdr->ref_seq_start + s->hdr->ref_seq_span - 1);
		return -1;
	}

	if (CRAM_MAJOR_VERS(fd->version) != 1
		&& (fd->required_fields & SAM_SEQ)
		&& s->hdr->ref_seq_id >= 0
		&& !fd->ignore_md5
		&& memcmp(s->hdr->md5, "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", 16)) {
		MD5_CTX md5;
		unsigned char digest[16];

		if (s->ref && s->hdr->ref_seq_id >= 0) {
			int start, len;

			if (s->hdr->ref_seq_start >= s->ref_start) {
				start = s->hdr->ref_seq_start - s->ref_start;
			}
			else {
				fprintf(stderr, "Slice starts before base 1.\n");
				start = 0;
			}

			if (s->hdr->ref_seq_span <= s->ref_end - s->ref_start + 1) {
				len = s->hdr->ref_seq_span;
			}
			else {
				fprintf(stderr, "Slice ends beyond reference end.\n");
				len = s->ref_end - s->ref_start + 1;
			}

			MD5_Init(&md5);
			if (start + len > s->ref_end - s->ref_start + 1)
				len = s->ref_end - s->ref_start + 1 - start;
			if (len >= 0)
				MD5_Update(&md5, s->ref + start, len);
			MD5_Final(digest, &md5);
		}
		else if (!s->ref && s->hdr->ref_base_id >= 0) {
			cram_block *b;
			if (s->block_by_id && (b = s->block_by_id[s->hdr->ref_base_id])) {
				MD5_Init(&md5);
				MD5_Update(&md5, b->data, b->uncomp_size);
				MD5_Final(digest, &md5);
			}
		}

		if ((!s->ref && s->hdr->ref_base_id < 0)
			|| memcmp(digest, s->hdr->md5, 16) != 0) {
			char M[33];
			fprintf(stderr, "ERROR: md5sum reference mismatch for ref "
				"%d pos %d..%d\n", ref_id, s->ref_start, s->ref_end);
			fprintf(stderr, "CRAM: %s\n", md5_print(s->hdr->md5, M));
			fprintf(stderr, "Ref : %s\n", md5_print(digest, M));
			return -1;
		}
	}

	if (ref_id == -2) {
		pthread_mutex_lock(&fd->ref_lock);
		pthread_mutex_lock(&fd->refs->lock);

		refs = (char**)calloc(fd->refs->nref, sizeof(char *));
		pthread_mutex_unlock(&fd->refs->lock);
		pthread_mutex_unlock(&fd->ref_lock);
		if (!refs)
			return -1;
	}

	for (rec = 0; rec < s->hdr->num_records; rec++) {
		cram_record *cr = &s->crecs[rec];

		//fprintf(stderr, "Decode seq %d, %d/%d\n", rec, blk->byte, blk->bit);

		cr->s = s;

		out_sz = 1; /* decode 1 item */
		if (ds & CRAM_BF) {
			if (!c->comp_hdr->codecs[DS_BF]) return -1;
			r |= c->comp_hdr->codecs[DS_BF]
				->decode(s, c->comp_hdr->codecs[DS_BF], blk,
				(char *)&bf, &out_sz);
			if (bf < 0 ||
				bf >= sizeof(fd->bam_flag_swap) / sizeof(*fd->bam_flag_swap))
				return -1;
			bf = fd->bam_flag_swap[bf];
			cr->flags = bf;
		}
		else {
			cr->flags = bf = 0x4; // unmapped
		}

		if (ds & CRAM_CF) {
			if (CRAM_MAJOR_VERS(fd->version) == 1) {
				/* CF is byte in 1.0, int32 in 2.0 */
				if (!c->comp_hdr->codecs[DS_CF]) return -1;
				r |= c->comp_hdr->codecs[DS_CF]
					->decode(s, c->comp_hdr->codecs[DS_CF], blk,
					(char *)&cf, &out_sz);
				cr->cram_flags = cf;
			}
			else {
				if (!c->comp_hdr->codecs[DS_CF]) return -1;
				r |= c->comp_hdr->codecs[DS_CF]
					->decode(s, c->comp_hdr->codecs[DS_CF], blk,
					(char *)&cr->cram_flags,
					&out_sz);
				cf = cr->cram_flags;
			}
		}

		if (CRAM_MAJOR_VERS(fd->version) != 1 && ref_id == -2) {
			if (ds & CRAM_RI) {
				if (!c->comp_hdr->codecs[DS_RI]) return -1;
				r |= c->comp_hdr->codecs[DS_RI]
					->decode(s, c->comp_hdr->codecs[DS_RI], blk,
					(char *)&cr->ref_id, &out_sz);
				if ((fd->required_fields & (SAM_SEQ | SAM_TLEN))
					&& cr->ref_id >= 0) {
					if (!fd->no_ref) {
						if (!refs[cr->ref_id])
							refs[cr->ref_id] = cram_get_ref(fd, cr->ref_id,
							1, 0);
						s->ref = refs[cr->ref_id];
					}
					s->ref_start = 1;
					pthread_mutex_lock(&fd->ref_lock);
					pthread_mutex_lock(&fd->refs->lock);
					s->ref_end = fd->refs->ref_id[cr->ref_id]->length;
					pthread_mutex_unlock(&fd->refs->lock);
					pthread_mutex_unlock(&fd->ref_lock);
				}
			}
			else {
				cr->ref_id = 0;
			}
		}
		else {
			cr->ref_id = ref_id; // Forced constant in CRAM 1.0
		}


		if (ds & CRAM_RL) {
			if (!c->comp_hdr->codecs[DS_RL]) return -1;
			r |= c->comp_hdr->codecs[DS_RL]
				->decode(s, c->comp_hdr->codecs[DS_RL], blk,
				(char *)&cr->len, &out_sz);
		}

		if (ds & CRAM_AP) {
			if (!c->comp_hdr->codecs[DS_AP]) return -1;
			r |= c->comp_hdr->codecs[DS_AP]
				->decode(s, c->comp_hdr->codecs[DS_AP], blk,
				(char *)&cr->apos, &out_sz);
			if (c->comp_hdr->AP_delta)
				cr->apos += s->last_apos;
			s->last_apos = cr->apos;
		}
		else {
			cr->apos = c->ref_seq_start;
		}

		if (ds & CRAM_RG) {
			if (!c->comp_hdr->codecs[DS_RG]) return -1;
			r |= c->comp_hdr->codecs[DS_RG]
				->decode(s, c->comp_hdr->codecs[DS_RG], blk,
				(char *)&cr->rg, &out_sz);
			if (cr->rg == unknown_rg)
				cr->rg = -1;
		}
		else {
			cr->rg = -1;
		}

		cr->name_len = 0;

		if (c->comp_hdr->read_names_included) {
			int32_t out_sz2 = 1;

			// Read directly into name cram_block
			cr->name = BLOCK_SIZE(s->name_blk);
			if (ds & CRAM_RN) {
				if (!c->comp_hdr->codecs[DS_RN]) return -1;
				r |= c->comp_hdr->codecs[DS_RN]
					->decode(s, c->comp_hdr->codecs[DS_RN], blk,
					(char *)s->name_blk, &out_sz2);
				cr->name_len = out_sz2;
			}
		}

		cr->mate_pos = 0;
		cr->mate_line = -1;
		cr->mate_ref_id = -1;
		if ((ds & CRAM_CF) && (cf & CRAM_FLAG_DETACHED)) {
			if (ds & CRAM_MF) {
				if (CRAM_MAJOR_VERS(fd->version) == 1) {
					/* MF is byte in 1.0, int32 in 2.0 */
					unsigned char mf;
					if (!c->comp_hdr->codecs[DS_MF]) return -1;
					r |= c->comp_hdr->codecs[DS_MF]
						->decode(s, c->comp_hdr->codecs[DS_MF],
						blk, (char *)&mf, &out_sz);
					cr->mate_flags = mf;
				}
				else {
					if (!c->comp_hdr->codecs[DS_MF]) return -1;
					r |= c->comp_hdr->codecs[DS_MF]
						->decode(s, c->comp_hdr->codecs[DS_MF],
						blk,
						(char *)&cr->mate_flags,
						&out_sz);
				}
			}
			else {
				cr->mate_flags = 0;
			}

			if (!c->comp_hdr->read_names_included) {
				int32_t out_sz2 = 1;

				// Read directly into name cram_block
				cr->name = BLOCK_SIZE(s->name_blk);
				if (ds & CRAM_RN) {
					if (!c->comp_hdr->codecs[DS_RN]) return -1;
					r |= c->comp_hdr->codecs[DS_RN]
						->decode(s, c->comp_hdr->codecs[DS_RN],
						blk, (char *)s->name_blk,
						&out_sz2);
					cr->name_len = out_sz2;
				}
			}

			if (ds & CRAM_NS) {
				if (!c->comp_hdr->codecs[DS_NS]) return -1;
				r |= c->comp_hdr->codecs[DS_NS]
					->decode(s, c->comp_hdr->codecs[DS_NS], blk,
					(char *)&cr->mate_ref_id, &out_sz);
			}

			// Skip as mate_ref of "*" is legit. It doesn't mean unmapped, just unknown.
			//	    if (cr->mate_ref_id == -1 && cr->flags & 0x01) {
			//		/* Paired, but unmapped */
			//		cr->flags |= BAM_FMUNMAP;
			//	    }

			if (ds & CRAM_NP) {
				if (!c->comp_hdr->codecs[DS_NP]) return -1;
				r |= c->comp_hdr->codecs[DS_NP]
					->decode(s, c->comp_hdr->codecs[DS_NP], blk,
					(char *)&cr->mate_pos, &out_sz);
			}

			if (ds & CRAM_TS) {
				if (!c->comp_hdr->codecs[DS_TS]) return -1;
				r |= c->comp_hdr->codecs[DS_TS]
					->decode(s, c->comp_hdr->codecs[DS_TS], blk,
					(char *)&cr->tlen, &out_sz);
			}
			else {
				cr->tlen = INT_MIN;
			}
		}
		else if ((ds & CRAM_CF) && (cf & CRAM_FLAG_MATE_DOWNSTREAM)) {
			if (ds & CRAM_NF) {
				if (!c->comp_hdr->codecs[DS_NF]) return -1;
				r |= c->comp_hdr->codecs[DS_NF]
					->decode(s, c->comp_hdr->codecs[DS_NF], blk,
					(char *)&cr->mate_line, &out_sz);
				cr->mate_line += rec + 1;

				//cr->name_len = sprintf(name, "%d", name_id++);
				//cr->name = DSTRING_LEN(name_ds);
				//dstring_nappend(name_ds, name, cr->name_len);

				cr->mate_ref_id = -1;
				cr->tlen = INT_MIN;
				cr->mate_pos = 0;
			}
			else  {
				cr->mate_flags = 0;
				cr->tlen = INT_MIN;
			}
		}
		else {
			cr->mate_flags = 0;
			cr->tlen = INT_MIN;
		}
		/*
		else if (!name[0]) {
		//name[0] = '?'; name[1] = 0;
		//cr->name_len = 1;
		//cr->name=  DSTRING_LEN(s->name_ds);
		//dstring_nappend(s->name_ds, "?", 1);

		cr->mate_ref_id = -1;
		cr->tlen = 0;
		cr->mate_pos = 0;
		}
		*/

		/* Auxiliary tags */
		if (CRAM_MAJOR_VERS(fd->version) == 1)
			r |= cram_decode_aux_1_0(c, s, blk, cr);
		else
			r |= cram_decode_aux(c, s, blk, cr);

		/* Fake up dynamic string growth and appending */
		if (ds & CRAM_RL) {
			cr->seq = BLOCK_SIZE(s->seqs_blk);
			BLOCK_GROW(s->seqs_blk, cr->len);
			seq = (char *)BLOCK_END(s->seqs_blk);
			BLOCK_SIZE(s->seqs_blk) += cr->len;

			if (!seq)
				return -1;

			cr->qual = BLOCK_SIZE(s->qual_blk);
			BLOCK_GROW(s->qual_blk, cr->len);
			qual = (char *)BLOCK_END(s->qual_blk);
			BLOCK_SIZE(s->qual_blk) += cr->len;

			if (!s->ref)
				memset(seq, '=', cr->len);
		}

		if (!(bf & BAM_FUNMAP)) {
			/* Decode sequence and generate CIGAR */
			if (ds & (CRAM_SEQ | CRAM_MQ)) {
				r |= cram_decode_seq(fd, c, s, blk, cr, bfd, cf, seq, qual);
			}
			else {
				cr->cigar = 0;
				cr->ncigar = 0;
				cr->aend = cr->apos;
				cr->mqual = 0;
			}
		}
		else {
			int out_sz2 = cr->len;

			//puts("Unmapped");
			cr->cigar = 0;
			cr->ncigar = 0;
			cr->aend = cr->apos;
			cr->mqual = 0;

			if (ds & CRAM_BA) {
				if (!c->comp_hdr->codecs[DS_BA]) return -1;
				r |= c->comp_hdr->codecs[DS_BA]
					->decode(s, c->comp_hdr->codecs[DS_BA], blk,
					(char *)seq, &out_sz2);
			}

			if ((ds & CRAM_CF) && (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES)) {
				out_sz2 = cr->len;
				if (ds & CRAM_QS) {
					if (!c->comp_hdr->codecs[DS_QS]) return -1;
					r |= c->comp_hdr->codecs[DS_QS]
						->decode(s, c->comp_hdr->codecs[DS_QS],
						blk, qual, &out_sz2);
				}
			}
			else {
				if (ds & CRAM_RL)
					memset(qual, 30, cr->len);
			}
		}
	}
	pthread_mutex_lock(&fd->ref_lock);
	if (refs) {
		int i;
		for (i = 0; i < fd->refs->nref; i++) {
			if (refs[i])
				cram_ref_decr(fd->refs, i);
		}
		free(refs);
	}
	else if (ref_id >= 0 && s->ref != fd->ref_free) {
		cram_ref_decr(fd->refs, ref_id);
	}
	pthread_mutex_unlock(&fd->ref_lock);
	/* Resolve mate pair cross-references between recs within this slice */
	cram_decode_slice_xref(s, fd->required_fields);

	return r;
}

/*
* Here be dragons! The multi-threading code in this is crufty beyond belief.
*/
static cram_slice *cram_next_slice(cram_fd *fd, cram_container **cp) {
	cram_container *c;
	cram_slice *s = NULL;

	if (!(c = fd->ctr)) {
		// Load first container.
		do {
			if (!(c = fd->ctr = cram_read_container(fd)))
				return NULL;
		} while (c->length == 0);

		/*
		* The first container may be a result of a sub-range query.
		* In which case it may still not be the optimal starting point
		* due to skipped containers/slices in the index.
		*/
		if (fd->range.refid != -2) {
			while (c->ref_seq_id != -2 &&
				(c->ref_seq_id < fd->range.refid ||
				c->ref_seq_start + c->ref_seq_span - 1 < fd->range.start)) {
				if (0 != cram_seek(fd, c->length, SEEK_CUR))
					return NULL;
				cram_free_container(fd->ctr);
				do {
					if (!(c = fd->ctr = cram_read_container(fd)))
						return NULL;
				} while (c->length == 0);
			}

			if (c->ref_seq_id != -2 && c->ref_seq_id != fd->range.refid)
				return NULL;
		}

		if (!(c->comp_hdr_block = cram_read_block(fd)))
			return NULL;
		if (c->comp_hdr_block->content_type != COMPRESSION_HEADER)
			return NULL;

		c->comp_hdr = cram_decode_compression_header(fd, c->comp_hdr_block);
		if (!c->comp_hdr)
			return NULL;
		if (!c->comp_hdr->AP_delta) {
			pthread_mutex_lock(&fd->ref_lock);
			fd->unsorted = 1;
			pthread_mutex_unlock(&fd->ref_lock);
		}
	}

	if ((s = c->slice)) {
		c->slice = NULL;
		cram_free_slice(s);
		s = NULL;
	}

	if (c->curr_slice == c->max_slice) {
		cram_free_container(c);
		c = NULL;
	}

	/* Sorry this is so contorted! */
	for (;;) {
		if (fd->job_pending) {
			cram_decode_job *j = (cram_decode_job *)fd->job_pending;
			c = j->c;
			s = j->s;
			free(fd->job_pending);
			fd->job_pending = NULL;
		}
		else if (!fd->ooc) {
		empty_container:
			if (!c || c->curr_slice == c->max_slice) {
				// new container
				do {
					if (!(c = fd->ctr = cram_read_container(fd))) {
						if (fd->pool) {
							fd->ooc = 1;
							break;
						}

						return NULL;
					}
				} while (c->length == 0);
				if (fd->ooc)
					break;

				/* Skip containers not yet spanning our range */
				if (fd->range.refid != -2 && c->ref_seq_id != -2) {
					fd->required_fields |= SAM_POS;

					if (c->ref_seq_id != fd->range.refid) {
						cram_free_container(c);
						fd->ctr = NULL;
						fd->ooc = 1;
						fd->eof = 1;
						break;
					}

					if (c->ref_seq_start > fd->range.end) {
						cram_free_container(c);
						fd->ctr = NULL;
						fd->ooc = 1;
						fd->eof = 1;
						break;
					}

					if (c->ref_seq_start + c->ref_seq_span - 1 <
						fd->range.start) {
						c->curr_rec = c->max_rec;
						c->curr_slice = c->max_slice;
						cram_seek(fd, c->length, SEEK_CUR);
						cram_free_container(c);
						c = NULL;
						continue;
					}
				}

				if (!(c->comp_hdr_block = cram_read_block(fd)))
					return NULL;
				if (c->comp_hdr_block->content_type != COMPRESSION_HEADER)
					return NULL;

				c->comp_hdr =
					cram_decode_compression_header(fd, c->comp_hdr_block);
				if (!c->comp_hdr)
					return NULL;

				if (!c->comp_hdr->AP_delta) {
					pthread_mutex_lock(&fd->ref_lock);
					fd->unsorted = 1;
					pthread_mutex_unlock(&fd->ref_lock);
				}
			}

			if (c->num_records == 0) {
				cram_free_container(c); c = NULL;
				goto empty_container;
			}


			if (!(s = c->slice = cram_read_slice(fd)))
				return NULL;
			c->curr_slice++;
			c->curr_rec = 0;
			c->max_rec = s->hdr->num_records;

			s->last_apos = s->hdr->ref_seq_start;

			/* Skip slices not yet spanning our range */
			if (fd->range.refid != -2 && s->hdr->ref_seq_id != -2) {
				if (s->hdr->ref_seq_id != fd->range.refid) {
					fd->eof = 1;
					cram_free_slice(s);
					c->slice = NULL;
					return NULL;
				}

				if (s->hdr->ref_seq_start > fd->range.end) {
					fd->eof = 1;
					cram_free_slice(s);
					c->slice = NULL;
					return NULL;
				}

				if (s->hdr->ref_seq_start + s->hdr->ref_seq_span - 1 <
					fd->range.start) {
					cram_free_slice(s);
					c->slice = NULL;
					cram_free_container(c);
					c = NULL;
					continue;
				}
			}
		}

		/* Test decoding of 1st seq */
		if (!c || !s)
			break;

		if (cram_decode_slice_mt(fd, c, s, fd->header) != 0) {
			//	if (cram_decode_slice(fd, c, s, fd->header) != 0) {
			fprintf(stderr, "Failure to decode slice\n");
			cram_free_slice(s);
			c->slice = NULL;
			return NULL;
		}

		if (!fd->pool || fd->job_pending)
			break;

		// Push it a bit far, to qsize in queue rather than pending arrival,
		// as cram tends to be a bit bursty in decode timings.
		if (t_pool_results_queue_len(fd->rqueue) > fd->pool->qsize)
			break;
	}

	if (fd->pool) {
		t_pool_result *res;
		cram_decode_job *j;

		//	fprintf(stderr, "Thread pool len = %d, %d\n",
		//		t_pool_results_queue_len(fd->rqueue),
		//		t_pool_results_queue_sz(fd->rqueue));

		if (fd->ooc && t_pool_results_queue_empty(fd->rqueue))
			return NULL;

		res = t_pool_next_result_wait(fd->rqueue);

		if (!res || !res->data) {
			fprintf(stderr, "t_pool_next_result failure\n");
			return NULL;
		}

		j = (cram_decode_job *)res->data;
		c = j->c;
		s = j->s;

		fd->ctr = c;

		t_pool_delete_result(res, 1);
	}

	*cp = c;
	return s;
}

#else // ~ !defined(_WIN32) && defined(USE_PTHREAD)

/*
 * Decode an entire slice from container blocks. Fills out s->crecs[] array.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_decode_slice(cram_fd *fd, cram_container *c, cram_slice *s,
		      SAM_hdr *bfd) {
    cram_block *blk = s->block[0];
    int32_t bf, ref_id;
    unsigned char cf;
    int out_sz, r = 0;
    int rec;
    char *seq = NULL, *qual = NULL;
    int unknown_rg = -1;
    int embed_ref;
    char **refs = NULL;
    uint32_t ds;

    if (cram_dependent_data_series(fd, c->comp_hdr, s) != 0)
	return -1;

    ds = c->comp_hdr->data_series;

    blk->bit = 7; // MSB first

    /* Look for unknown RG, added as last by Java CRAM? */
    if (bfd->nrg > 0 &&
	!strcmp(bfd->rg[bfd->nrg-1].name, "UNKNOWN"))
	unknown_rg = bfd->nrg-1;

    if (blk->content_type != CORE)
	return -1;

    if (s->crecs)
	free(s->crecs);
    if (!(s->crecs = (cram_record*)malloc(s->hdr->num_records * sizeof(*s->crecs))))
	return -1;

    ref_id = s->hdr->ref_seq_id;
    embed_ref = s->hdr->ref_base_id >= 0 ? 1 : 0;

    if (ref_id >= 0) {
	if (embed_ref) {
	    cram_block *b;
	    if (s->hdr->ref_base_id < 0) {
		fprintf(stderr, "No reference specified and "
			"no embedded reference is available.\n");
		return -1;
	    }
	    if (!s->block_by_id ||
		!(b = s->block_by_id[s->hdr->ref_base_id]))
		return -1;
	    cram_uncompress_block(b);
	    s->ref = (char *)BLOCK_DATA(b);
	    s->ref_start = s->hdr->ref_seq_start;
	    s->ref_end   = s->hdr->ref_seq_start + s->hdr->ref_seq_span-1;
	} else if (!fd->no_ref) {
	    //// Avoid Java cramtools bug by loading entire reference seq 
	    //s->ref = cram_get_ref(fd, s->hdr->ref_seq_id, 1, 0);
	    //s->ref_start = 1;

	    if (fd->required_fields & SAM_SEQ)
		s->ref =
		cram_get_ref(fd, s->hdr->ref_seq_id,
			     s->hdr->ref_seq_start,
			     s->hdr->ref_seq_start + s->hdr->ref_seq_span -1);
	    s->ref_start = s->hdr->ref_seq_start;
	    s->ref_end   = s->hdr->ref_seq_start + s->hdr->ref_seq_span-1;

	    /* Sanity check */
	    if (s->ref_start < 0) {
		fprintf(stderr, "Slice starts before base 1.\n");
		s->ref_start = 0;
	    }

		fd->ref_lock.lock();
		fd->refs->lock.lock();

	    if ((fd->required_fields & SAM_SEQ) &&
		s->ref_end > fd->refs->ref_id[ref_id]->length) {
		fprintf(stderr, "Slice ends beyond reference end.\n");
		s->ref_end = fd->refs->ref_id[ref_id]->length;
	    }

		fd->refs->lock.unlock();
		fd->ref_lock.unlock();
	}
    }

    if ((fd->required_fields & SAM_SEQ) &&
	s->ref == NULL && s->hdr->ref_seq_id >= 0 && !fd->no_ref) {
	fprintf(stderr, "Unable to fetch reference #%d %d..%d\n",
		s->hdr->ref_seq_id, s->hdr->ref_seq_start,
		s->hdr->ref_seq_start + s->hdr->ref_seq_span-1);
	return -1;
    }

    if (CRAM_MAJOR_VERS(fd->version) != 1
	&& (fd->required_fields & SAM_SEQ)
	&& s->hdr->ref_seq_id >= 0
	&& !fd->ignore_md5
	&& memcmp(s->hdr->md5, "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", 16)) {
	MD5_CTX md5;
	unsigned char digest[16];

	if (s->ref && s->hdr->ref_seq_id >= 0) {
	    int start, len;

	    if (s->hdr->ref_seq_start >= s->ref_start) {
		start = s->hdr->ref_seq_start - s->ref_start;
	    } else {
		fprintf(stderr, "Slice starts before base 1.\n");
		start = 0;
	    }

	    if (s->hdr->ref_seq_span <= s->ref_end - s->ref_start + 1) {
		len = s->hdr->ref_seq_span;
	    } else {
		fprintf(stderr, "Slice ends beyond reference end.\n");
		len = s->ref_end - s->ref_start + 1;
	    }

	    MD5_Init(&md5);
	    if (start + len > s->ref_end - s->ref_start + 1)
		len = s->ref_end - s->ref_start + 1 - start;
	    if (len >= 0)
		MD5_Update(&md5, s->ref + start, len);
	    MD5_Final(digest, &md5);
	} else if (!s->ref && s->hdr->ref_base_id >= 0) {
	    cram_block *b;
	    if (s->block_by_id && (b = s->block_by_id[s->hdr->ref_base_id])) {
		MD5_Init(&md5);
		MD5_Update(&md5, b->data, b->uncomp_size);
		MD5_Final(digest, &md5);
	    }
	}

	if ((!s->ref && s->hdr->ref_base_id < 0)
	    || memcmp(digest, s->hdr->md5, 16) != 0) {
	    char M[33];
	    fprintf(stderr, "ERROR: md5sum reference mismatch for ref "
		    "%d pos %d..%d\n", ref_id, s->ref_start, s->ref_end);
	    fprintf(stderr, "CRAM: %s\n", md5_print(s->hdr->md5, M));
	    fprintf(stderr, "Ref : %s\n", md5_print(digest, M));
	    return -1;
	}
    }

    if (ref_id == -2) {

		fd->ref_lock.lock();
		fd->refs->lock.lock();
	
	refs = (char**)calloc(fd->refs->nref, sizeof(char *));

	fd->refs->lock.unlock();
	fd->ref_lock.unlock();

	if (!refs)
	    return -1;
    }

    for (rec = 0; rec < s->hdr->num_records; rec++) {
	cram_record *cr = &s->crecs[rec];

	//fprintf(stderr, "Decode seq %d, %d/%d\n", rec, blk->byte, blk->bit);

	cr->s = s;

	out_sz = 1; /* decode 1 item */
	if (ds & CRAM_BF) {
	    if (!c->comp_hdr->codecs[DS_BF]) return -1;
	    r |= c->comp_hdr->codecs[DS_BF]
		            ->decode(s, c->comp_hdr->codecs[DS_BF], blk,
				     (char *)&bf, &out_sz);
	    if (bf < 0 ||
		bf >= sizeof(fd->bam_flag_swap)/sizeof(*fd->bam_flag_swap))
		return -1;
	    bf = fd->bam_flag_swap[bf];
	    cr->flags = bf;
	} else {
	    cr->flags = bf = 0x4; // unmapped
	}

	if (ds & CRAM_CF) {
	    if (CRAM_MAJOR_VERS(fd->version) == 1) {
		/* CF is byte in 1.0, int32 in 2.0 */
		if (!c->comp_hdr->codecs[DS_CF]) return -1;
		r |= c->comp_hdr->codecs[DS_CF]
		                ->decode(s, c->comp_hdr->codecs[DS_CF], blk,
					 (char *)&cf, &out_sz);
		cr->cram_flags = cf;
	    } else {
		if (!c->comp_hdr->codecs[DS_CF]) return -1;
		r |= c->comp_hdr->codecs[DS_CF]
		                ->decode(s, c->comp_hdr->codecs[DS_CF], blk,
					 (char *)&cr->cram_flags,
					 &out_sz);
		cf = cr->cram_flags;
	    }
	}

	if (CRAM_MAJOR_VERS(fd->version) != 1 && ref_id == -2) {
	    if (ds & CRAM_RI) {
		if (!c->comp_hdr->codecs[DS_RI]) return -1;
		r |= c->comp_hdr->codecs[DS_RI]
		                ->decode(s, c->comp_hdr->codecs[DS_RI], blk,
					 (char *)&cr->ref_id, &out_sz);
		if ((fd->required_fields & (SAM_SEQ|SAM_TLEN))
		    && cr->ref_id >= 0) {
		    if (!fd->no_ref) {
			if (!refs[cr->ref_id])
			    refs[cr->ref_id] = cram_get_ref(fd, cr->ref_id,
							    1, 0);
			s->ref = refs[cr->ref_id];
		    }
		    s->ref_start = 1;

			fd->ref_lock.lock();
			fd->refs->lock.lock();

		    s->ref_end = fd->refs->ref_id[cr->ref_id]->length;

			fd->refs->lock.unlock();
			fd->ref_lock.unlock();

		}
	    } else {
		cr->ref_id = 0;
	    }
	} else {
	    cr->ref_id = ref_id; // Forced constant in CRAM 1.0
	}


	if (ds & CRAM_RL) {
	    if (!c->comp_hdr->codecs[DS_RL]) return -1;
	    r |= c->comp_hdr->codecs[DS_RL]
		            ->decode(s, c->comp_hdr->codecs[DS_RL], blk,
				     (char *)&cr->len, &out_sz);
	}

	if (ds & CRAM_AP) {
	    if (!c->comp_hdr->codecs[DS_AP]) return -1;
	    r |= c->comp_hdr->codecs[DS_AP]
		            ->decode(s, c->comp_hdr->codecs[DS_AP], blk,
				     (char *)&cr->apos, &out_sz);
	    if (c->comp_hdr->AP_delta)
		cr->apos += s->last_apos;
	    s->last_apos=  cr->apos;
	} else {
	    cr->apos = c->ref_seq_start;
	}
		    
	if (ds & CRAM_RG) {
	    if (!c->comp_hdr->codecs[DS_RG]) return -1;
	    r |= c->comp_hdr->codecs[DS_RG]
		           ->decode(s, c->comp_hdr->codecs[DS_RG], blk,
				    (char *)&cr->rg, &out_sz);
	    if (cr->rg == unknown_rg)
		cr->rg = -1;
	} else {
	    cr->rg = -1;
	}

	cr->name_len = 0;

	if (c->comp_hdr->read_names_included) {
	    int32_t out_sz2 = 1;

	    // Read directly into name cram_block
	    cr->name = BLOCK_SIZE(s->name_blk);
	    if (ds & CRAM_RN) {
		if (!c->comp_hdr->codecs[DS_RN]) return -1;
		r |= c->comp_hdr->codecs[DS_RN]
		                ->decode(s, c->comp_hdr->codecs[DS_RN], blk,
					 (char *)s->name_blk, &out_sz2);
		cr->name_len = out_sz2;
	    }
	}

	cr->mate_pos = 0;
	cr->mate_line = -1;
	cr->mate_ref_id = -1;
	if ((ds & CRAM_CF) && (cf & CRAM_FLAG_DETACHED)) {
	    if (ds & CRAM_MF) {
		if (CRAM_MAJOR_VERS(fd->version) == 1) {
		    /* MF is byte in 1.0, int32 in 2.0 */
		    unsigned char mf;
		    if (!c->comp_hdr->codecs[DS_MF]) return -1;
		    r |= c->comp_hdr->codecs[DS_MF]
			            ->decode(s, c->comp_hdr->codecs[DS_MF],
					     blk, (char *)&mf, &out_sz);
		    cr->mate_flags = mf;
		} else {
		    if (!c->comp_hdr->codecs[DS_MF]) return -1;
		    r |= c->comp_hdr->codecs[DS_MF]
			            ->decode(s, c->comp_hdr->codecs[DS_MF],
					     blk,
					     (char *)&cr->mate_flags,
					     &out_sz);
		}
	    } else {
		cr->mate_flags = 0;
	    }

	    if (!c->comp_hdr->read_names_included) {
		int32_t out_sz2 = 1;
	    
		// Read directly into name cram_block
		cr->name = BLOCK_SIZE(s->name_blk);
		if (ds & CRAM_RN) {
		    if (!c->comp_hdr->codecs[DS_RN]) return -1;
		    r |= c->comp_hdr->codecs[DS_RN]
			            ->decode(s, c->comp_hdr->codecs[DS_RN],
					     blk, (char *)s->name_blk,
					     &out_sz2);
		    cr->name_len = out_sz2;
		}
	    }
		    
	    if (ds & CRAM_NS) {
		if (!c->comp_hdr->codecs[DS_NS]) return -1;
		r |= c->comp_hdr->codecs[DS_NS]
		                ->decode(s, c->comp_hdr->codecs[DS_NS], blk,
					 (char *)&cr->mate_ref_id, &out_sz);
	    }

// Skip as mate_ref of "*" is legit. It doesn't mean unmapped, just unknown.
//	    if (cr->mate_ref_id == -1 && cr->flags & 0x01) {
//		/* Paired, but unmapped */
//		cr->flags |= BAM_FMUNMAP;
//	    }

	    if (ds & CRAM_NP) {
		if (!c->comp_hdr->codecs[DS_NP]) return -1;
		r |= c->comp_hdr->codecs[DS_NP]
		                ->decode(s, c->comp_hdr->codecs[DS_NP], blk,
					 (char *)&cr->mate_pos, &out_sz);
	    }

	    if (ds & CRAM_TS) {
		if (!c->comp_hdr->codecs[DS_TS]) return -1;
		r |= c->comp_hdr->codecs[DS_TS]
		                ->decode(s, c->comp_hdr->codecs[DS_TS], blk,
					 (char *)&cr->tlen, &out_sz);
	    } else {
		cr->tlen = INT_MIN;
	    }
	} else if ((ds & CRAM_CF) && (cf & CRAM_FLAG_MATE_DOWNSTREAM)) {
	    if (ds & CRAM_NF) {
		if (!c->comp_hdr->codecs[DS_NF]) return -1;
		r |= c->comp_hdr->codecs[DS_NF]
		                ->decode(s, c->comp_hdr->codecs[DS_NF], blk,
					 (char *)&cr->mate_line, &out_sz);
		cr->mate_line += rec + 1;

		//cr->name_len = sprintf(name, "%d", name_id++);
		//cr->name = DSTRING_LEN(name_ds);
		//dstring_nappend(name_ds, name, cr->name_len);

		cr->mate_ref_id = -1;
		cr->tlen = INT_MIN;
		cr->mate_pos = 0;
	    } else  {
		cr->mate_flags = 0;
		cr->tlen = INT_MIN;
	    }
	} else {
	    cr->mate_flags = 0;
	    cr->tlen = INT_MIN;
	}
	/*
	else if (!name[0]) {
	    //name[0] = '?'; name[1] = 0;
	    //cr->name_len = 1;
	    //cr->name=  DSTRING_LEN(s->name_ds);
	    //dstring_nappend(s->name_ds, "?", 1);

	    cr->mate_ref_id = -1;
	    cr->tlen = 0;
	    cr->mate_pos = 0;
	}
	*/

	/* Auxiliary tags */
	if (CRAM_MAJOR_VERS(fd->version) == 1)
	    r |= cram_decode_aux_1_0(c, s, blk, cr);
	else
	    r |= cram_decode_aux(c, s, blk, cr);

	/* Fake up dynamic string growth and appending */
	if (ds & CRAM_RL) {
	    cr->seq = BLOCK_SIZE(s->seqs_blk);
	    BLOCK_GROW(s->seqs_blk, cr->len);
	    seq = (char *)BLOCK_END(s->seqs_blk);
	    BLOCK_SIZE(s->seqs_blk) += cr->len;

	    if (!seq)
		return -1;
	
	    cr->qual = BLOCK_SIZE(s->qual_blk);
	    BLOCK_GROW(s->qual_blk, cr->len);
	    qual = (char *)BLOCK_END(s->qual_blk);
	    BLOCK_SIZE(s->qual_blk) += cr->len;

	    if (!s->ref)
		memset(seq, '=', cr->len);
	}

	if (!(bf & BAM_FUNMAP)) {
	    /* Decode sequence and generate CIGAR */
	    if (ds & (CRAM_SEQ | CRAM_MQ)) {
		r |= cram_decode_seq(fd, c, s, blk, cr, bfd, cf, seq, qual);
	    } else {
		cr->cigar = 0;
		cr->ncigar = 0;
		cr->aend = cr->apos;
		cr->mqual = 0;
	    }
	} else {
	    int out_sz2 = cr->len;

	    //puts("Unmapped");
	    cr->cigar = 0;
	    cr->ncigar = 0;
	    cr->aend = cr->apos;
	    cr->mqual = 0;

	    if (ds & CRAM_BA) {
		if (!c->comp_hdr->codecs[DS_BA]) return -1;
		r |= c->comp_hdr->codecs[DS_BA]
		                ->decode(s, c->comp_hdr->codecs[DS_BA], blk,
					 (char *)seq, &out_sz2);
	    }

	    if ((ds & CRAM_CF) && (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES)) {
		out_sz2 = cr->len;
		if (ds & CRAM_QS) {
		    if (!c->comp_hdr->codecs[DS_QS]) return -1;
		    r |= c->comp_hdr->codecs[DS_QS]
			            ->decode(s, c->comp_hdr->codecs[DS_QS],
					     blk, qual, &out_sz2);
		}
	    } else {
		if (ds & CRAM_RL)
		    memset(qual, 30, cr->len);
	    }
	}
    }

	fd->ref_lock.lock();

    if (refs) {
	int i;
	for (i = 0; i < fd->refs->nref; i++) {
	    if (refs[i])
		cram_ref_decr(fd->refs, i);
	}
	free(refs);
    } else if (ref_id >= 0 && s->ref != fd->ref_free) {
	cram_ref_decr(fd->refs, ref_id);
    }
	fd->ref_lock.unlock(); 

    /* Resolve mate pair cross-references between recs within this slice */
    cram_decode_slice_xref(s, fd->required_fields);

    return r;
}

/*
 * Here be dragons! The multi-threading code in this is crufty beyond belief.
 */
static cram_slice *cram_next_slice(cram_fd *fd, cram_container **cp) {
    cram_container *c;
    cram_slice *s = NULL;

    if (!(c = fd->ctr)) {
	// Load first container.
	do {
	    if (!(c = fd->ctr = cram_read_container(fd)))
		return NULL;
	} while (c->length == 0);

	/*
	 * The first container may be a result of a sub-range query.
	 * In which case it may still not be the optimal starting point
	 * due to skipped containers/slices in the index. 
	 */
	if (fd->range.refid != -2) {
	    while (c->ref_seq_id != -2 &&
		   (c->ref_seq_id < fd->range.refid ||
		    c->ref_seq_start + c->ref_seq_span-1 < fd->range.start)) {
		if (0 != cram_seek(fd, c->length, SEEK_CUR))
		    return NULL;
		cram_free_container(fd->ctr);
		do {
		    if (!(c = fd->ctr = cram_read_container(fd)))
			return NULL;
		} while (c->length == 0);
	    }

	    if (c->ref_seq_id != -2 && c->ref_seq_id != fd->range.refid)
		return NULL;
	}

	if (!(c->comp_hdr_block = cram_read_block(fd)))
	    return NULL;
	if (c->comp_hdr_block->content_type != COMPRESSION_HEADER)
	    return NULL;

	c->comp_hdr = cram_decode_compression_header(fd, c->comp_hdr_block);
	if (!c->comp_hdr)
	    return NULL;
	if (!c->comp_hdr->AP_delta) {
		fd->ref_lock.lock(); 
	    fd->unsorted = 1;
		fd->ref_lock.unlock(); 
	}
    }

    if ((s = c->slice)) {
	c->slice = NULL;
	cram_free_slice(s);
	s = NULL;
    }

    if (c->curr_slice == c->max_slice) {
	cram_free_container(c);
	c = NULL;
    }

    /* Sorry this is so contorted! */
    for (;;) {
	if (fd->job_pending) {
	    cram_decode_job *j = (cram_decode_job *)fd->job_pending;
	    c = j->c;
	    s = j->s;
	    free(fd->job_pending);
	    fd->job_pending = NULL;
	} else if (!fd->ooc) {
	empty_container:
	    if (!c || c->curr_slice == c->max_slice) {
		// new container
		do {
		    if (!(c = fd->ctr = cram_read_container(fd))) {
			if (fd->pool) {
			    fd->ooc = 1;
			    break;
			}

			return NULL;
		    }
		} while (c->length == 0);
		if (fd->ooc)
		    break;

		/* Skip containers not yet spanning our range */
		if (fd->range.refid != -2 && c->ref_seq_id != -2) {
		    fd->required_fields |= SAM_POS;

		    if (c->ref_seq_id != fd->range.refid) {
			cram_free_container(c);
			fd->ctr = NULL;
			fd->ooc = 1;
			fd->eof = 1;
			break;
		    }

		    if (c->ref_seq_start > fd->range.end) {
			cram_free_container(c);
			fd->ctr = NULL;
			fd->ooc = 1;
			fd->eof = 1;
			break;
		    }

		    if (c->ref_seq_start + c->ref_seq_span-1 <
			fd->range.start) {
			c->curr_rec = c->max_rec;
			c->curr_slice = c->max_slice;
			cram_seek(fd, c->length, SEEK_CUR);
			cram_free_container(c);
			c = NULL;
			continue;
		    }
		}

		if (!(c->comp_hdr_block = cram_read_block(fd)))
		    return NULL;
		if (c->comp_hdr_block->content_type != COMPRESSION_HEADER)
		    return NULL;

		c->comp_hdr =
		    cram_decode_compression_header(fd, c->comp_hdr_block);
		if (!c->comp_hdr)
		    return NULL;

		if (!c->comp_hdr->AP_delta) {
			fd->ref_lock.lock(); 
		    fd->unsorted = 1;
			fd->ref_lock.unlock(); 
		}
	    }

	    if (c->num_records == 0) {
		cram_free_container(c); c = NULL;
		goto empty_container;
	    }


	    if (!(s = c->slice = cram_read_slice(fd)))
		return NULL;
	    c->curr_slice++;
	    c->curr_rec = 0;
	    c->max_rec = s->hdr->num_records;

	    s->last_apos = s->hdr->ref_seq_start;
	    
	    /* Skip slices not yet spanning our range */
	    if (fd->range.refid != -2 && s->hdr->ref_seq_id != -2) {
		if (s->hdr->ref_seq_id != fd->range.refid) {
		    fd->eof = 1;
		    cram_free_slice(s);
		    c->slice = NULL;
		    return NULL;
		}

		if (s->hdr->ref_seq_start > fd->range.end) {
		    fd->eof = 1;
		    cram_free_slice(s);
		    c->slice = NULL;
		    return NULL;
		}

		if (s->hdr->ref_seq_start + s->hdr->ref_seq_span-1 <
		    fd->range.start) {
		    cram_free_slice(s);
		    c->slice = NULL;
		    cram_free_container(c);
		    c = NULL;
		    continue;
		}
	    }
	}

	/* Test decoding of 1st seq */
	if (!c || !s)
	    break;

	if (cram_decode_slice_mt(fd, c, s, fd->header) != 0) {
	    //	if (cram_decode_slice(fd, c, s, fd->header) != 0) {
	    fprintf(stderr, "Failure to decode slice\n");
	    cram_free_slice(s);
	    c->slice = NULL;
	    return NULL;
	}

	if (!fd->pool || fd->job_pending)
	    break;

	// Push it a bit far, to qsize in queue rather than pending arrival,
	// as cram tends to be a bit bursty in decode timings.
	if (t_pool_results_queue_len(fd->rqueue) > fd->pool->qsize)
	    break;
    }

    if (fd->pool) {
	t_pool_result *res;
	cram_decode_job *j;
	
//	fprintf(stderr, "Thread pool len = %d, %d\n",
//		t_pool_results_queue_len(fd->rqueue),
//		t_pool_results_queue_sz(fd->rqueue));

	if (fd->ooc && t_pool_results_queue_empty(fd->rqueue))
	    return NULL;

	res = t_pool_next_result_wait(fd->rqueue);

	if (!res || !res->data) {
	    fprintf(stderr, "t_pool_next_result failure\n");
	    return NULL;
	}

	j = (cram_decode_job *)res->data;
	c = j->c;
	s = j->s;

	fd->ctr = c;

	t_pool_delete_result(res, 1);
    }

    *cp = c;
    return s;
}

#endif // ~ !defined(_WIN32) && defined(USE_PTHREAD)

/*
* Read the next cram record and return it.
* Note that to decode cram_record the caller will need to look up some data
* in the current slice, pointed to by fd->ctr->slice. This is valid until
* the next call to cram_get_seq (which may invalidate it).
*
* Returns record pointer on success (do not free)
*        NULL on failure
*/
cram_record *cram_get_seq(cram_fd *fd) {
	cram_container *c;
	cram_slice *s;

	for (;;) {
		c = fd->ctr;
		if (c && c->slice && c->curr_rec < c->max_rec) {
			s = c->slice;
		}
		else {
			if (!(s = cram_next_slice(fd, &c)))
				return NULL;
		}

		if (fd->range.refid != -2) {
			if (s->crecs[c->curr_rec].ref_id < fd->range.refid) {
				c->curr_rec++;
				continue;
			}

			if (s->crecs[c->curr_rec].ref_id != fd->range.refid) {
				fd->eof = 1;
				cram_free_slice(s);
				c->slice = NULL;
				return NULL;
			}

			if (s->crecs[c->curr_rec].apos > fd->range.end) {
				fd->eof = 1;
				cram_free_slice(s);
				c->slice = NULL;
				return NULL;
			}

			if (s->crecs[c->curr_rec].aend < fd->range.start) {
				c->curr_rec++;
				continue;
			}
		}

		break;
	}

	fd->ctr = c;
	c->slice = s;
	return &s->crecs[c->curr_rec++];
}