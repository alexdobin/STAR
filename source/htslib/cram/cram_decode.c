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

    if (!(b = cram_new_block(0, 0)))
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
    if (!(h->TL = calloc(h->nTL, sizeof(unsigned char *))))
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
    cram_block_compression_hdr *hdr = calloc(1, sizeof(*hdr));
    int i;
    int32_t map_size, map_count;

    if (!hdr)
	return NULL;

    if (b->method != RAW) {
	if (cram_uncompress_block(b))
	    return NULL;
    }

    cp = (char *)b->data;

    if (fd->version == CRAM_1_VERS) {
	cp += itf8_get(cp, &hdr->ref_seq_id);
	cp += itf8_get(cp, &hdr->ref_seq_start);
	cp += itf8_get(cp, &hdr->ref_seq_span);
	cp += itf8_get(cp, &hdr->num_records);
	cp += itf8_get(cp, &hdr->num_landmarks);
	if (!(hdr->landmark = malloc(hdr->num_landmarks * sizeof(int32_t)))) {
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
	cram_map *m = malloc(sizeof(*m)); // FIXME: use pooled_alloc

	if (!m) {
	    cram_free_compression_header(hdr);
	    return NULL;
	}

	cp += 2;
	cp += itf8_get(cp, &encoding);
	cp += itf8_get(cp, &size);

	// Fill out cram_map purely for cram_dump to dump out.
	m->key = (key[0]<<8)|key[1];
	m->encoding = encoding;
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
	    if (!(hdr->BF_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'C' && key[1] == 'F') {
	    if (!(hdr->CF_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'I') {
	    if (!(hdr->RI_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'L') {
	    if (!(hdr->RL_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'A' && key[1] == 'P') {
	    if (!(hdr->AP_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'G') {
	    if (!(hdr->RG_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'M' && key[1] == 'F') {
	    if (!(hdr->MF_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'N' && key[1] == 'S') {
	    if (!(hdr->NS_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'N' && key[1] == 'P') {
	    if (!(hdr->NP_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'S') {
	    if (!(hdr->TS_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'N' && key[1] == 'F') {
	    if (!(hdr->NF_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'C') {
	    if (!(hdr->TC_codec = cram_decoder_init(encoding, cp, size, E_BYTE,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'N') {
	    if (!(hdr->TN_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'F' && key[1] == 'N') {
	    if (!(hdr->FN_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'F' && key[1] == 'C') {
	    if (!(hdr->FC_codec = cram_decoder_init(encoding, cp, size, E_BYTE,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'F' && key[1] == 'P') {
	    if (!(hdr->FP_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'B' && key[1] == 'S') {
	    if (!(hdr->BS_codec = cram_decoder_init(encoding, cp, size, E_BYTE,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'I' && key[1] == 'N') {
	    if (!(hdr->IN_codec = cram_decoder_init(encoding, cp, size,
						    E_BYTE_ARRAY,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'S' && key[1] == 'C') {
	    if (!(hdr->SC_codec = cram_decoder_init(encoding, cp, size,
						    E_BYTE_ARRAY,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'D' && key[1] == 'L') {
	    if (!(hdr->DL_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'B' && key[1] == 'A') {
	    if (!(hdr->BA_codec = cram_decoder_init(encoding, cp, size, E_BYTE,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'S') {
	    if (!(hdr->RS_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'P' && key[1] == 'D') {
	    if (!(hdr->PD_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'H' && key[1] == 'C') {
	    if (!(hdr->HC_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'M' && key[1] == 'Q') {
	    if (!(hdr->MQ_codec = cram_decoder_init(encoding, cp, size, E_INT,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'R' && key[1] == 'N') {
	    if (!(hdr->RN_codec = cram_decoder_init(encoding, cp, size,
						    E_BYTE_ARRAY_BLOCK,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'Q' && key[1] == 'S') {
	    if (!(hdr->QS_codec = cram_decoder_init(encoding, cp, size, E_BYTE,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	    if (!(hdr->Qs_codec = cram_decoder_init(encoding, cp, size,
						    E_BYTE_ARRAY,
						    fd->version))) {
		cram_free_compression_header(hdr);
		return NULL;
	    }
	} else if (key[0] == 'T' && key[1] == 'L') {
	    if (!(hdr->TL_codec = cram_decoder_init(encoding, cp, size, E_INT,
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
	cram_map *m = malloc(sizeof(*m)); // FIXME: use pooled_alloc
	char *key = cp+1;

	if (!m) {
	    cram_free_compression_header(hdr);
	    return NULL;
	}

	m->key = (key[0]<<16)|(key[1]<<8)|key[2];

	cp += 4; // Strictly ITF8, but this suffices
	cp += itf8_get(cp, &encoding);
	cp += itf8_get(cp, &size);

	m->encoding = encoding;
	m->size     = size;
	m->offset   = cp - (char *)b->data;
	if (!(m->codec = cram_decoder_init(encoding, cp, size,
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

    if (!(hdr  = calloc(1, sizeof(*hdr))))
	return NULL;

    hdr->content_type = b->content_type;

    if (b->content_type == MAPPED_SLICE) {
	cp += itf8_get(cp, &hdr->ref_seq_id);
	cp += itf8_get(cp, &hdr->ref_seq_start);
	cp += itf8_get(cp, &hdr->ref_seq_span);
    }
    cp += itf8_get(cp, &hdr->num_records);
    if (fd->version != CRAM_1_VERS)
	cp += itf8_get(cp, &hdr->record_counter);
    cp += itf8_get(cp, &hdr->num_blocks);

    cp += itf8_get(cp, &hdr->num_content_ids);
    hdr->block_content_ids = malloc(hdr->num_content_ids * sizeof(int32_t));
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

    if (fd->version != CRAM_1_VERS) {
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
    enum cigar_op cig_op = BAM_CMATCH;
    uint32_t *cigar = s->cigar;
    uint32_t ncigar = s->ncigar;
    uint32_t cigar_alloc = s->cigar_alloc;
    uint32_t nm = 0, md_dist = 0;
    int orig_aux = 0;
    int decode_md = fd->decode_md;
    char buf[20];

    if (!(cf & CRAM_FLAG_PRESERVE_QUAL_SCORES)) {
	memset(qual, 30, cr->len);
    }

    if (decode_md) {
	orig_aux = BLOCK_SIZE(s->aux_blk);
	BLOCK_APPEND(s->aux_blk, "MDZ", 3);
    }
    
    if (!c->comp_hdr->FN_codec) return -1;
    r |= c->comp_hdr->FN_codec->decode(s,c->comp_hdr->FN_codec, blk,
				       (char *)&fn, &out_sz);

    ref_pos--; // count from 0
    cr->cigar = ncigar;
    for (f = 0; f < fn; f++) {
	int32_t pos;
	char op;

	if (ncigar+2 >= cigar_alloc) {
	    cigar_alloc = cigar_alloc ? cigar_alloc*2 : 1024;
	    s->cigar = cigar;
	    if (!(cigar = realloc(cigar, cigar_alloc * sizeof(*cigar))))
		return -1;
	}

	if (!c->comp_hdr->FC_codec) return -1;
	r |= c->comp_hdr->FC_codec->decode(s, c->comp_hdr->FC_codec, blk,
					   &op,  &out_sz);
	if (!c->comp_hdr->FP_codec) return -1;
	r |= c->comp_hdr->FP_codec->decode(s, c->comp_hdr->FP_codec, blk,
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
	    cig_op = BAM_CMATCH;
#endif
	    cig_len += pos - seq_pos;
	    ref_pos += pos - seq_pos;
	    md_dist += pos - seq_pos;
	    seq_pos = pos;
	}

	prev_pos = pos;

	switch(op) {
	case 'S': { // soft clip: IN
	    int32_t out_sz2 = 1;

	    if (cig_len) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (fd->version == CRAM_1_VERS) {
		r |= c->comp_hdr->IN_codec
		    ? c->comp_hdr->IN_codec->decode(s, c->comp_hdr->IN_codec,
						    blk, &seq[pos-1], &out_sz2)
		    : (seq[pos-1] = 'N', out_sz2 = 1, 0);
	    } else {
		r |= c->comp_hdr->SC_codec
		    ? c->comp_hdr->SC_codec->decode(s, c->comp_hdr->SC_codec,
						    blk, &seq[pos-1], &out_sz2)
		    : (seq[pos-1] = 'N', out_sz2 = 1, 0);
	    }
	    cigar[ncigar++] = (out_sz2<<4) + BAM_CSOFT_CLIP;
	    cig_op = BAM_CSOFT_CLIP;
	    seq_pos += out_sz2;
	    break;
	}

	case 'X': { // Substitution; BS
	    unsigned char base;
#ifdef USE_X
	    if (cig_len && cig_op != BAM_CBASE_MISMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (!c->comp_hdr->BS_codec) return -1;
	    r |= c->comp_hdr->BS_codec->decode(s, c->comp_hdr->BS_codec, blk,
					       (char *)&base, &out_sz);
	    seq[pos-1] = 'N'; // FIXME look up BS=base value
	    cig_op = BAM_CBASE_MISMATCH;
#else
	    int ref_base;
	    if (cig_len && cig_op != BAM_CMATCH) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (!c->comp_hdr->BS_codec) return -1;
	    r |= c->comp_hdr->BS_codec->decode(s, c->comp_hdr->BS_codec, blk,
					       (char *)&base, &out_sz);
	    if (ref_pos >= bfd->ref[cr->ref_id].len || !s->ref) {
		seq[pos-1] = 'N';
	    } else {
		ref_base = fd->L1[(uc)s->ref[ref_pos - s->ref_start +1]];
		seq[pos-1] = c->comp_hdr->substitution_matrix[ref_base][base];
		if (decode_md) {
		    BLOCK_APPENDF_2(s->aux_blk, buf, "%d%c",
				    md_dist, s->ref[ref_pos-s->ref_start +1]);
		    md_dist = 0;
		}
	    }
	    cig_op = BAM_CMATCH;
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
	    if (!c->comp_hdr->DL_codec) return -1;
	    r |= c->comp_hdr->DL_codec->decode(s, c->comp_hdr->DL_codec, blk,
					       (char *)&i32, &out_sz);
	    if (decode_md) {
		BLOCK_APPENDF_1(s->aux_blk, buf, "%d^", md_dist);
		BLOCK_APPEND(s->aux_blk, &s->ref[ref_pos - s->ref_start +1],
			     i32);
		md_dist = 0;
	    }
	    cig_op = BAM_CDEL;
	    cig_len += i32;
	    ref_pos += i32;
	    nm      += i32;
	    //printf("  %d: DL = %d (ret %d)\n", f, i32, r);
	    break;
	}

	case 'I': { // Insertion (several bases); IN
	    int32_t out_sz2 = 1;

	    if (cig_len && cig_op != BAM_CINS) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }

	    if (!c->comp_hdr->IN_codec) return -1;
	    r |= c->comp_hdr->IN_codec->decode(s, c->comp_hdr->IN_codec, blk,
					       &seq[pos-1], &out_sz2);
	    cig_op = BAM_CINS;
	    cig_len += out_sz2;
	    seq_pos += out_sz2;
	    nm      += out_sz2;
	    //printf("  %d: IN(I) = %.*s (ret %d, out_sz %d)\n", f, out_sz2, dat, r, out_sz2);
	    break;
	}

	case 'i': { // Insertion (single base); BA
	    if (cig_len && cig_op != BAM_CINS) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (!c->comp_hdr->BA_codec) return -1;
	    r |= c->comp_hdr->BA_codec->decode(s, c->comp_hdr->BA_codec, blk,
					       (char *)&seq[pos-1], &out_sz);
	    cig_op = BAM_CINS;
	    cig_len++;
	    seq_pos++;
	    nm++;
	    //printf("  %d: BA = %c (ret %d)\n", f, seq[pos-1], r);
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
	    if (!c->comp_hdr->BA_codec) return -1;
	    r |= c->comp_hdr->BA_codec->decode(s, c->comp_hdr->BA_codec, blk,
					       (char *)&seq[pos-1], &out_sz);
	    if (!c->comp_hdr->QS_codec) return -1;
	    r |= c->comp_hdr->QS_codec->decode(s, c->comp_hdr->QS_codec, blk,
					       (char *)&qual[pos-1], &out_sz);
#ifdef USE_X
	    cig_op = BAM_CBASE_MISMATCH;
#else
	    cig_op = BAM_CMATCH;
#endif
	    cig_len++;
	    seq_pos++;
	    ref_pos++;
	    //printf("  %d: BA/QS(B) = %c/%d (ret %d)\n", f, i32, qc, r);
	    break;
	}

	case 'Q': { // Quality score; QS
	    if (!c->comp_hdr->QS_codec) return -1;
	    r |= c->comp_hdr->QS_codec->decode(s, c->comp_hdr->QS_codec, blk,
					       (char *)&qual[pos-1], &out_sz);
	    //printf("  %d: QS = %d (ret %d)\n", f, qc, r);
	    break;
	}

	case 'H': { // hard clip; HC
	    if (cig_len && cig_op != BAM_CHARD_CLIP) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (!c->comp_hdr->HC_codec) return -1;
	    r |= c->comp_hdr->HC_codec->decode(s, c->comp_hdr->HC_codec, blk,
					       (char *)&i32, &out_sz);
	    cig_op = BAM_CHARD_CLIP;
	    cig_len += i32;
	    nm      += i32;
	    break;
	}

	case 'P': { // padding; PD
	    if (cig_len && cig_op != BAM_CPAD) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (!c->comp_hdr->PD_codec) return -1;
	    r |= c->comp_hdr->PD_codec->decode(s, c->comp_hdr->PD_codec, blk,
					       (char *)&i32, &out_sz);
	    cig_op = BAM_CPAD;
	    cig_len += i32;
	    nm      += i32;
	    break;
	}

	case 'N': { // Ref skip; RS
	    if (cig_len && cig_op != BAM_CREF_SKIP) {
		cigar[ncigar++] = (cig_len<<4) + cig_op;
		cig_len = 0;
	    }
	    if (!c->comp_hdr->RS_codec) return -1;
	    r |= c->comp_hdr->RS_codec->decode(s, c->comp_hdr->RS_codec, blk,
					       (char *)&i32, &out_sz);
	    cig_op = BAM_CREF_SKIP;
	    cig_len += i32;
	    ref_pos += i32;
	    nm      += i32;
	    break;
	}

	default:
	    abort();
	}
    }

    /* An implement match op for any unaccounted for bases */
    if (cr->len >= seq_pos) {
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
	    if (!(cigar = realloc(cigar, cigar_alloc * sizeof(*cigar))))
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
	cig_op = BAM_CMATCH;
#endif
	cig_len += cr->len - seq_pos+1;
    }
    if (decode_md) {
	BLOCK_APPENDF_1(s->aux_blk, buf, "%d", md_dist);
    }

    if (cig_len) {
	if (ncigar >= cigar_alloc) {
	    cigar_alloc = cigar_alloc ? cigar_alloc*2 : 1024;
	    s->cigar = cigar;
	    if (!(cigar = realloc(cigar, cigar_alloc * sizeof(*cigar))))
		return -1;
	}

	cigar[ncigar++] = (cig_len<<4) + cig_op;
    }

    cr->ncigar = ncigar - cr->cigar;
    cr->aend = ref_pos;

    //printf("2: %.*s %d .. %d\n", cr->name_len, DSTRING_STR(name_ds) + cr->name, cr->apos, ref_pos);

    if (!c->comp_hdr->MQ_codec) return -1;
    r |= c->comp_hdr->MQ_codec->decode(s, c->comp_hdr->MQ_codec, blk,
				       (char *)&cr->mqual, &out_sz);

    if (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES) {
	int32_t out_sz2 = cr->len;

	if (!c->comp_hdr->Qs_codec) return -1;
	r |= c->comp_hdr->Qs_codec->decode(s, c->comp_hdr->Qs_codec, blk,
					   qual, &out_sz2);
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
	    
    if (!c->comp_hdr->TC_codec) return -1;
    r |= c->comp_hdr->TC_codec->decode(s, c->comp_hdr->TC_codec, blk,
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
	if (!c->comp_hdr->TN_codec) return -1;
	r |= c->comp_hdr->TN_codec->decode(s, c->comp_hdr->TN_codec,
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
    int32_t TL;
    unsigned char *TN;
	    
    if (!c->comp_hdr->TL_codec) return -1;
    r |= c->comp_hdr->TL_codec->decode(s, c->comp_hdr->TL_codec, blk,
				       (char *)&TL, &out_sz);
    if (r || TL < 0 || TL >= c->comp_hdr->nTL)
	return -1;

    TN = c->comp_hdr->TL[TL];
    cr->ntags = strlen((char *)TN)/3; // optimise to remove strlen

    //printf("TC=%d\n", cr->ntags);
    cr->aux_size = 0;
    cr->aux = BLOCK_SIZE(s->aux_blk);

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
static void cram_decode_slice_xref(cram_slice *s) {
    int rec;

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

		    do {
			if (aleft > s->crecs[id2].apos)
			    aleft = s->crecs[id2].apos;
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
			    if (s->crecs[id2].aend != aright)
				s->crecs[id2].tlen = tlen;
			    else if (s->crecs[id2].flags & BAM_FREAD1)
				s->crecs[id2].tlen = tlen;
			    else
				s->crecs[id2].tlen = -tlen;
			} else {
			    s->crecs[id2].tlen = -tlen;
			}

			id2 = s->crecs[id2].mate_line;
			while (id2 != id1) {
			    if (s->crecs[id2].apos == aleft) {
				if (s->crecs[id2].aend != aright)
				    s->crecs[id2].tlen = tlen;
				else if (s->crecs[id2].flags & BAM_FREAD1)
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
    char *seq, *qual;
    int unknown_rg = -1;
    int id, embed_ref;
    char **refs = NULL;

    for (id = 0; id < s->hdr->num_blocks; id++) {
	if (cram_uncompress_block(s->block[id]))
	    return -1;
    }

    blk->bit = 7; // MSB first

    /* Look for unknown RG, added as last by Java CRAM? */
    if (bfd->nrg > 0 &&
	!strcmp(bfd->rg[bfd->nrg-1].name, "UNKNOWN"))
	unknown_rg = bfd->nrg-1;

    if (blk->content_type != CORE)
	return -1;

    if (s->crecs)
	free(s->crecs);
    if (!(s->crecs = malloc(s->hdr->num_records * sizeof(*s->crecs))))
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
	    s->ref = (char *)BLOCK_DATA(b);
	    s->ref_start = s->hdr->ref_seq_start;
	    s->ref_end   = s->hdr->ref_seq_start + s->hdr->ref_seq_span-1;
	} else if (!fd->no_ref) {
	    //// Avoid Java cramtools bug by loading entire reference seq 
	    //s->ref = cram_get_ref(fd, s->hdr->ref_seq_id, 1, 0);
	    //s->ref_start = 1;

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
	    pthread_mutex_lock(&fd->ref_lock);
	    pthread_mutex_lock(&fd->refs->lock);
	    if (s->ref_end > fd->refs->ref_id[ref_id]->length) {
		fprintf(stderr, "Slice ends beyond reference end.\n");
		s->ref_end = fd->refs->ref_id[ref_id]->length;
	    }
	    pthread_mutex_unlock(&fd->refs->lock);
	    pthread_mutex_unlock(&fd->ref_lock);
	}
    }

    if (s->ref == NULL && s->hdr->ref_seq_id >= 0 && !fd->no_ref) {
	fprintf(stderr, "Unable to fetch reference #%d %d..%d\n",
		s->hdr->ref_seq_id, s->hdr->ref_seq_start,
		s->hdr->ref_seq_start + s->hdr->ref_seq_span-1);
	return -1;
    }

    if (fd->version != CRAM_1_VERS && s->hdr->ref_seq_id >= 0
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
	pthread_mutex_lock(&fd->ref_lock);
	pthread_mutex_lock(&fd->refs->lock);
	refs = calloc(fd->refs->nref, sizeof(char *));
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
	if (!c->comp_hdr->BF_codec) return -1;
	r |= c->comp_hdr->BF_codec->decode(s, c->comp_hdr->BF_codec, blk,
					   (char *)&bf, &out_sz);
	if (bf < 0 ||
	    bf >= sizeof(fd->bam_flag_swap)/sizeof(*fd->bam_flag_swap))
	    return -1;
	bf = fd->bam_flag_swap[bf];
	cr->flags = bf;

	if (fd->version == CRAM_1_VERS) {
	    /* CF is byte in 1.0, int32 in 2.0 */
	    if (!c->comp_hdr->CF_codec) return -1;
	    r |= c->comp_hdr->CF_codec->decode(s, c->comp_hdr->CF_codec, blk,
					       (char *)&cf, &out_sz);
	    cr->cram_flags = cf;
	} else {
	    if (!c->comp_hdr->CF_codec) return -1;
	    r |= c->comp_hdr->CF_codec->decode(s, c->comp_hdr->CF_codec, blk,
					       (char *)&cr->cram_flags,
					       &out_sz);
	    cf = cr->cram_flags;
	}

	if (fd->version != CRAM_1_VERS && ref_id == -2) {
	    if (!c->comp_hdr->RI_codec) return -1;
	    r |= c->comp_hdr->RI_codec->decode(s, c->comp_hdr->RI_codec, blk,
					       (char *)&cr->ref_id, &out_sz);
	    if (cr->ref_id >= 0) {
		if (!fd->no_ref) {
		    if (!refs[cr->ref_id])
			refs[cr->ref_id] = cram_get_ref(fd, cr->ref_id, 1, 0);
		    s->ref = refs[cr->ref_id];
		}
		s->ref_start = 1;
		pthread_mutex_lock(&fd->ref_lock);
		pthread_mutex_lock(&fd->refs->lock);
		s->ref_end = fd->refs->ref_id[cr->ref_id]->length;
		pthread_mutex_unlock(&fd->refs->lock);
		pthread_mutex_unlock(&fd->ref_lock);
	    }
	} else {
	    cr->ref_id = ref_id; // Forced constant in CRAM 1.0
	}


	if (!c->comp_hdr->RL_codec) return -1;
	r |= c->comp_hdr->RL_codec->decode(s, c->comp_hdr->RL_codec, blk,
					   (char *)&cr->len, &out_sz);

	if (!c->comp_hdr->AP_codec) return -1;
	r |= c->comp_hdr->AP_codec->decode(s, c->comp_hdr->AP_codec, blk,
					   (char *)&cr->apos, &out_sz);
	if (c->comp_hdr->AP_delta)
	    cr->apos += s->last_apos;
	s->last_apos=  cr->apos;
		    
	if (!c->comp_hdr->RG_codec) return -1;
	r |= c->comp_hdr->RG_codec->decode(s, c->comp_hdr->RG_codec, blk,
					   (char *)&cr->rg, &out_sz);
	if (cr->rg == unknown_rg)
	    cr->rg = -1;

	cr->name_len = 0;

	if (c->comp_hdr->read_names_included) {
	    int32_t out_sz2 = 1;

	    // Read directly into name cram_block
	    cr->name = BLOCK_SIZE(s->name_blk);
	    if (!c->comp_hdr->RN_codec) return -1;
	    r |= c->comp_hdr->RN_codec->decode(s, c->comp_hdr->RN_codec, blk,
					       (char *)s->name_blk, &out_sz2);
	    cr->name_len = out_sz2;
	}

	cr->mate_line = -1;
	cr->mate_ref_id = -1;
	if (cf & CRAM_FLAG_DETACHED) {
	    if (fd->version == CRAM_1_VERS) {
		/* MF is byte in 1.0, int32 in 2.0 */
		unsigned char mf;
		if (!c->comp_hdr->MF_codec) return -1;
		r |= c->comp_hdr->MF_codec->decode(s, c->comp_hdr->MF_codec,
						   blk, (char *)&mf, &out_sz);
		cr->mate_flags = mf;
	    } else {
		if (!c->comp_hdr->MF_codec) return -1;
		r |= c->comp_hdr->MF_codec->decode(s, c->comp_hdr->MF_codec,
						   blk,
						   (char *)&cr->mate_flags,
						   &out_sz);
	    }

	    if (!c->comp_hdr->read_names_included) {
		int32_t out_sz2 = 1;
	    
		// Read directly into name cram_block
		cr->name = BLOCK_SIZE(s->name_blk);
		if (!c->comp_hdr->RN_codec) return -1;
		r |= c->comp_hdr->RN_codec->decode(s, c->comp_hdr->RN_codec,
						   blk, (char *)s->name_blk,
						   &out_sz2);
		cr->name_len = out_sz2;
	    }
		    
	    if (!c->comp_hdr->NS_codec) return -1;
	    r |= c->comp_hdr->NS_codec->decode(s, c->comp_hdr->NS_codec, blk,
					       (char *)&cr->mate_ref_id, &out_sz);

// Skip as mate_ref of "*" is legit. It doesn't mean unmapped, just unknown.
//	    if (cr->mate_ref_id == -1 && cr->flags & 0x01) {
//		/* Paired, but unmapped */
//		cr->flags |= BAM_FMUNMAP;
//	    }

	    if (!c->comp_hdr->NP_codec) return -1;
	    r |= c->comp_hdr->NP_codec->decode(s, c->comp_hdr->NP_codec, blk,
					       (char *)&cr->mate_pos, &out_sz);
	    if (!c->comp_hdr->TS_codec) return -1;
	    r |= c->comp_hdr->TS_codec->decode(s, c->comp_hdr->TS_codec, blk,
					       (char *)&cr->tlen, &out_sz);
	} else if (cf & CRAM_FLAG_MATE_DOWNSTREAM) {
	    if (!c->comp_hdr->NF_codec) return -1;
	    r |= c->comp_hdr->NF_codec->decode(s, c->comp_hdr->NF_codec, blk,
					       (char *)&cr->mate_line, &out_sz);
	    cr->mate_line += rec + 1;

	    //cr->name_len = sprintf(name, "%d", name_id++);
	    //cr->name = DSTRING_LEN(name_ds);
	    //dstring_nappend(name_ds, name, cr->name_len);

	    cr->mate_ref_id = -1;
	    cr->tlen = INT_MIN;
	    cr->mate_pos = 0;
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
	if (fd->version == CRAM_1_VERS)
	    r |= cram_decode_aux_1_0(c, s, blk, cr);
	else
	    r |= cram_decode_aux(c, s, blk, cr);

	/* Fake up dynamic string growth and appending */
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

	if (!(bf & BAM_FUNMAP)) {
	    /* Decode sequence and generate CIGAR */
	    r |= cram_decode_seq(fd, c, s, blk, cr, bfd, cf, seq, qual);
	} else {
	    int out_sz2 = cr->len;

	    //puts("Unmapped");
	    cr->cigar = 0;
	    cr->ncigar = 0;
	    cr->aend = cr->apos;
	    cr->mqual = 0;

	    if (!c->comp_hdr->BA_codec) return -1;
	    r |= c->comp_hdr->BA_codec->decode(s, c->comp_hdr->BA_codec, blk,
					       (char *)seq, &out_sz2);

	    if (cf & CRAM_FLAG_PRESERVE_QUAL_SCORES) {
		out_sz2 = cr->len;
		if (!c->comp_hdr->Qs_codec) return -1;
		r |= c->comp_hdr->Qs_codec->decode(s, c->comp_hdr->Qs_codec,
						   blk, qual, &out_sz2);
	    } else {
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
    } else if (ref_id >= 0 && s->ref != fd->ref_free) {
	cram_ref_decr(fd->refs, ref_id);
    }
    pthread_mutex_unlock(&fd->ref_lock);

    /* Resolve mate pair cross-references between recs within this slice */
    cram_decode_slice_xref(s);

    return r;
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

    if (!(j = malloc(sizeof(*j))))
	return -1;

    j->fd = fd;
    j->c  = c;
    j->s  = s;
    j->h  = bfd;
    
    nonblock = t_pool_results_queue_len(fd->rqueue) ? 0 : 1;

    if (-1 == t_pool_dispatch2(fd->pool, fd->rqueue, cram_decode_slice_thread,
			       j, nonblock)) {
	/* Would block */
	fd->job_pending = j;
    } else {
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

    /* Assign names if not explicitly set */
    if (cr->name_len) {
	name = (char *)BLOCK_DATA(s->name_blk) + cr->name;
	name_len = cr->name_len;
    } else {
	// FIXME: add prefix, container number, slice number, etc
	name = name_a;

	if (cr->mate_line >= 0 && cr->mate_line < rec)
	    name_len = sprintf(name_a, "%s:%"PRId64":%d",
			       fd->prefix, s->id, cr->mate_line);
	else
	    name_len = sprintf(name_a, "%s:%"PRId64":%d",
			       fd->prefix, s->id, rec);
    }

    /* Generate BAM record */
    if (cr->rg < -1 || cr->rg >= bfd->nrg)
	return -1;
    rg_len = (cr->rg != -1) ? bfd->rg[cr->rg].name_len + 4 : 0;

    if (!BLOCK_DATA(s->seqs_blk))
	return -1;
    if (!BLOCK_DATA(s->qual_blk))
	return -1;

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
                                (char *)BLOCK_DATA(s->seqs_blk) + cr->seq,
				(char *)BLOCK_DATA(s->qual_blk) + cr->qual);
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
    
#ifndef SAMTOOLS
    bam_set_blk_size(*bam, bam_blk_size(*bam) + (aux - aux_orig));
#endif

    *aux++ = 0;

    return bam_idx + (aux - aux_orig);
}

/*
 * Here be dragons! The multi-threading code in this is crufty beyond belief.
 */
static cram_slice *cram_next_slice(cram_fd *fd, cram_container **cp) {
    cram_container *c;
    cram_slice *s = NULL;

    fd->eof = 0;

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
	    pthread_mutex_lock(&fd->ref_lock);
	    fd->unsorted = 1;
	    pthread_mutex_unlock(&fd->ref_lock);
	}
    }

    if ((s = c->slice))
	cram_free_slice(s);

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
		    if (c->ref_seq_id != fd->range.refid) {
			fd->eof = 1;
			return NULL;
		    }

		    if (c->ref_seq_start > fd->range.end) {
			fd->eof = 1;
			return NULL;
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

	if (t_pool_results_queue_sz(fd->rqueue) > fd->pool->qsize)
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

	t_pool_delete_result(res, 1);
    }

    *cp = c;
    return s;
}

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
	} else {
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

    return cram_to_bam(fd->header, fd, s, cr, c->curr_rec-1, bam);
}
