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

#ifdef SAMTOOLS
#    define bam_copy(dst, src) bam_copy1(*(dst), (src))
#else
void bam_copy(bam_seq_t **bt, bam_seq_t *bf) {
    size_t a;

    if (bf->alloc > (*bt)->alloc) {
	a = ((int)((bf->alloc+15)/16))*16;
	*bt = realloc(*bt, a);
	memcpy(*bt, bf, bf->alloc);
    } else {
	a = (*bt)->alloc;
	memcpy(*bt, bf, bf->alloc);
    }

    (*bt)->alloc = a;
}
#endif

#define Z_CRAM_STRAT Z_FILTERED
//#define Z_CRAM_STRAT Z_RLE
//#define Z_CRAM_STRAT Z_HUFFMAN_ONLY
//#define Z_CRAM_STRAT Z_DEFAULT_STRATEGY

static int process_one_read(cram_fd *fd, cram_container *c,
			    cram_slice *s, cram_record *cr,
			    bam_seq_t *b, int rnum);

/*
 * Returns index of val into key.
 * Basically strchr(key, val)-key;
 */
static int sub_idx(char *key, char val) {
    int i;

    for (i = 0; *key && *key++ != val; i++);
    return i;
}

/*
 * Encodes a compression header block into a generic cram_block structure.
 *
 * Returns cram_block ptr on success
 *         NULL on failure
 */
cram_block *cram_encode_compression_header(cram_fd *fd, cram_container *c,
					   cram_block_compression_hdr *h) {
    cram_block *cb  = cram_new_block(COMPRESSION_HEADER, 0);
    cram_block *map = cram_new_block(COMPRESSION_HEADER, 0);
    int i, mc;

    if (!cb || !map)
	return NULL;

    /*
     * This is a concatenation of several blocks of data:
     * header + landmarks, preservation map, read encoding map, and the tag
     * encoding map.
     * All 4 are variable sized and we need to know how large these are
     * before creating the compression header itself as this starts with
     * the total size (stored as a variable length string).
     */

    // Duplicated from container itself, and removed in 1.1
    if (fd->version == CRAM_1_VERS) {
	itf8_put_blk(cb, h->ref_seq_id);
	itf8_put_blk(cb, h->ref_seq_start);
	itf8_put_blk(cb, h->ref_seq_span);
	itf8_put_blk(cb, h->num_records);
	itf8_put_blk(cb, h->num_landmarks);
	for (i = 0; i < h->num_landmarks; i++) {
	    itf8_put_blk(cb, h->landmark[i]);
	}
    }

    /* Create in-memory preservation map */
    /* FIXME: should create this when we create the container */
    {
	khint_t k;
	int r;

	if (!(h->preservation_map = kh_init(map)))
	    return NULL;

	k = kh_put(map, h->preservation_map, "RN", &r);
	if (-1 == r) return NULL;
	kh_val(h->preservation_map, k).i = 1;

	if (fd->version == CRAM_1_VERS) {
	    k = kh_put(map, h->preservation_map, "PI", &r);
	    if (-1 == r) return NULL;
	    kh_val(h->preservation_map, k).i = 0;

	    k = kh_put(map, h->preservation_map, "UI", &r);
	    if (-1 == r) return NULL;
	    kh_val(h->preservation_map, k).i = 1;

	    k = kh_put(map, h->preservation_map, "MI", &r);
	    if (-1 == r) return NULL;
	    kh_val(h->preservation_map, k).i = 1;

	} else {
	    // Technically SM was in 1.0, but wasn't in Java impl.
	    k = kh_put(map, h->preservation_map, "SM", &r);
	    if (-1 == r) return NULL;
	    kh_val(h->preservation_map, k).i = 0;

	    k = kh_put(map, h->preservation_map, "TD", &r);
	    if (-1 == r) return NULL;
	    kh_val(h->preservation_map, k).i = 0;

	    k = kh_put(map, h->preservation_map, "AP", &r);
	    if (-1 == r) return NULL;
	    kh_val(h->preservation_map, k).i = c->pos_sorted;

	    if (fd->no_ref || fd->embed_ref) {
		// Reference Required == No
		k = kh_put(map, h->preservation_map, "RR", &r);
		if (-1 == r) return NULL;
		kh_val(h->preservation_map, k).i = 0;
	    }
	}
    }

    /* Encode preservation map; could collapse this and above into one */
    mc = 0;
    BLOCK_SIZE(map) = 0;
    if (h->preservation_map) {
	khint_t k;

	for (k = kh_begin(h->preservation_map);
	     k != kh_end(h->preservation_map);
	     k++) {
	    const char *key;
	    khash_t(map) *pmap = h->preservation_map;


	    if (!kh_exist(pmap, k))
		continue;

	    key = kh_key(pmap, k);
	    BLOCK_APPEND(map, key, 2);

	    switch(CRAM_KEY(key[0], key[1])) {
	    case CRAM_KEY('M','I'):
		BLOCK_APPEND_CHAR(map, kh_val(pmap, k).i);
		break;

	    case CRAM_KEY('U','I'):
		BLOCK_APPEND_CHAR(map, kh_val(pmap, k).i);
		break;

	    case CRAM_KEY('P','I'):
		BLOCK_APPEND_CHAR(map, kh_val(pmap, k).i);
		break;

	    case CRAM_KEY('A','P'):
		BLOCK_APPEND_CHAR(map, kh_val(pmap, k).i);
		break;

	    case CRAM_KEY('R','N'):
		BLOCK_APPEND_CHAR(map, kh_val(pmap, k).i);
		break;

	    case CRAM_KEY('R','R'):
		BLOCK_APPEND_CHAR(map, kh_val(pmap, k).i);
		break;

	    case CRAM_KEY('S','M'): {
		char smat[5], *mp = smat;
		*mp++ =
		    (sub_idx("CGTN", h->substitution_matrix[0][0]) << 6) |
		    (sub_idx("CGTN", h->substitution_matrix[0][1]) << 4) |
		    (sub_idx("CGTN", h->substitution_matrix[0][2]) << 2) |
		    (sub_idx("CGTN", h->substitution_matrix[0][3]) << 0);
		*mp++ =
		    (sub_idx("AGTN", h->substitution_matrix[1][0]) << 6) |
		    (sub_idx("AGTN", h->substitution_matrix[1][1]) << 4) |
		    (sub_idx("AGTN", h->substitution_matrix[1][2]) << 2) |
		    (sub_idx("AGTN", h->substitution_matrix[1][3]) << 0);
		*mp++ =
		    (sub_idx("ACTN", h->substitution_matrix[2][0]) << 6) |
		    (sub_idx("ACTN", h->substitution_matrix[2][1]) << 4) |
		    (sub_idx("ACTN", h->substitution_matrix[2][2]) << 2) |
		    (sub_idx("ACTN", h->substitution_matrix[2][3]) << 0);
		*mp++ =
		    (sub_idx("ACGN", h->substitution_matrix[3][0]) << 6) |
		    (sub_idx("ACGN", h->substitution_matrix[3][1]) << 4) |
		    (sub_idx("ACGN", h->substitution_matrix[3][2]) << 2) |
		    (sub_idx("ACGN", h->substitution_matrix[3][3]) << 0);
		*mp++ =
		    (sub_idx("ACGT", h->substitution_matrix[4][0]) << 6) |
		    (sub_idx("ACGT", h->substitution_matrix[4][1]) << 4) |
		    (sub_idx("ACGT", h->substitution_matrix[4][2]) << 2) |
		    (sub_idx("ACGT", h->substitution_matrix[4][3]) << 0);
		BLOCK_APPEND(map, smat, 5);
		break;
	    }

	    case CRAM_KEY('T','D'): {
		itf8_put_blk(map, BLOCK_SIZE(h->TD_blk));
		BLOCK_APPEND(map,
			     BLOCK_DATA(h->TD_blk),
			     BLOCK_SIZE(h->TD_blk));
		break;
	    }

	    default:
		fprintf(stderr, "Unknown preservation key '%.2s'\n", key);
		break;
	    }

	    mc++;
        }
    }
    itf8_put_blk(cb, BLOCK_SIZE(map) + itf8_size(mc));
    itf8_put_blk(cb, mc);    
    BLOCK_APPEND(cb, BLOCK_DATA(map), BLOCK_SIZE(map));
    
    /* rec encoding map */
    mc = 0;
    BLOCK_SIZE(map) = 0;
    if (h->BF_codec) {
	if (-1 == h->BF_codec->store(h->BF_codec, map, "BF", fd->version))
	    return NULL;
	mc++;
    }
    if (h->CF_codec) {
	if (-1 == h->CF_codec->store(h->CF_codec, map, "CF", fd->version))
	    return NULL;
	mc++;
    }
    if (h->RL_codec) {
	if (-1 == h->RL_codec->store(h->RL_codec, map, "RL", fd->version))
	    return NULL;
	mc++;
    }
    if (h->AP_codec) {
	if (-1 == h->AP_codec->store(h->AP_codec, map, "AP", fd->version))
	    return NULL;
	mc++;
    }
    if (h->RG_codec) {
	if (-1 == h->RG_codec->store(h->RG_codec, map, "RG", fd->version))
	    return NULL;
	mc++;
    }
    if (h->MF_codec) {
	if (-1 == h->MF_codec->store(h->MF_codec, map, "MF", fd->version))
	    return NULL;
	mc++;
    }
    if (h->NS_codec) {
	if (-1 == h->NS_codec->store(h->NS_codec, map, "NS", fd->version))
	    return NULL;
	mc++;
    }
    if (h->NP_codec) {
	if (-1 == h->NP_codec->store(h->NP_codec, map, "NP", fd->version))
	    return NULL;
	mc++;
    }
    if (h->TS_codec) {
	if (-1 == h->TS_codec->store(h->TS_codec, map, "TS", fd->version))
	    return NULL;
	mc++;
    }
    if (h->NF_codec) {
	if (-1 == h->NF_codec->store(h->NF_codec, map, "NF", fd->version))
	    return NULL;
	mc++;
    }
    if (h->TC_codec) {
	if (-1 == h->TC_codec->store(h->TC_codec, map, "TC", fd->version))
	    return NULL;
	mc++;
    }
    if (h->TN_codec) {
	if (-1 == h->TN_codec->store(h->TN_codec, map, "TN", fd->version))
	    return NULL;
	mc++;
    }
    if (h->TL_codec) {
	if (-1 == h->TL_codec->store(h->TL_codec, map, "TL", fd->version))
	    return NULL;
	mc++;
    }
    if (h->FN_codec) {
	if (-1 == h->FN_codec->store(h->FN_codec, map, "FN", fd->version))
	    return NULL;
	mc++;
    }
    if (h->FC_codec) {
	if (-1 == h->FC_codec->store(h->FC_codec, map, "FC", fd->version))
	    return NULL;
	mc++;
    }
    if (h->FP_codec) {
	if (-1 == h->FP_codec->store(h->FP_codec, map, "FP", fd->version))
	    return NULL;
	mc++;
    }
    if (h->BS_codec) {
	if (-1 == h->BS_codec->store(h->BS_codec, map, "BS", fd->version))
	    return NULL;
	mc++;
    }
    if (h->IN_codec) {
	if (-1 == h->IN_codec->store(h->IN_codec, map, "IN", fd->version))
	    return NULL;
	mc++;
    }
    if (h->DL_codec) {
	if (-1 == h->DL_codec->store(h->DL_codec, map, "DL", fd->version))
	    return NULL;
	mc++;
    }
    if (h->BA_codec) {
	if (-1 == h->BA_codec->store(h->BA_codec, map, "BA", fd->version))
	    return NULL;
	mc++;
    }
    if (h->MQ_codec) {
	if (-1 == h->MQ_codec->store(h->MQ_codec, map, "MQ", fd->version))
	    return NULL;
	mc++;
    }
    if (h->RN_codec) {
	if (-1 == h->RN_codec->store(h->RN_codec, map, "RN", fd->version))
	    return NULL;
	mc++;
    }
    if (h->QS_codec) {
	if (-1 == h->QS_codec->store(h->QS_codec, map, "QS", fd->version))
	    return NULL;
	mc++;
    }
    if (h->Qs_codec) {
	if (-1 == h->Qs_codec->store(h->Qs_codec, map, "Qs", fd->version))
	    return NULL;
	mc++;
    }
    if (h->RI_codec) {
	if (-1 == h->RI_codec->store(h->RI_codec, map, "RI", fd->version))
	    return NULL;
	mc++;
    }
    if (fd->version != CRAM_1_VERS) {
	if (h->SC_codec) {
	    if (-1 == h->SC_codec->store(h->SC_codec, map, "SC", fd->version))
		return NULL;
	    mc++;
	}
	if (h->RS_codec) {
	    if (-1 == h->RS_codec->store(h->RS_codec, map, "RS", fd->version))
		return NULL;
	    mc++;
	}
	if (h->PD_codec) {
	    if (-1 == h->PD_codec->store(h->PD_codec, map, "PD", fd->version))
		return NULL;
	    mc++;
	}
	if (h->HC_codec) {
	    if (-1 == h->HC_codec->store(h->HC_codec, map, "HC", fd->version))
		return NULL;
	    mc++;
	}
    }
    if (h->TM_codec) {
	if (-1 == h->TM_codec->store(h->TM_codec, map, "TM", fd->version))
	    return NULL;
	mc++;
    }
    if (h->TV_codec) {
	if (-1 == h->TV_codec->store(h->TV_codec, map, "TV", fd->version))
	    return NULL;
	mc++;
    }
    itf8_put_blk(cb, BLOCK_SIZE(map) + itf8_size(mc));
    itf8_put_blk(cb, mc);    
    BLOCK_APPEND(cb, BLOCK_DATA(map), BLOCK_SIZE(map));

    /* tag encoding map */
#if 0
    mp = map; mc = 0;
    if (h->tag_encoding_map) {
        HashItem *hi;
        HashIter *iter = HashTableIterCreate();
	if (!iter)
	    return NULL;

        while ((hi = HashTableIterNext(h->tag_encoding_map, iter))) {
            cram_map *m = hi->data.p;
	    int sz;

	    mp += itf8_put(mp, (hi->key[0]<<16)|(hi->key[1]<<8)|hi->key[2]);
	    if (-1 == (sz = m->codec->store(m->codec, mp, NULL, fd->version)))
		return NULL;
	    mp += sz;
	    mc++;
        }

        HashTableIterDestroy(iter);
    }
#else
    mc = 0;
    BLOCK_SIZE(map) = 0;
    if (c->tags_used) {
	khint_t k;

	for (k = kh_begin(c->tags_used); k != kh_end(c->tags_used); k++) {
	    if (!kh_exist(c->tags_used, k))
		continue;

	    mc++;
	    itf8_put_blk(map, kh_key(c->tags_used, k));

	    // use block content id 4
	    switch(kh_key(c->tags_used, k) & 0xff) {
	    case 'Z': case 'H':
		// string as byte_array_stop
		if (fd->version == CRAM_1_VERS) {
		    BLOCK_APPEND(map,
				 "\005" // BYTE_ARRAY_STOP
				 "\005" // len
				 "\t"   // stop-byte is also SAM separator
				 CRAM_EXT_TAG_S "\000\000\000",
				 7);
		} else {
		    BLOCK_APPEND(map,
				 "\005" // BYTE_ARRAY_STOP
				 "\002" // len
				 "\t"   // stop-byte is also SAM separator
				 CRAM_EXT_TAG_S,
				 4);
		}
		break;

	    case 'A': case 'c': case 'C':
		// byte array len, 1 byte
		BLOCK_APPEND(map,
			     "\004" // BYTE_ARRAY_LEN
			     "\011" // length
			     "\003" // HUFFMAN (len)
			     "\004" // huffman-len
			     "\001" // 1 symbol
			     "\001" // symbol=1 byte value
			     "\001" // 1 length
			     "\000" // length=0
			     "\001" // EXTERNAL (val)
			     "\001" // external-len
			     CRAM_EXT_TAG_S,// content-id
			     11);
		break;

	    case 's': case 'S':
		// byte array len, 2 byte
		BLOCK_APPEND(map,
			     "\004" // BYTE_ARRAY_LEN
			     "\011" // length
			     "\003" // HUFFMAN (len)
			     "\004" // huffman-len
			     "\001" // 1 symbol
			     "\002" // symbol=2 byte value
			     "\001" // 1 length
			     "\000" // length=0
			     "\001" // EXTERNAL (val)
			     "\001" // external-len
			     CRAM_EXT_TAG_S,// content-id
			     11);
		break;

	    case 'i': case 'I': case 'f':
		// byte array len, 4 byte
		BLOCK_APPEND(map,
			     "\004" // BYTE_ARRAY_LEN
			     "\011" // length
			     "\003" // HUFFMAN (len)
			     "\004" // huffman-len
			     "\001" // 1 symbol
			     "\004" // symbol=4 byte value
			     "\001" // 1 length
			     "\000" // length=0
			     "\001" // EXTERNAL (val)
			     "\001" // external-len
			     CRAM_EXT_TAG_S,// content-id
			     11);
		break;

	    case 'B':
		// Byte array of variable size, but we generate our tag
		// byte stream at the wrong stage (during reading and not
		// after slice header construction). So we use
		// BYTE_ARRAY_LEN with the length codec being external
		// too.
		BLOCK_APPEND(map,
			     "\004" // BYTE_ARRAY_LEN
			     "\006" // length
			     "\001" // EXTERNAL (len)
			     "\001" // external-len
			     "\004" // content-id
			     "\001" // EXTERNAL (val)
			     "\001" // external-len
			     CRAM_EXT_TAG_S,// content-id
			     8);
		break;

	    default:
		fprintf(stderr, "Unsupported SAM aux type '%c'\n",
			kh_key(c->tags_used, k) & 0xff);
	    }
	    //mp += m->codec->store(m->codec, mp, NULL, fd->version);
	}
    }
#endif
    itf8_put_blk(cb, BLOCK_SIZE(map) + itf8_size(mc));
    itf8_put_blk(cb, mc);    
    BLOCK_APPEND(cb, BLOCK_DATA(map), BLOCK_SIZE(map));

    if (fd->verbose)
	fprintf(stderr, "Wrote compression block header in %d bytes\n",
		(int)BLOCK_SIZE(cb));

    BLOCK_UPLEN(cb);

    cram_free_block(map);

    return cb;
}


/*
 * Encodes a slice compression header. 
 *
 * Returns cram_block on success
 *         NULL on failure
 */
cram_block *cram_encode_slice_header(cram_fd *fd, cram_slice *s) {
    char *buf;
    char *cp;
    cram_block *b = cram_new_block(MAPPED_SLICE, 0);
    int j;

    if (!b)
	return NULL;

    if (NULL == (cp = buf = malloc(16+5*(8+s->hdr->num_blocks)))) {
	cram_free_block(b);
	return NULL;
    }

    cp += itf8_put(cp, s->hdr->ref_seq_id);
    cp += itf8_put(cp, s->hdr->ref_seq_start);
    cp += itf8_put(cp, s->hdr->ref_seq_span);
    cp += itf8_put(cp, s->hdr->num_records);
    if (fd->version != CRAM_1_VERS)
	cp += itf8_put(cp, s->hdr->record_counter);
    cp += itf8_put(cp, s->hdr->num_blocks);
    cp += itf8_put(cp, s->hdr->num_content_ids);
    for (j = 0; j < s->hdr->num_content_ids; j++) {
	cp += itf8_put(cp, s->hdr->block_content_ids[j]);
    }
    if (s->hdr->content_type == MAPPED_SLICE)
	cp += itf8_put(cp, s->hdr->ref_base_id);

    if (fd->version != CRAM_1_VERS) {
	memcpy(cp, s->hdr->md5, 16); cp += 16;
    }
    
    assert(cp-buf <= 16+5*(8+s->hdr->num_blocks));

    b->data = (unsigned char *)buf;
    b->comp_size = b->uncomp_size = cp-buf;

    return b;
}


/*
 * Encodes a single slice from a container
 * FIXME: break into smaller components.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int cram_encode_slice(cram_fd *fd, cram_container *c,
			     cram_block_compression_hdr *h, cram_slice *s) {
    int rec, r = 0, last_pos;
    cram_block *core;
    int nblk, embed_ref;

    embed_ref = fd->embed_ref && s->hdr->ref_seq_id != -1 ? 1 : 0;

    /*
     * Slice external blocks:
     * ID 0 => base calls (insertions, soft-clip)
     * ID 1 => qualities
     * ID 2 => names
     * ID 3 => TS (insert size), NP (next frag)
     * ID 4 => tag values
     * ID 5 => BA, ifdef BA_external
     * ID 6 => tag IDs (TN), ifdef TN_external and CRAM_1_VERS
     * ID 7 => TD tag dictionary, if !CRAM_1_VERS
     */

    /* Create cram slice header, num_blocks etc */
    s->hdr->ref_base_id = embed_ref ? CRAM_EXT_REF : -1;
    s->hdr->record_counter = c->num_records + c->record_counter;
    c->num_records += s->hdr->num_records;
    nblk = (fd->version == CRAM_1_VERS) ? 5 : 6;
#ifdef BA_external
    nblk++;
#endif
#ifdef TN_external
    if (fd->version == CRAM_1_VERS) {
	nblk++;
    }
#endif
    if (embed_ref)
	nblk++;

    s->hdr->num_content_ids = nblk;
    s->hdr->num_blocks = s->hdr->num_content_ids+1;
    s->block = calloc(s->hdr->num_blocks, sizeof(s->block[0]));
    s->hdr->block_content_ids = malloc(s->hdr->num_content_ids *
				       sizeof(int32_t));
    if (!s->block || !s->hdr->block_content_ids)
	return -1;
    s->hdr->block_content_ids[0] = 0; // core
    s->hdr->block_content_ids[1] = CRAM_EXT_QUAL;
    s->hdr->block_content_ids[2] = CRAM_EXT_NAME;
    s->hdr->block_content_ids[3] = CRAM_EXT_TS_NP;
    s->hdr->block_content_ids[4] = CRAM_EXT_TAG;
    s->hdr->block_content_ids[5] = CRAM_EXT_SC;
    nblk = (fd->version == CRAM_1_VERS) ? 5 : 6;
#ifdef BA_external
    s->hdr->block_content_ids[(s->ba_id = ++nblk)-1] = CRAM_EXT_BA;
#endif
#ifdef TN_external
    if (fd->version == CRAM_1_VERS) {
	s->hdr->block_content_ids[(s->tn_id = ++nblk)-1] = CRAM_EXT_TN;
    }
#endif
    if (embed_ref)
	s->hdr->block_content_ids[(s->ref_id = ++nblk)-1] = CRAM_EXT_REF;

    if (!(s->block[0] = cram_new_block(CORE, 0)))                  return -1;
    if (!(s->block[1] = cram_new_block(EXTERNAL, CRAM_EXT_IN)))    return -1;
    if (!(s->block[2] = cram_new_block(EXTERNAL, CRAM_EXT_QUAL)))  return -1;
    if (!(s->block[3] = cram_new_block(EXTERNAL, CRAM_EXT_NAME)))  return -1;
    if (!(s->block[4] = cram_new_block(EXTERNAL, CRAM_EXT_TS_NP))) return -1;
    if (!(s->block[5] = cram_new_block(EXTERNAL, CRAM_EXT_TAG)))   return -1;
    if (fd->version != CRAM_1_VERS) {
	if (!(s->block[6] = cram_new_block(EXTERNAL, CRAM_EXT_SC)))
	    return -1;
    }
#ifdef BA_external
    if (!(s->block[s->ba_id] = cram_new_block(EXTERNAL, CRAM_EXT_BA)))
	return -1;
#endif
#ifdef TN_external
    if (fd->version == CRAM_1_VERS) {
	if (!(s->block[s->tn_id] = cram_new_block(EXTERNAL, CRAM_EXT_TN)))
	    return -1;
    }
#endif
    if (embed_ref) {
	if (!(s->block[s->ref_id] = cram_new_block(EXTERNAL, CRAM_EXT_REF)))
	    return -1;
	BLOCK_APPEND(s->block[s->ref_id],
		     c->ref + c->first_base - c->ref_start,
		     c->last_base - c->first_base + 1);
    }

    core = s->block[0];
		 
    /* Create a formal method for stealing from dstrings! */
    s->block[4]->data = calloc(10, s->hdr->num_records); // NP TS
    if (!s->block[4]->data)
	return -1;
    s->block[4]->comp_size = s->block[4]->uncomp_size = 0;

#ifdef BA_external
    s->block[s->ba_id]->data = calloc(1, s->BA_len);
    if (!s->block[s->ba_id]->data)
	return -1;
    s->block[s->ba_id]->comp_size = s->block[s->ba_id]->uncomp_size = 0;
#endif

    /* Generate core block */
    if (!(s->hdr_block = cram_encode_slice_header(fd, s)))
	return -1;

    last_pos = s->hdr->ref_seq_start;
    for (rec = 0; rec < s->hdr->num_records; rec++) {
	cram_record *cr = &s->crecs[rec];
	int32_t i32;
	unsigned char uc;

	//fprintf(stderr, "Encode seq %d, %d/%d FN=%d, %s\n", rec, core->byte, core->bit, cr->nfeature, s->name_ds->str + cr->name);

	//printf("BF=0x%x\n", cr->flags);
	//	    bf = cram_flag_swap[cr->flags];
	i32 = fd->cram_flag_swap[cr->flags & 0xfff];
	r |= h->BF_codec->encode(s, h->BF_codec, core, (char *)&i32, 1);

	i32 = cr->cram_flags;
	r |= h->CF_codec->encode(s, h->CF_codec, core,
				 (char *)&i32, 1);

	if (fd->version != CRAM_1_VERS)
	    r |= h->RI_codec->encode(s, h->RI_codec, core,
				     (char *)&cr->ref_id, 1);

	r |= h->RL_codec->encode(s, h->RL_codec, core,
				 (char *)&cr->len, 1);

	if (c->pos_sorted) {
	    i32 = cr->apos - last_pos;
	    r |= h->AP_codec->encode(s, h->AP_codec, core, (char *)&i32, 1);
	    last_pos = cr->apos;
	} else {
	    i32 = cr->apos;
	    r |= h->AP_codec->encode(s, h->AP_codec, core, (char *)&i32, 1);
	}

	r |= h->RG_codec->encode(s, h->RG_codec, core,
				 (char *)&cr->rg, 1);

	if (c->comp_hdr->read_names_included) {
	    // RN codec: Already stored in block[3].
	}

	if (cr->cram_flags & CRAM_FLAG_DETACHED) {
	    i32 = cr->mate_flags;
	    r |= h->MF_codec->encode(s, h->MF_codec, core, (char *)&i32, 1);

	    if (!c->comp_hdr->read_names_included) {
		// RN codec: Already stored in block[3].
	    }

#ifndef NS_external
	    r |= h->NS_codec->encode(s, h->NS_codec, core,
				     (char *)&cr->mate_ref_id, 1);
#else
	    s->block[4]->uncomp_size +=
		itf8_put(&s->block[4]->data[s->block[4]->uncomp_size],
			 cr->mate_ref_id);
#endif

#ifndef TS_external
	    r |= h->NP_codec->encode(s, h->NP_codec, core,
				     (char *)&cr->mate_pos, 1);

	    r |= h->TS_codec->encode(s, h->TS_codec, core,
				     (char *)&cr->tlen, 1);
#else
	    s->block[4]->uncomp_size +=
		itf8_put((char *)&s->block[4]->data[s->block[4]->uncomp_size],
			 cr->mate_pos);
	    s->block[4]->uncomp_size +=
		itf8_put((char *)&s->block[4]->data[s->block[4]->uncomp_size],
			 cr->tlen);
#endif
	} else if (cr->cram_flags & CRAM_FLAG_MATE_DOWNSTREAM) {
	    r |= h->NF_codec->encode(s, h->NF_codec, core,
				     (char *)&cr->mate_line, 1);
	}

	/* Aux tags */
	if (fd->version == CRAM_1_VERS) {
	    uc = cr->ntags;
	    r |= h->TC_codec->encode(s, h->TC_codec, core, (char *)&uc, 1);
#ifndef TN_external
	    {
		int j;
		for (j = 0; j < cr->ntags; j++) {
		    uint32_t i32 = s->TN[cr->TN_idx + j]; // id
		    r |= h->TN_codec->encode(s, h->TN_codec, core,
					     (char *)&i32, 1);
		}
	    }
#endif
	} else {
	    r |= h->TL_codec->encode(s, h->TL_codec, core, (char *)&cr->TL, 1);
	}

	// qual
	// QS codec : Already stored in block[2].

	// features (diffs)
	if (!(cr->flags & BAM_FUNMAP)) {
	    int prev_pos = 0, j;

	    r |= h->FN_codec->encode(s, h->FN_codec, core,
				     (char *)&cr->nfeature, 1);
	    for (j = 0; j < cr->nfeature; j++) {
		cram_feature *f = &s->features[cr->feature + j];

		uc = f->X.code;
		r |= h->FC_codec->encode(s, h->FC_codec, core,
					 (char *)&uc, 1);
		i32 = f->X.pos - prev_pos;
		r |= h->FP_codec->encode(s, h->FP_codec, core,
					 (char *)&i32, 1);
		prev_pos = f->X.pos;

		switch(f->X.code) {
		    //char *seq;

		case 'X':
		    //fprintf(stderr, "    FC=%c FP=%d base=%d\n", f->X.code, i32, f->X.base);
		
		    uc = f->X.base;
		    r |= h->BS_codec->encode(s, h->BS_codec, core,
					     (char *)&uc, 1);
		    break;
		case 'S':
		    //seq = DSTRING_STR(s->seqs_ds) + f->S.seq_idx;
		    //r |= h->SC_codec->encode(s, h->SC_codec, core,
		    //			     seq, f->S.len);
		    break;
		case 'I':
		    //seq = DSTRING_STR(s->seqs_ds) + f->S.seq_idx;
		    //r |= h->IN_codec->encode(s, h->IN_codec, core,
		    //			     seq, f->S.len);
		    break;
		case 'i':
		    uc = f->i.base;
#ifdef BA_external
		    s->block[s->ba_id]->data[s->block[s->ba_id]->uncomp_size++] = uc;
#else
		    r |= h->BA_codec->encode(s, h->BA_codec, core,
					     (char *)&uc, 1);
#endif
		    //seq = DSTRING_STR(s->seqs_ds) + f->S.seq_idx;
		    //r |= h->IN_codec->encode(s, h->IN_codec, core,
		    //			     seq, 1);
		    break;
		case 'D':
		    i32 = f->D.len;
		    r |= h->DL_codec->encode(s, h->DL_codec, core,
					     (char *)&i32, 1);
		    break;

		case 'B':
//		    // Used when we try to store a non ACGTN base or an N
//		    // that aligns against a non ACGTN reference

		    uc  = f->B.base;
#ifdef BA_external
		    s->block[s->ba_id]->data[s->block[s->ba_id]->uncomp_size++] = uc;
#else
		    r |= h->BA_codec->encode(s, h->BA_codec, core,
					     (char *)&uc, 1);
#endif

//                  Already added
//		    uc  = f->B.qual;
//		    r |= h->QS_codec->encode(s, h->QS_codec, core,
//					     (char *)&uc, 1);
		    break;

		case 'Q':
//                  Already added
//		    uc  = f->B.qual;
//		    r |= h->QS_codec->encode(s, h->QS_codec, core,
//					     (char *)&uc, 1);
		    break;

		case 'N':
		    i32 = f->N.len;
		    r |= h->RS_codec->encode(s, h->RS_codec, core,
					     (char *)&i32, 1);
		    break;
		    
		case 'P':
		    i32 = f->P.len;
		    r |= h->PD_codec->encode(s, h->PD_codec, core,
					     (char *)&i32, 1);
		    break;
		    
		case 'H':
		    i32 = f->H.len;
		    r |= h->HC_codec->encode(s, h->HC_codec, core,
					     (char *)&i32, 1);
		    break;
		    

		default:
		    fprintf(stderr, "unhandled feature code %c\n",
			    f->X.code);
		    return -1;
		}
	    }

	    r |= h->MQ_codec->encode(s, h->MQ_codec, core,
				     (char *)&cr->mqual, 1);
	} else {
	    char *seq = (char *)BLOCK_DATA(s->seqs_blk) + cr->seq;
#ifdef BA_external
	    memcpy(&s->block[s->ba_id]->data[s->block[s->ba_id]->uncomp_size],
		   seq, cr->len);
	    s->block[s->ba_id]->uncomp_size += cr->len;
#else
	    r |= h->BA_codec->encode(s, h->BA_codec, core, seq, cr->len);
#endif
	}

	if (r)
	    return -1;
    }
    s->block[0]->uncomp_size = s->block[0]->byte + (s->block[0]->bit < 7);
    s->block[0]->comp_size = s->block[0]->uncomp_size;

    // FIXME: we should avoid creating these in the first place and just
    // point them to s->base_blk et al.
    cram_free_block(s->block[1]);
    cram_free_block(s->block[2]);
    cram_free_block(s->block[3]);
    cram_free_block(s->block[5]);
    if (fd->version != CRAM_1_VERS) {
	cram_free_block(s->block[6]);
	BLOCK_UPLEN(s->soft_blk);
	s->block[6] = s->soft_blk;
	s->soft_blk = NULL;
    }
    BLOCK_UPLEN(s->base_blk); s->block[1] = s->base_blk; s->base_blk = NULL;
    BLOCK_UPLEN(s->qual_blk); s->block[2] = s->qual_blk; s->qual_blk = NULL;
    BLOCK_UPLEN(s->name_blk); s->block[3] = s->name_blk; s->name_blk = NULL;
    BLOCK_UPLEN(s->aux_blk);  s->block[5] = s->aux_blk;  s->aux_blk  = NULL;

#ifdef TN_external
    if (fd->version == CRAM_1_VERS) {
	cram_free_block(s->block[s->tn_id]);
	BLOCK_UPLEN(s->tn_blk); s->block[s->tn_id] = s->tn_blk;
	s->tn_blk = NULL;
    }
#endif

    s->block[4]->comp_size = s->block[4]->uncomp_size;
    
#ifdef BA_external
    s->block[s->ba_id]->comp_size = s->block[s->ba_id]->uncomp_size;
#endif

    /* Compress the CORE Block too, with minimal zlib level */
    if (fd->level > 5)
	cram_compress_block(fd, s->block[0], NULL, 1, Z_CRAM_STRAT, -1, -1);

#define USE_METRICS

#ifdef USE_METRICS
#  define LEVEL2 1
#  define STRAT2 Z_RLE
#else
#  define LEVEL2 -1
#  define STRAT2 -1
#endif

    /* Compress the other blocks */
    if (cram_compress_block(fd, s->block[1], NULL, //IN (seq)
			    fd->level, Z_CRAM_STRAT,
			    -1, -1))
	return -1;

    if (fd->level == 0) {
	/* Do nothing */
    } else if (fd->level == 1) {
	if (cram_compress_block(fd, s->block[2], fd->m[1], //qual
				1, Z_RLE, -1, -1))
	    return -1;
	if (cram_compress_block(fd, s->block[5], fd->m[4], //Tags
				1, Z_RLE, -1, -1))
	    return -1;
    } else if (fd->level < 3) {
	if (cram_compress_block(fd, s->block[2], fd->m[1], //qual
				1, Z_RLE,
				1, Z_HUFFMAN_ONLY))
	    return -1;
	if (cram_compress_block(fd, s->block[5], fd->m[4], //Tags
				1, Z_RLE,
				1, Z_HUFFMAN_ONLY))
	    return -1;
    } else {
	if (cram_compress_block(fd, s->block[2], fd->m[1], //qual
				fd->level, Z_CRAM_STRAT, 
				LEVEL2, STRAT2))
	    return -1;
	if (cram_compress_block(fd, s->block[5], fd->m[4], //Tags
				fd->level, Z_CRAM_STRAT,
				LEVEL2, STRAT2))
	    return -1;
    }
    if (cram_compress_block(fd, s->block[3], NULL, //Name
			    fd->level, Z_CRAM_STRAT,
			    -1, -1))
	return -1;
    if (cram_compress_block(fd, s->block[4], NULL, //TS, NP
			    fd->level, Z_CRAM_STRAT,
			    -1, -1))
	return -1;
    if (fd->version != CRAM_1_VERS) {
	if (cram_compress_block(fd, s->block[6], NULL, //SC (seq)
				fd->level, Z_CRAM_STRAT,
				-1, -1))
	    return -1;
    }
#ifdef BA_external
    if (cram_compress_block(fd, s->block[s->ba_id], NULL,
			    fd->level, Z_CRAM_STRAT, -1, -1))
	return -1;
#endif
#ifdef TN_external
    if (fd->version == CRAM_1_VERS) {
	if (cram_compress_block(fd, s->block[s->tn_id], NULL,
				fd->level, Z_DEFAULT_STRATEGY, -1, -1))
	    return -1;
    }
#endif
    if (embed_ref) {
	BLOCK_UPLEN(s->block[s->ref_id]);
	if (cram_compress_block(fd, s->block[s->ref_id], NULL,
				fd->level, Z_DEFAULT_STRATEGY, -1, -1))
	    return -1;
    }

    return r ? -1 : 0;
}

/*
 * Encodes all slices in a container into blocks.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_encode_container(cram_fd *fd, cram_container *c) {
    int i, j, slice_offset;
    cram_block_compression_hdr *h = c->comp_hdr;
    cram_block *c_hdr;
    int multi_ref = 0;
    int r1, r2, sn, nref;
    spare_bams *spares;

    /* Cache references up-front if we have unsorted access patterns */
    pthread_mutex_lock(&fd->ref_lock);
    nref = fd->refs->nref;
    pthread_mutex_unlock(&fd->ref_lock);

    if (c->refs_used) {
	for (i = 0; i < nref; i++) {
	    if (c->refs_used[i]) {	
		cram_get_ref(fd, i, 1, 0);
	    }
	}
    }

    /* Fetch reference sequence */
    if (!fd->no_ref) {
	bam_seq_t *b = c->bams[0];
	char *ref;

	ref = cram_get_ref(fd, bam_ref(b), 1, 0);
	if (!ref && bam_ref(b) >= 0) {
	    fprintf(stderr, "Failed to load reference #%d\n", bam_ref(b));
	    return -1;
	}
	if ((c->ref_id = bam_ref(b)) >= 0) {
	    c->ref_seq_id = c->ref_id;
	    c->ref       = fd->refs->ref_id[c->ref_seq_id]->seq;
	    c->ref_start = 1;
	    c->ref_end   = fd->refs->ref_id[c->ref_seq_id]->length;
	} else {
	    c->ref_seq_id = c->ref_id; // FIXME remove one var!
	}
    } else {
	    c->ref_seq_id = c->ref_id; // FIXME remove one var!
    }

    /* Turn bams into cram_records and gather basic stats */
    for (r1 = sn = 0; r1 < c->curr_c_rec; sn++) {
	cram_slice *s = c->slices[sn];
	int first_base = INT_MAX, last_base = INT_MIN;

	assert(sn < c->curr_slice);

	/* FIXME: we could create our slice objects here too instead of
	 * in cram_put_bam_seq. It's more natural here and also this is
	 * bit is threaded so it's less work in the main thread.
	 */

	for (r2 = 0; r1 < c->curr_c_rec && r2 < c->max_rec; r1++, r2++) {
	    cram_record *cr = &s->crecs[r2];
	    bam_seq_t *b = c->bams[r1];

	    /* If multi-ref we need to cope with changing reference per seq */
	    if (c->multi_seq && !fd->no_ref) {
		if (bam_ref(b) != c->ref_seq_id && bam_ref(b) >= 0) {
		    if (c->ref_seq_id >= 0)
			cram_ref_decr(fd->refs, c->ref_seq_id);

		    if (!cram_get_ref(fd, bam_ref(b), 1, 0)) {
			fprintf(stderr, "Failed to load reference #%d\n",
				bam_ref(b));
			return -1;
		    }

		    c->ref_seq_id = bam_ref(b); // overwritten later by -2
		    assert(fd->refs->ref_id[c->ref_seq_id]->seq);
		    c->ref       = fd->refs->ref_id[c->ref_seq_id]->seq;
		    c->ref_start = 1;
		    c->ref_end   = fd->refs->ref_id[c->ref_seq_id]->length;
		}
	    }

	    process_one_read(fd, c, s, cr, b, r2);

	    if (first_base > cr->apos)
		first_base = cr->apos;

	    if (last_base < cr->aend)
		last_base = cr->aend;
	}

	if (c->multi_seq) {
	    s->hdr->ref_seq_id    = -2;
	    s->hdr->ref_seq_start = 0;
	    s->hdr->ref_seq_span  = 0;
	} else {
	    s->hdr->ref_seq_id    = c->ref_id;
	    s->hdr->ref_seq_start = first_base;
	    s->hdr->ref_seq_span  = last_base - first_base + 1;
	}
	s->hdr->num_records = r2;
    }

    /* Link our bams[] array onto the spare bam list for reuse */
    spares = malloc(sizeof(*spares));
    pthread_mutex_lock(&fd->bam_list_lock);
    spares->bams = c->bams;
    spares->next = fd->bl;
    fd->bl = spares;
    pthread_mutex_unlock(&fd->bam_list_lock);
    c->bams = NULL;

    /* Detect if a multi-seq container */
    cram_stats_encoding(fd, c->RI_stats);
    multi_ref = c->RI_stats->nvals > 1;

    if (multi_ref) {
	if (fd->verbose)
	    fprintf(stderr, "Multi-ref container\n");
	c->ref_seq_id = -2;
	c->ref_seq_start = 0;
	c->ref_seq_span = 0;
    }


    /* Compute MD5s */
    for (i = 0; i < c->curr_slice; i++) {
	cram_slice *s = c->slices[i];
	
	if (fd->version != CRAM_1_VERS) {
	    if (s->hdr->ref_seq_id >= 0 && c->multi_seq == 0 && !fd->no_ref) {
		MD5_CTX md5;
		MD5_Init(&md5);
		MD5_Update(&md5,
			   c->ref + s->hdr->ref_seq_start - c->ref_start,
			   s->hdr->ref_seq_span);
		MD5_Final(s->hdr->md5, &md5);
	    } else {
		memset(s->hdr->md5, 0, 16);
	    }
	}
    }

    c->num_records = 0;
    c->num_blocks = 0;
    c->length = 0;

    //fprintf(stderr, "=== BF ===\n");
    h->BF_codec = cram_encoder_init(cram_stats_encoding(fd, c->BF_stats),
				    c->BF_stats, E_INT, NULL,
				    fd->version);

    //fprintf(stderr, "=== CF ===\n");
    h->CF_codec = cram_encoder_init(cram_stats_encoding(fd, c->CF_stats),
				    c->CF_stats, E_INT, NULL,
				    fd->version);
//    fprintf(stderr, "=== RN ===\n");
//    h->RN_codec = cram_encoder_init(cram_stats_encoding(fd, c->RN_stats),
//				    c->RN_stats, E_BYTE_ARRAY, NULL,
//				    fd->version);

    //fprintf(stderr, "=== AP ===\n");
    if (c->pos_sorted) {
	h->AP_codec = cram_encoder_init(cram_stats_encoding(fd, c->AP_stats),
					c->AP_stats, E_INT, NULL,
					fd->version);
    } else {
	int p[2] = {0, c->max_apos};
	h->AP_codec = cram_encoder_init(E_BETA, NULL, E_INT, p, fd->version);
    }

    //fprintf(stderr, "=== RG ===\n");
    h->RG_codec = cram_encoder_init(cram_stats_encoding(fd, c->RG_stats),
				    c->RG_stats, E_INT, NULL,
				    fd->version);

    //fprintf(stderr, "=== MQ ===\n");
    h->MQ_codec = cram_encoder_init(cram_stats_encoding(fd, c->MQ_stats),
				    c->MQ_stats, E_INT, NULL,
				    fd->version);

    //fprintf(stderr, "=== NS ===\n");
#ifdef NS_external
    h->NS_codec = cram_encoder_init(E_EXTERNAL, NULL, E_INT,
				    (void *)CRAM_EXT_NS,
				    fd->version);
#else
    h->NS_codec = cram_encoder_init(cram_stats_encoding(fd, c->NS_stats),
				    c->NS_stats, E_INT, NULL,
				    fd->version);
#endif

    //fprintf(stderr, "=== MF ===\n");
    h->MF_codec = cram_encoder_init(cram_stats_encoding(fd, c->MF_stats),
				    c->MF_stats, E_INT, NULL,
				    fd->version);

#ifdef TS_external
    h->TS_codec = cram_encoder_init(E_EXTERNAL, NULL, E_INT,
				    (void *)CRAM_EXT_TS_NP,
				    fd->version);
    h->NP_codec = cram_encoder_init(E_EXTERNAL, NULL, E_INT,
				    (void *)CRAM_EXT_TS_NP,
				    fd->version);
#else
    //fprintf(stderr, "=== TS ===\n");
    h->TS_codec = cram_encoder_init(cram_stats_encoding(fd, c->TS_stats),
				    c->TS_stats, E_INT, NULL,
				    fd->version);
    //fprintf(stderr, "=== NP ===\n");
    h->NP_codec = cram_encoder_init(cram_stats_encoding(fd, c->NP_stats),
				    c->NP_stats, E_INT, NULL,
				    fd->version);
#endif

    //fprintf(stderr, "=== NF ===\n");
    h->NF_codec = cram_encoder_init(cram_stats_encoding(fd, c->NF_stats),
				    c->NF_stats, E_INT, NULL,
				    fd->version);

    //fprintf(stderr, "=== RL ===\n");
    h->RL_codec = cram_encoder_init(cram_stats_encoding(fd, c->RL_stats),
				    c->RL_stats, E_INT, NULL,
				    fd->version);

    //fprintf(stderr, "=== FN ===\n");
    h->FN_codec = cram_encoder_init(cram_stats_encoding(fd, c->FN_stats),
				    c->FN_stats, E_INT, NULL,
				    fd->version);

    //fprintf(stderr, "=== FC ===\n");
    h->FC_codec = cram_encoder_init(cram_stats_encoding(fd, c->FC_stats),
				    c->FC_stats, E_BYTE, NULL,
				    fd->version);

    //fprintf(stderr, "=== FP ===\n");
    h->FP_codec = cram_encoder_init(cram_stats_encoding(fd, c->FP_stats),
				    c->FP_stats, E_INT, NULL,
				    fd->version);

    //fprintf(stderr, "=== DL ===\n");
    h->DL_codec = cram_encoder_init(cram_stats_encoding(fd, c->DL_stats),
				    c->DL_stats, E_INT, NULL,
				    fd->version);

#ifdef BA_external
    h->BA_codec = cram_encoder_init(E_EXTERNAL, NULL, E_BYTE,
				    (void *)CRAM_EXT_BA,
				    fd->version);
#else
    //fprintf(stderr, "=== BA ===\n");
    h->BA_codec = cram_encoder_init(cram_stats_encoding(fd, c->BA_stats),
				    c->BA_stats, E_BYTE, NULL,
				    fd->version);
#endif

    //fprintf(stderr, "=== BS ===\n");
    h->BS_codec = cram_encoder_init(cram_stats_encoding(fd, c->BS_stats),
				    c->BS_stats, E_BYTE, NULL,
				    fd->version);

    if (fd->version == CRAM_1_VERS) {
	h->TL_codec = NULL;
	h->RI_codec = NULL;
	h->RS_codec = NULL;
	h->PD_codec = NULL;
	h->HC_codec = NULL;
	h->SC_codec = NULL;

	//fprintf(stderr, "=== TC ===\n");
	h->TC_codec = cram_encoder_init(cram_stats_encoding(fd, c->TC_stats),
					c->TC_stats, E_BYTE, NULL,
					fd->version);

    //fprintf(stderr, "=== TN ===\n");
#ifdef TN_external
	h->TN_codec = cram_encoder_init(E_EXTERNAL, NULL, E_INT,
					(void *)CRAM_EXT_TN,
					fd->version);
#else
	h->TN_codec = cram_encoder_init(cram_stats_encoding(fd, c->TN_stats),
					c->TN_stats, E_INT, NULL,
					fd->version);
#endif
    } else {
	int i2[2] = {0, CRAM_EXT_SC};

	h->TC_codec = NULL;
	h->TN_codec = NULL;

	//fprintf(stderr, "=== TL ===\n");
	h->TL_codec = cram_encoder_init(cram_stats_encoding(fd, c->TL_stats),
					c->TL_stats, E_INT, NULL,
					fd->version);


	//fprintf(stderr, "=== RI ===\n");
	h->RI_codec = cram_encoder_init(cram_stats_encoding(fd, c->RI_stats),
					c->RI_stats, E_INT, NULL,
					fd->version);

	//fprintf(stderr, "=== RS ===\n");
	h->RS_codec = cram_encoder_init(cram_stats_encoding(fd, c->RS_stats),
					c->RS_stats, E_INT, NULL,
					fd->version);

	//fprintf(stderr, "=== PD ===\n");
	h->PD_codec = cram_encoder_init(cram_stats_encoding(fd, c->PD_stats),
					c->PD_stats, E_INT, NULL,
					fd->version);

	//fprintf(stderr, "=== HC ===\n");
	h->HC_codec = cram_encoder_init(cram_stats_encoding(fd, c->HC_stats),
					c->HC_stats, E_INT, NULL,
					fd->version);

	//fprintf(stderr, "=== SC ===\n");
	h->SC_codec = cram_encoder_init(E_BYTE_ARRAY_STOP, NULL,
					E_BYTE_ARRAY, (void *)i2,
					fd->version);
    }
    
    //fprintf(stderr, "=== IN ===\n");
    {
	int i2[2] = {0, CRAM_EXT_IN};
	h->IN_codec = cram_encoder_init(E_BYTE_ARRAY_STOP, NULL,
					E_BYTE_ARRAY, (void *)i2,
					fd->version);
    }

    {
	//int i2[2] = {0, 1};
	//h->QS_codec = cram_encoder_init(E_BYTE_ARRAY_STOP, NULL, (void *)i2,
	//				    fd->version);
	h->QS_codec = cram_encoder_init(E_EXTERNAL, NULL, E_BYTE,
					(void *)CRAM_EXT_QUAL,
					fd->version);
    }
    {
	int i2[2] = {0, CRAM_EXT_NAME};
	h->RN_codec = cram_encoder_init(E_BYTE_ARRAY_STOP, NULL,
					E_BYTE_ARRAY, (void *)i2,
					fd->version);
    }


    /* Encode slices */
    for (i = 0; i < c->curr_slice; i++) {
	if (fd->verbose)
	    fprintf(stderr, "Encode slice %d\n", i);
	if (cram_encode_slice(fd, c, h, c->slices[i]) != 0)
	    return -1;
    }

    /* Create compression header */
    {
	h->ref_seq_id    = c->ref_seq_id;
	h->ref_seq_start = c->ref_seq_start;
	h->ref_seq_span  = c->ref_seq_span;
	h->num_records   = c->num_records;
	
	h->mapped_qs_included = 0;   // fixme
	h->unmapped_qs_included = 0; // fixme
	// h->...  fixme
	memcpy(h->substitution_matrix, CRAM_SUBST_MATRIX, 20);

	if (!(c_hdr = cram_encode_compression_header(fd, c, h)))
	    return -1;
    }

    /* Compute landmarks */
    /* Fill out slice landmarks */
    c->num_landmarks = c->curr_slice;
    c->landmark = malloc(c->num_landmarks * sizeof(*c->landmark));
    if (!c->landmark)
	return -1;

    /*
     * Slice offset starts after the first block, so we need to simulate
     * writing it to work out the correct offset
     */
    {
	slice_offset = c_hdr->method == RAW
	    ? c_hdr->uncomp_size
	    : c_hdr->comp_size;
	slice_offset += 2 +
	    itf8_size(c_hdr->content_id) +
	    itf8_size(c_hdr->comp_size) +
	    itf8_size(c_hdr->uncomp_size);
    }

    c->ref_seq_id    = c->slices[0]->hdr->ref_seq_id;
    c->ref_seq_start = c->slices[0]->hdr->ref_seq_start;
    c->ref_seq_span  = c->slices[0]->hdr->ref_seq_span;
    for (i = 0; i < c->curr_slice; i++) {
	cram_slice *s = c->slices[i];
	
	c->num_blocks += s->hdr->num_blocks + 2;
	c->landmark[i] = slice_offset;

	if (s->hdr->ref_seq_start + s->hdr->ref_seq_span >
	    c->ref_seq_start + c->ref_seq_span) {
	    c->ref_seq_span = s->hdr->ref_seq_start + s->hdr->ref_seq_span
		- c->ref_seq_start;
	}
	
	slice_offset += s->hdr_block->method == RAW
	    ? s->hdr_block->uncomp_size
	    : s->hdr_block->comp_size;

	slice_offset += 2 + 
	    itf8_size(s->hdr_block->content_id) +
	    itf8_size(s->hdr_block->comp_size) +
	    itf8_size(s->hdr_block->uncomp_size);

	for (j = 0; j < s->hdr->num_blocks; j++) {
	    slice_offset += 2 + 
		itf8_size(s->block[j]->content_id) +
		itf8_size(s->block[j]->comp_size) +
		itf8_size(s->block[j]->uncomp_size);

	    slice_offset += s->block[j]->method == RAW
		? s->block[j]->uncomp_size
		: s->block[j]->comp_size;
	}
    }
    c->length += slice_offset; // just past the final slice

    c->comp_hdr_block = c_hdr;

    if (c->ref_seq_id >= 0) {
	cram_ref_decr(fd->refs, c->ref_seq_id);
    }

    /* Cache references up-front if we have unsorted access patterns */
    if (c->refs_used) {
	for (i = 0; i < fd->refs->nref; i++) {
	    if (c->refs_used[i])
		cram_ref_decr(fd->refs, i);
	}
    }

    return 0;
}


/*
 * Adds a feature code to a read within a slice. For purposes of minimising
 * memory allocations and fragmentation we have one array of features for all
 * reads within the slice. We return the index into this array for this new
 * feature.
 *
 * Returns feature index on success
 *         -1 on failure.
 */
static int cram_add_feature(cram_container *c, cram_slice *s,
			    cram_record *r, cram_feature *f) {
    if (s->nfeatures >= s->afeatures) {
	s->afeatures = s->afeatures ? s->afeatures*2 : 1024;
	s->features = realloc(s->features, s->afeatures*sizeof(*s->features));
	if (!s->features)
	    return -1;
    }

    if (!r->nfeature++) {
	r->feature = s->nfeatures;
	cram_stats_add(c->FP_stats, f->X.pos);
    } else {
	cram_stats_add(c->FP_stats,
		       f->X.pos - s->features[r->feature + r->nfeature-2].X.pos);
    }
    cram_stats_add(c->FC_stats, f->X.code);

    s->features[s->nfeatures++] = *f;

    return 0;
}

static int cram_add_substitution(cram_fd *fd, cram_container *c,
				 cram_slice *s, cram_record *r,
				 int pos, char base, char qual, char ref) {
    cram_feature f;

    // seq=ACGTN vs ref=ACGT or seq=ACGT vs ref=ACGTN
    if (fd->L2[(uc)base]<4 || (fd->L2[(uc)base]<5 && fd->L2[(uc)ref]<4)) {
	f.X.pos = pos+1;
	f.X.code = 'X';
	f.X.base = fd->cram_sub_matrix[ref&0x1f][base&0x1f];
	cram_stats_add(c->BS_stats, f.X.base);
    } else {
	f.B.pos = pos+1;
	f.B.code = 'B';
	f.B.base = base;
	f.B.qual = qual;
	cram_stats_add(c->BA_stats, f.B.base);
	cram_stats_add(c->QS_stats, f.B.qual);
	BLOCK_APPEND_CHAR(s->qual_blk, qual);
    }
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_base(cram_fd *fd, cram_container *c,
			 cram_slice *s, cram_record *r,
			 int pos, char base, char qual) {
    cram_feature f;
    f.B.pos = pos+1;
    f.B.code = 'B';
    f.B.base = base;
    f.B.qual = qual;
#ifdef BA_external
    s->BA_len++;
#else
    cram_stats_add(c->BA_stats, base);
#endif
    cram_stats_add(c->QS_stats, qual);
    BLOCK_APPEND_CHAR(s->qual_blk, qual);
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_quality(cram_fd *fd, cram_container *c,
			    cram_slice *s, cram_record *r,
			    int pos, char qual) {
    cram_feature f;
    f.Q.pos = pos+1;
    f.Q.code = 'Q';
    f.Q.qual = qual;
    cram_stats_add(c->QS_stats, qual);
    BLOCK_APPEND_CHAR(s->qual_blk, qual);
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_deletion(cram_container *c, cram_slice *s, cram_record *r,
			     int pos, int len, char *base) {
    cram_feature f;
    f.D.pos = pos+1;
    f.D.code = 'D';
    f.D.len = len;
    cram_stats_add(c->DL_stats, len);
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_softclip(cram_container *c, cram_slice *s, cram_record *r,
			     int pos, int len, char *base, int version) {
    cram_feature f;
    f.S.pos = pos+1;
    f.S.code = 'S';
    f.S.len = len;
    if (version == CRAM_1_VERS) {
	f.S.seq_idx = BLOCK_SIZE(s->base_blk);
	BLOCK_APPEND(s->base_blk, base, len);
	BLOCK_APPEND_CHAR(s->base_blk, '\0');
    } else {
	f.S.seq_idx = BLOCK_SIZE(s->soft_blk);
	if (base) {
	    BLOCK_APPEND(s->soft_blk, base, len);
	} else {
	    int i;
	    for (i = 0; i < len; i++)
		BLOCK_APPEND_CHAR(s->soft_blk, 'N');
	}
	BLOCK_APPEND_CHAR(s->soft_blk, '\0');
    }
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_hardclip(cram_container *c, cram_slice *s, cram_record *r,
			     int pos, int len, char *base) {
    cram_feature f;
    f.S.pos = pos+1;
    f.S.code = 'H';
    f.S.len = len;
    cram_stats_add(c->HC_stats, len);
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_skip(cram_container *c, cram_slice *s, cram_record *r,
			     int pos, int len, char *base) {
    cram_feature f;
    f.S.pos = pos+1;
    f.S.code = 'N';
    f.S.len = len;
    cram_stats_add(c->RS_stats, len);
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_pad(cram_container *c, cram_slice *s, cram_record *r,
			     int pos, int len, char *base) {
    cram_feature f;
    f.S.pos = pos+1;
    f.S.code = 'P';
    f.S.len = len;
    cram_stats_add(c->PD_stats, len);
    return cram_add_feature(c, s, r, &f);
}

static int cram_add_insertion(cram_container *c, cram_slice *s, cram_record *r,
			      int pos, int len, char *base) {
    cram_feature f;
    f.I.pos = pos+1;
    if (len == 1) {
	char b = base ? *base : 'N';
	f.i.code = 'i';
	f.i.base = b;
#ifdef BA_external
	s->BA_len++;
#else
	cram_stats_add(c->BA_stats, b);
#endif
    } else {
	f.I.code = 'I';
	f.I.len = len;
	f.S.seq_idx = BLOCK_SIZE(s->base_blk);
	if (base) {
	    BLOCK_APPEND(s->base_blk, base, len);
	} else {
	    int i;
	    for (i = 0; i < len; i++)
		BLOCK_APPEND_CHAR(s->base_blk, 'N');
	}
	BLOCK_APPEND_CHAR(s->base_blk, '\0');
    }
    return cram_add_feature(c, s, r, &f);
}

/*
 * Encodes auxiliary data.
 * Returns the read-group parsed out of the BAM aux fields on success
 *         NULL on failure or no rg present (FIXME)
 */
static char *cram_encode_aux_1_0(cram_fd *fd, bam_seq_t *b, cram_container *c,
				 cram_slice *s, cram_record *cr) {
    char *aux, *tmp, *rg = NULL, *tmp_tn;
    int aux_size = bam_blk_size(b) -
	((char *)bam_aux(b) - (char *)&bam_ref(b));
	
    /* Worst case is 1 nul char on every ??:Z: string, so +33% */
    BLOCK_GROW(s->aux_blk, aux_size*1.34+1);
    tmp = (char *)BLOCK_END(s->aux_blk);

#ifdef TN_external
    BLOCK_GROW(s->tn_blk, aux_size);
    tmp_tn = (char *)BLOCK_END(s->tn_blk);
#endif

    aux = (char *)bam_aux(b);
#ifndef TN_external
    cr->TN_idx = s->nTN;
#endif
    while (aux[0] != 0) {
	int32_t i32;
	int r;

	if (aux[0] == 'R' && aux[1] == 'G' && aux[2] == 'Z') {
	    rg = &aux[3];
	    while (*aux++);
	    continue;
	}
	if (aux[0] == 'M' && aux[1] == 'D' && aux[2] == 'Z') {
	    while (*aux++);
	    continue;
	}
	if (aux[0] == 'N' && aux[1] == 'M') {
	    switch(aux[2]) {
	    case 'A': case 'C': case 'c': aux+=4; break;
	    case 'I': case 'i': case 'f': aux+=7; break;
	    default:
		fprintf(stderr, "Unhandled type code for NM tag\n");
		return NULL;
	    }
	    continue;
	}

	cr->ntags++;

	i32 = (aux[0]<<16) | (aux[1]<<8) | aux[2];
	kh_put(s_i2i, c->tags_used, i32, &r);
	if (-1 == r)
	    return NULL;

#ifndef TN_external
	if (s->nTN >= s->aTN) {
	    s->aTN = s->aTN ? s->aTN*2 : 1024;
	    if (!(s->TN = realloc(s->TN, s->aTN * sizeof(*s->TN))))
		return NULL;
	}
	s->TN[s->nTN++] = i32;
	cram_stats_add(c->TN_stats, i32);
#else
	tmp_tn += itf8_put(tmp_tn, i32);
#endif

	switch(aux[2]) {
	case 'A': case 'C': case 'c':
	    aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    *tmp++=*aux++;
	    break;

	case 'S': case 's':
	    aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    *tmp++=*aux++; *tmp++=*aux++;
	    break;

	case 'I': case 'i': case 'f':
	    aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    break;

	case 'd':
	    aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    break;

	case 'Z': case 'H':
	    aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    while ((*tmp++=*aux++));
	    *tmp++ = '\t'; // stop byte
	    break;

	case 'B': {
	    int type = aux[3], blen;
	    uint32_t count = (uint32_t)((((unsigned char *)aux)[4]<< 0) +
					(((unsigned char *)aux)[5]<< 8) +
					(((unsigned char *)aux)[6]<<16) +
					(((unsigned char *)aux)[7]<<24));
	    // skip TN field
	    aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;

	    // We use BYTE_ARRAY_LEN with external length, so store that first
	    switch (type) {
	    case 'c': case 'C':
		blen = count;
		break;
	    case 's': case 'S':
		blen = 2*count;
		break;
	    case 'i': case 'I': case 'f':
		blen = 4*count;
		break;
	    default:
		fprintf(stderr, "Unknown sub-type '%c' for aux type 'B'\n",
			type);
		return NULL;
		    
	    }

	    tmp += itf8_put(tmp, blen+5);

	    *tmp++=*aux++; // sub-type & length
	    *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;

	    // The tag data itself
	    memcpy(tmp, aux, blen); tmp += blen; aux += blen;

	    //cram_stats_add(c->aux_B_stats, blen);
	    break;
	}
	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", aux[2]);
	    return NULL;
	}
    }
    cram_stats_add(c->TC_stats, cr->ntags);

    cr->aux = BLOCK_SIZE(s->aux_blk);
    cr->aux_size = (uc *)tmp - (BLOCK_DATA(s->aux_blk) + cr->aux);
    BLOCK_SIZE(s->aux_blk) = (uc *)tmp - BLOCK_DATA(s->aux_blk);
    assert(s->aux_blk->byte <= s->aux_blk->alloc);

#ifdef TN_external
    cr->tn = BLOCK_SIZE(s->tn_blk);
    BLOCK_SIZE(s->tn_blk) = (uc *)tmp_tn - BLOCK_DATA(s->tn_blk);
    assert(s->tn_blk->byte <= s->tn_blk->alloc);
#endif

    return rg;
}

/*
 * Encodes auxiliary data. Largely duplicated from above, but done so to
 * keep it simple and avoid a myriad of version ifs.
 *
 * Returns the read-group parsed out of the BAM aux fields on success
 *         NULL on failure or no rg present (FIXME)
 */
static char *cram_encode_aux(cram_fd *fd, bam_seq_t *b, cram_container *c,
			     cram_slice *s, cram_record *cr) {
    char *aux, *orig, *tmp, *rg = NULL;
#ifdef SAMTOOLS
    int aux_size = bam_get_l_aux(b);
#else
    int aux_size = bam_blk_size(b) -
	((char *)bam_aux(b) - (char *)&bam_ref(b));
#endif
    cram_block *td_b = c->comp_hdr->TD_blk;
    int TD_blk_size = BLOCK_SIZE(td_b), new;
    char *key;
    khint_t k;


    /* Worst case is 1 nul char on every ??:Z: string, so +33% */
    BLOCK_GROW(s->aux_blk, aux_size*1.34+1);
    tmp = (char *)BLOCK_END(s->aux_blk);


    orig = aux = (char *)bam_aux(b);

    // Copy aux keys to td_b and aux values to s->aux_blk
    while (aux - orig < aux_size && aux[0] != 0) {
	uint32_t i32;
	int r;

	if (aux[0] == 'R' && aux[1] == 'G' && aux[2] == 'Z') {
	    rg = &aux[3];
	    while (*aux++);
	    continue;
	}
	if (aux[0] == 'M' && aux[1] == 'D' && aux[2] == 'Z') {
	    while (*aux++);
	    continue;
	}
	if (aux[0] == 'N' && aux[1] == 'M') {
	    switch(aux[2]) {
	    case 'A': case 'C': case 'c': aux+=4; break;
	    case 'S': case 's':           aux+=5; break;
	    case 'I': case 'i': case 'f': aux+=7; break;
	    default:
		fprintf(stderr, "Unhandled type code for NM tag\n");
		return NULL;
	    }
	    continue;
	}

	BLOCK_APPEND(td_b, aux, 3);

	i32 = (aux[0]<<16) | (aux[1]<<8) | aux[2];
	kh_put(s_i2i, c->tags_used, i32, &r);
	if (-1 == r)
	    return NULL;

	switch(aux[2]) {
	case 'A': case 'C': case 'c':
	    aux+=3;
	    *tmp++=*aux++;
	    break;

	case 'S': case 's':
	    aux+=3;
	    *tmp++=*aux++; *tmp++=*aux++;
	    break;

	case 'I': case 'i': case 'f':
	    aux+=3;
	    *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    break;

	case 'd':
	    aux+=3; //*tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;
	    *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;

	case 'Z': case 'H':
	    aux+=3;
	    while ((*tmp++=*aux++));
	    *tmp++ = '\t'; // stop byte
	    break;

	case 'B': {
	    int type = aux[3], blen;
	    uint32_t count = (uint32_t)((((unsigned char *)aux)[4]<< 0) +
					(((unsigned char *)aux)[5]<< 8) +
					(((unsigned char *)aux)[6]<<16) +
					(((unsigned char *)aux)[7]<<24));
	    // skip TN field
	    aux+=3;

	    // We use BYTE_ARRAY_LEN with external length, so store that first
	    switch (type) {
	    case 'c': case 'C':
		blen = count;
		break;
	    case 's': case 'S':
		blen = 2*count;
		break;
	    case 'i': case 'I': case 'f':
		blen = 4*count;
		break;
	    default:
		fprintf(stderr, "Unknown sub-type '%c' for aux type 'B'\n",
			type);
		return NULL;
		    
	    }

	    tmp += itf8_put(tmp, blen+5);

	    *tmp++=*aux++; // sub-type & length
	    *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++; *tmp++=*aux++;

	    // The tag data itself
	    memcpy(tmp, aux, blen); tmp += blen; aux += blen;

	    //cram_stats_add(c->aux_B_stats, blen);
	    break;
	}
	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", aux[2]);
	    return NULL;
	}
    }

    // FIXME: sort BLOCK_DATA(td_b) by char[3] triples
    
    // And and increment TD hash entry
    BLOCK_APPEND_CHAR(td_b, 0);

    // Duplicate key as BLOCK_DATA() can be realloced to a new pointer.
    key = string_ndup(c->comp_hdr->TD_keys, 
		      (char *)BLOCK_DATA(td_b) + TD_blk_size,
		      BLOCK_SIZE(td_b) - TD_blk_size);
    k = kh_put(m_s2i, c->comp_hdr->TD_hash, key, &new);
    if (new < 0) {
	return NULL;
    } else if (new == 0) {
	BLOCK_SIZE(td_b) = TD_blk_size;
    } else {
	kh_val(c->comp_hdr->TD_hash, k) = c->comp_hdr->nTL;
	c->comp_hdr->nTL++;
    }

    cr->TL = kh_val(c->comp_hdr->TD_hash, k);
    cram_stats_add(c->TL_stats, cr->TL);

    cr->aux = BLOCK_SIZE(s->aux_blk);
    cr->aux_size = (uc *)tmp - (BLOCK_DATA(s->aux_blk) + cr->aux);
    BLOCK_SIZE(s->aux_blk) = (uc *)tmp - BLOCK_DATA(s->aux_blk);
    assert(s->aux_blk->byte <= s->aux_blk->alloc);

    return rg;
}


/*
 * Handles creation of a new container or new slice, flushing any
 * existing containers when appropriate. 
 *
 * Really this is next slice, which may or may not lead to a new container.
 *
 * Returns cram_container pointer on success
 *         NULL on failure.
 */
static cram_container *cram_next_container(cram_fd *fd, bam_seq_t *b) {
    cram_container *c = fd->ctr;
    cram_slice *s;
    int i;

    /* First occurence */
    if (c->curr_ref == -2)
	c->curr_ref = bam_ref(b);

    if (c->slice) {
	s = c->slice;
	if (c->multi_seq) {
	    s->hdr->ref_seq_id    = -2;
	    s->hdr->ref_seq_start = 0;
	    s->hdr->ref_seq_span  = 0;
	} else {
	    s->hdr->ref_seq_id    = c->curr_ref;
	    s->hdr->ref_seq_start = c->first_base;
	    s->hdr->ref_seq_span  = c->last_base - c->first_base + 1;
	}
	s->hdr->num_records   = c->curr_rec;

	if (c->curr_slice == 0) {
	    if (c->ref_seq_id != s->hdr->ref_seq_id)
		c->ref_seq_id  = s->hdr->ref_seq_id;
	    c->ref_seq_start = c->first_base;
	}

	c->curr_slice++;
    }

    /* Flush container */
    if (c->curr_slice == c->max_slice ||
	(bam_ref(b) != c->curr_ref && !c->multi_seq)) {
	c->ref_seq_span = fd->last_base - c->ref_seq_start + 1;
	if (fd->verbose)
	    fprintf(stderr, "Flush container %d/%d..%d\n",
		    c->ref_seq_id, c->ref_seq_start,
		    c->ref_seq_start + c->ref_seq_span -1);

	/* Encode slices */
	if (fd->pool) {
	    if (-1 == cram_flush_container_mt(fd, c))
		return NULL;
	} else {
	    if (-1 == cram_flush_container(fd, c))
		return NULL;

	    // Move to sep func, as we need cram_flush_container for
	    // the closing phase to flush the partial container.
	    for (i = 0; i < c->max_slice; i++) {
		cram_free_slice(c->slices[i]);
		c->slices[i] = NULL;
	    }

	    c->slice = NULL;
	    c->curr_slice = 0;

	    /* Easy approach for purposes of freeing stats */
	    cram_free_container(c);
	}

	c = fd->ctr = cram_new_container(fd->seqs_per_slice,
					 fd->slices_per_container);
	if (!c)
	    return NULL;
	c->record_counter = fd->record_counter;
	c->curr_ref = bam_ref(b);
    }

    c->last_pos = c->first_base = c->last_base = bam_pos(b)+1;

    /* New slice */
    c->slice = c->slices[c->curr_slice] =
	cram_new_slice(MAPPED_SLICE, c->max_rec);
    if (!c->slice)
	return NULL;

    if (c->multi_seq) {
	c->slice->hdr->ref_seq_id = -2;
	c->slice->hdr->ref_seq_start = 0;
	c->slice->last_apos = 1;
    } else {
	c->slice->hdr->ref_seq_id = bam_ref(b);
	// wrong for unsorted data, will fix during encoding.
	c->slice->hdr->ref_seq_start = bam_pos(b)+1;
	c->slice->last_apos = bam_pos(b)+1;
    }

    c->curr_rec = 0;

    return c;
}

/*
 * Converts a single bam record into a cram record.
 * Possibly used within a thread.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
static int process_one_read(cram_fd *fd, cram_container *c,
			    cram_slice *s, cram_record *cr,
			    bam_seq_t *b, int rnum) {
    int i, fake_qual = 0;
    char *cp, *rg;
    char *ref, *seq, *qual;

    // FIXME: multi-ref containers

    ref = c->ref;

    //fprintf(stderr, "%s => %d\n", rg ? rg : "\"\"", cr->rg);

    // Fields to resolve later
    //cr->mate_line;    // index to another cram_record
    //cr->mate_flags;   // MF
    //cr->ntags;        // TC
    cr->ntags      = 0; //cram_stats_add(c->TC_stats, cr->ntags);
    if (fd->version == CRAM_1_VERS)
	rg = cram_encode_aux_1_0(fd, b, c, s, cr);
    else
	rg = cram_encode_aux(fd, b, c, s, cr);

    //cr->aux_size = b->blk_size - ((char *)bam_aux(b) - (char *)&bam_ref(b));
    //cr->aux = DSTRING_LEN(s->aux_ds);
    //dstring_nappend(s->aux_ds, bam_aux(b), cr->aux_size);

    /* Read group, identified earlier */
    if (rg) {
	SAM_RG *brg = sam_hdr_find_rg(fd->header, rg);
	cr->rg = brg ? brg->id : -1;
    } else if (fd->version == CRAM_1_VERS) {
	SAM_RG *brg = sam_hdr_find_rg(fd->header, "UNKNOWN");
	assert(brg);
    } else {
	cr->rg = -1;
    }
    cram_stats_add(c->RG_stats, cr->rg);

    
    cr->ref_id      = bam_ref(b);  cram_stats_add(c->RI_stats, cr->ref_id);
    cr->flags       = bam_flag(b);
    if (bam_cigar_len(b) == 0)
	cr->flags |= BAM_FUNMAP;
    cram_stats_add(c->BF_stats, fd->cram_flag_swap[cr->flags & 0xfff]);

    if (!fd->no_ref)
	cr->cram_flags = CRAM_FLAG_PRESERVE_QUAL_SCORES;
    else
	cr->cram_flags = 0;
    //cram_stats_add(c->CF_stats, cr->cram_flags);

    cr->len         = bam_seq_len(b); cram_stats_add(c->RL_stats, cr->len);
    c->num_bases   += cr->len;
    cr->apos        = bam_pos(b)+1;
    if (c->pos_sorted) {
	if (cr->apos < s->last_apos) {
	    c->pos_sorted = 0;
	} else {
	    cram_stats_add(c->AP_stats, cr->apos - s->last_apos);
	    s->last_apos = cr->apos;
	}
    } else {
	//cram_stats_add(c->AP_stats, cr->apos);
    }
    c->max_apos += (cr->apos > c->max_apos) * (cr->apos - c->max_apos);

    cr->name        = BLOCK_SIZE(s->name_blk);
    cr->name_len    = bam_name_len(b);
    cram_stats_add(c->RN_stats, cr->name_len);

    BLOCK_APPEND(s->name_blk, bam_name(b), bam_name_len(b));


    /*
     * This seqs_ds is largely pointless and it could reuse the same memory
     * over and over.
     * s->base_ds is what we need for encoding.
     */
    cr->seq         = BLOCK_SIZE(s->seqs_blk);
    cr->qual        = BLOCK_SIZE(s->qual_blk);
    BLOCK_GROW(s->seqs_blk, cr->len+1);
    BLOCK_GROW(s->qual_blk, cr->len);
    seq = cp = (char *)BLOCK_END(s->seqs_blk);

    *seq = 0;
    for (i = 0; i < cr->len; i++) {
	// FIXME: do 2 char at a time for efficiency
#ifdef SAMTOOLS
	cp[i] = seq_nt16_str[bam_seqi(bam_seq(b), i)];
#else
	cp[i] = bam_nt16_rev_table[bam_seqi(bam_seq(b), i)];
#endif
    }
    BLOCK_SIZE(s->seqs_blk) += cr->len;

    qual = cp = (char *)bam_qual(b);

    /* Copy and parse */
    if (!(cr->flags & BAM_FUNMAP)) {
	int32_t *cig_to, *cig_from;
	int apos = cr->apos-1, spos = 0;

	cr->cigar       = s->ncigar;
	cr->ncigar      = bam_cigar_len(b);
	while (cr->cigar + cr->ncigar >= s->cigar_alloc) {
	    s->cigar_alloc = s->cigar_alloc ? s->cigar_alloc*2 : 1024;
	    s->cigar = realloc(s->cigar, s->cigar_alloc * sizeof(*s->cigar));
	    if (!s->cigar)
		return -1;
	}

	cig_to = (int32_t *)s->cigar;
	cig_from = (int32_t *)bam_cigar(b);

	cr->feature = 0;
	cr->nfeature = 0;
	for (i = 0; i < cr->ncigar; i++) {
	    enum cigar_op cig_op = cig_from[i] & BAM_CIGAR_MASK;
	    int cig_len = cig_from[i] >> BAM_CIGAR_SHIFT;
	    cig_to[i] = cig_from[i];

	    /* Can also generate events from here for CRAM diffs */

	    switch (cig_op) {
		int l;

		// Don't trust = and X ops to be correct.
	    case BAM_CMATCH:
	    case BAM_CBASE_MATCH:
	    case BAM_CBASE_MISMATCH:
		//fprintf(stderr, "\nBAM_CMATCH\nR: %.*s\nS: %.*s\n",
		//	cig_len, &ref[apos], cig_len, &seq[spos]);
		l = 0;
		if (!fd->no_ref && cr->len) {
		    int end = cig_len+apos < c->ref_end
			? cig_len : c->ref_end - apos;
		    for (l = 0; l < end && seq[spos]; l++, apos++, spos++) {
			if (ref[apos] != seq[spos]) {
			    //fprintf(stderr, "Subst: %d; %c vs %c\n",
			    //	spos, ref[apos], seq[spos]);
			    if (cram_add_substitution(fd, c, s, cr, spos,
						      seq[spos], qual[spos],
						      ref[apos]))
				return -1;
			}
		    }
		}

		if (l < cig_len && cr->len) {
		    /* off end of sequence or non-ref based output */
		    for (; l < cig_len && seq[spos]; l++, spos++) {
			if (cram_add_base(fd, c, s, cr, spos,
					  seq[spos], qual[spos]))
			    return -1;
		    }
		    apos += cig_len;
		} else if (!cr->len) {
		    /* Seq "*" */
		    apos += cig_len;
		    spos += cig_len;
		}
		break;
		
	    case BAM_CDEL:
		if (cram_add_deletion(c, s, cr, spos, cig_len, &seq[spos]))
		    return -1;
		apos += cig_len;
		break;

	    case BAM_CREF_SKIP:
		if (cram_add_skip(c, s, cr, spos, cig_len, &seq[spos]))
		    return -1;
		apos += cig_len;
		break;

	    case BAM_CINS:
		if (cram_add_insertion(c, s, cr, spos, cig_len,
				       cr->len ? &seq[spos] : NULL))
		    return -1;
		if (fd->no_ref && cr->len) {
		    for (l = 0; l < cig_len; l++, spos++) {
			cram_add_quality(fd, c, s, cr, spos, qual[spos]);
		    }
		} else {
		    spos += cig_len;
		}
		break;

	    case BAM_CSOFT_CLIP:
		if (cram_add_softclip(c, s, cr, spos, cig_len,
				      cr->len ? &seq[spos] : NULL,
				      fd->version))
		    return -1;
		if (fd->no_ref) {
		    if (cr->len) {
			for (l = 0; l < cig_len; l++, spos++) {
			    cram_add_quality(fd, c, s, cr, spos, qual[spos]);
			}
		    } else {
			for (l = 0; l < cig_len; l++, spos++) {
			    cram_add_quality(fd, c, s, cr, spos, -1);
			}
		    }
		} else {
		    spos += cig_len;
		}
		break;

	    case BAM_CHARD_CLIP:
		if (cram_add_hardclip(c, s, cr, spos, cig_len, &seq[spos]))
		    return -1;
		break;
	
	    case BAM_CPAD:
		if (cram_add_pad(c, s, cr, spos, cig_len, &seq[spos]))
		    return -1;
		break;
	    }
	}
	fake_qual = spos;
	cr->aend = MIN(apos, c->ref_end);
	cram_stats_add(c->FN_stats, cr->nfeature);
    } else {
	// Unmapped
	cr->cram_flags |= CRAM_FLAG_PRESERVE_QUAL_SCORES;
	cr->cigar  = 0;
	cr->ncigar = 0;
	cr->nfeature = 0;
	cr->aend = cr->apos;
#ifdef BA_external
	s->BA_len += cr->len;
#else
	for (i = 0; i < cr->len; i++)
	    cram_stats_add(c->BA_stats, seq[i]);
#endif
    }

    /*
     * Append to the qual block now. We do this here as
     * cram_add_substitution() can generate BA/QS events which need to 
     * be in the qual block before we append the rest of the data.
     */
    if (cr->cram_flags & CRAM_FLAG_PRESERVE_QUAL_SCORES) {
	/* Special case of seq "*" */
	if (cr->len == 0) {
	    cram_stats_add(c->RL_stats, cr->len = fake_qual);
	    BLOCK_GROW(s->qual_blk, cr->len);
	    cp = (char *)BLOCK_END(s->qual_blk);
	    memset(cp, 255, cr->len);
	} else {
	    BLOCK_GROW(s->qual_blk, cr->len);
	    cp = (char *)BLOCK_END(s->qual_blk);
	    char *from = (char *)&bam_qual(b)[0];
	    char *to = &cp[0];
	    memcpy(to, from, cr->len);
	    //for (i = 0; i < cr->len; i++) cp[i] = from[i];
	}
	BLOCK_SIZE(s->qual_blk) += cr->len;
    } else {
	if (cr->len == 0) {
	    cram_stats_add(c->RL_stats, cr->len = cr->aend - cr->apos + 1);
	}
    }

    /* Now we know apos and aend both, update mate-pair information */
    {
	int new;
	khint_t k;

	//fprintf(stderr, "Checking %"PRId64"/%.*s\t", rnum,
	//	cr->name_len, DSTRING_STR(s->name_ds)+cr->name);
	if (cr->flags & BAM_FPAIRED) {
	    char *key = string_ndup(s->pair_keys,
				    (char *)BLOCK_DATA(s->name_blk)+cr->name,
				    cr->name_len);
	    if (!key)
		return -1;

	    k = kh_put(m_s2i, s->pair, key, &new);
	    if (-1 == new)
		return -1;
	    else if (new > 0)
		kh_val(s->pair, k) = rnum;
	} else {
	    new = 1;
	}

	if (new == 0) {
	    cram_record *p = &s->crecs[kh_val(s->pair, k)];
	    
	    //fprintf(stderr, "paired %"PRId64"\n", kh_val(s->pair, k));

	    // copy from p to cr
	    cr->mate_pos = p->apos;
	    cram_stats_add(c->NP_stats, cr->mate_pos);

	    cr->tlen = cr->aend - p->apos;
	    cram_stats_add(c->TS_stats, cr->tlen);

	    cr->mate_flags =
		((p->flags & BAM_FMUNMAP)   == BAM_FMUNMAP)   * CRAM_M_UNMAP +
		((p->flags & BAM_FMREVERSE) == BAM_FMREVERSE) * CRAM_M_REVERSE;
	    cram_stats_add(c->MF_stats, cr->mate_flags);

	    // copy from cr to p
	    cram_stats_del(c->NP_stats, p->mate_pos);
	    p->mate_pos = cr->apos;
	    cram_stats_add(c->NP_stats, p->mate_pos);

	    cram_stats_del(c->MF_stats, p->mate_flags);
	    p->mate_flags =
		((cr->flags & BAM_FMUNMAP)   == BAM_FMUNMAP)  * CRAM_M_UNMAP +
		((cr->flags & BAM_FMREVERSE) == BAM_FMREVERSE)* CRAM_M_REVERSE;
	    cram_stats_add(c->MF_stats, p->mate_flags);

	    cram_stats_del(c->TS_stats, p->tlen);
	    p->tlen = p->apos - cr->aend;
	    cram_stats_add(c->TS_stats, p->tlen);

	    // Clear detached from cr flags
	    //cram_stats_del(c->CF_stats, cr->cram_flags);
	    cr->cram_flags &= ~CRAM_FLAG_DETACHED;
	    cram_stats_add(c->CF_stats, cr->cram_flags);

	    // Clear detached from p flags and set downstream
	    cram_stats_del(c->CF_stats, p->cram_flags);
	    p->cram_flags  &= ~CRAM_FLAG_DETACHED;
	    p->cram_flags  |=  CRAM_FLAG_MATE_DOWNSTREAM;
	    cram_stats_add(c->CF_stats, p->cram_flags);

	    p->mate_line = rnum - (kh_val(s->pair, k) + 1);
	    cram_stats_add(c->NF_stats, p->mate_line);

	    kh_val(s->pair, k) = rnum;
	} else {
	    //fprintf(stderr, "unpaired\n");

	    /* Derive mate flags from this flag */
	    cr->mate_flags = 0;
	    if (bam_flag(b) & BAM_FMUNMAP)
		cr->mate_flags |= CRAM_M_UNMAP;
	    if (bam_flag(b) & BAM_FMREVERSE)
		cr->mate_flags |= CRAM_M_REVERSE;

	    cram_stats_add(c->MF_stats, cr->mate_flags);

	    cr->mate_pos    = MAX(bam_mate_pos(b)+1, 0);
	    cram_stats_add(c->NP_stats, cr->mate_pos);

	    cr->tlen        = bam_ins_size(b);
	    cram_stats_add(c->TS_stats, cr->tlen);

	    cr->cram_flags |= CRAM_FLAG_DETACHED;
	    cram_stats_add(c->CF_stats, cr->cram_flags);
	}
    }

    cr->mqual       = bam_map_qual(b);
    cram_stats_add(c->MQ_stats, cr->mqual);

    cr->mate_ref_id = bam_mate_ref(b);
    cram_stats_add(c->NS_stats, cr->mate_ref_id);

    if (!(bam_flag(b) & BAM_FUNMAP)) {
	if (c->first_base > cr->apos)
	    c->first_base = cr->apos;

	if (c->last_base < cr->aend)
	    c->last_base = cr->aend;
    }

    return 0;
}

/*
 * Write iterator: put BAM format sequences into a CRAM file.
 * We buffer up a containers worth of data at a time.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_put_bam_seq(cram_fd *fd, bam_seq_t *b) {
    cram_container *c;

    if (!fd->ctr) {
	fd->ctr = cram_new_container(fd->seqs_per_slice,
				     fd->slices_per_container);
	if (!fd->ctr)
	    return -1;
	fd->ctr->record_counter = fd->record_counter;
    }
    c = fd->ctr;

    if (!c->slice || c->curr_rec == c->max_rec ||
	(bam_ref(b) != c->curr_ref && c->curr_ref >= -1)) {
	int slice_rec, curr_rec, multi_seq = fd->multi_seq == 1;
	int curr_ref = c->slice ? c->curr_ref : bam_ref(b);


	/*
	 * Start packing slices when we routinely have under 1/4tr full.
	 *
	 * This option isn't available if we choose to embed references
	 * since we can only have one per slice.
	 */
	if (fd->multi_seq == -1 && c->curr_rec < c->max_rec/4+10 &&
	    fd->last_slice && fd->last_slice < c->max_rec/4+10 &&
	    !fd->embed_ref) {
	    if (fd->verbose && !c->multi_seq)
		fprintf(stderr, "Multi-ref enabled for this container\n");
	    multi_seq = 1;
	}

	slice_rec = c->slice_rec;
	curr_rec  = c->curr_rec;

	if (fd->version == CRAM_1_VERS ||
	    c->curr_rec == c->max_rec || fd->multi_seq != 1 || !c->slice)
	    if (NULL == (c = cram_next_container(fd, b)))
		return -1;

	/*
	 * Due to our processing order, some things we've already done we
	 * cannot easily undo. So when we first notice we should be packing
	 * multiple sequences per container we emit the small partial
	 * container as-is and then start a fresh one in a different mode.
	 */
	if (multi_seq) {
	    fd->multi_seq = 1;
	    c->multi_seq = 1;
	    c->pos_sorted = 0; // required atm for multi_seq slices

	    if (!c->refs_used) {
		pthread_mutex_lock(&fd->ref_lock);
		c->refs_used = calloc(fd->refs->nref, sizeof(int));
		pthread_mutex_unlock(&fd->ref_lock);
		if (!c->refs_used)
		    return -1;
	    }
	}

	fd->last_slice = curr_rec - slice_rec;
	c->slice_rec = c->curr_rec;

	// Have we seen this reference before?
	if (bam_ref(b) >= 0 && bam_ref(b) != curr_ref && !fd->embed_ref &&
	    !fd->unsorted) {
	    
	    if (!c->refs_used) {
		pthread_mutex_lock(&fd->ref_lock);
		c->refs_used = calloc(fd->refs->nref, sizeof(int));
		pthread_mutex_unlock(&fd->ref_lock);
		if (!c->refs_used)
		    return -1;
	    } else if (c->refs_used && c->refs_used[bam_ref(b)]) {
		fprintf(stderr, "Unsorted mode enabled\n");
		pthread_mutex_lock(&fd->ref_lock);
		fd->unsorted = 1;
		pthread_mutex_unlock(&fd->ref_lock);
		fd->multi_seq = 1;
	    }
	}

	c->curr_ref = bam_ref(b);
	if (c->refs_used && c->curr_ref >= 0) c->refs_used[c->curr_ref]++;
    }

    if (!c->bams) {
	/* First time through, allocate a set of bam pointers */
	pthread_mutex_lock(&fd->bam_list_lock);
	if (fd->bl) {
	    spare_bams *spare = fd->bl;
	    c->bams = spare->bams;
	    fd->bl = spare->next;
	    free(spare);
	} else {
	    c->bams = calloc(c->max_c_rec, sizeof(bam_seq_t *));
	    if (!c->bams)
		return -1;
	}
	pthread_mutex_unlock(&fd->bam_list_lock);
    }

    /* Copy or alloc+copy the bam record, for later encoding */
    if (c->bams[c->curr_c_rec])
	bam_copy(&c->bams[c->curr_c_rec], b);
    else
	c->bams[c->curr_c_rec] = bam_dup(b);

    c->curr_rec++;
    c->curr_c_rec++;
    fd->record_counter++;

    return 0;
}
