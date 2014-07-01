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
 * FIXME: add checking of cram_external_type to return NULL on unsupported
 * {codec,type} tuples.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include "cram/cram.h"

static char *codec2str(enum cram_encoding codec) {
    switch (codec) {
    case E_NULL:            return "NULL";
    case E_EXTERNAL:        return "EXTERNAL";
    case E_GOLOMB:          return "GOLOMB";
    case E_HUFFMAN:         return "HUFFMAN";
    case E_BYTE_ARRAY_LEN:  return "BYTE_ARRAY_LEN";
    case E_BYTE_ARRAY_STOP: return "BYTE_ARRAY_STOP";
    case E_BETA:            return "BETA";
    case E_SUBEXP:          return "SUBEXP";
    case E_GOLOMB_RICE:     return "GOLOMB_RICE";
    case E_GAMMA:           return "GAMMA";
    }

    return "(unknown)";
}

/*
 * ---------------------------------------------------------------------------
 * Block bit-level I/O functions.
 * All defined static here to promote easy inlining by the compiler.
 */

#if 0
/* Get a single bit, MSB first */
static signed int get_bit_MSB(cram_block *block) {
    unsigned int val;

    if (block->byte > block->alloc)
	return -1;

    val = block->data[block->byte] >> block->bit;
    if (--block->bit == -1) {
	block->bit = 7;
	block->byte++;
	//printf("(%02X)", block->data[block->byte]);
    }

    //printf("-B%d-", val&1);

    return val & 1;
}
#endif

/*
 * Count number of successive 0 and 1 bits
 */
static int get_one_bits_MSB(cram_block *block) {
    int n = 0, b;
    do {
	b = block->data[block->byte] >> block->bit;
	if (--block->bit == -1) {
	    block->bit = 7;
	    block->byte++;
	}
	n++;
    } while (b&1);

    return n-1;
}

static int get_zero_bits_MSB(cram_block *block) {
    int n = 0, b;
    do {
	b = block->data[block->byte] >> block->bit;
	if (--block->bit == -1) {
	    block->bit = 7;
	    block->byte++;
	}
	n++;
    } while (!(b&1));

    return n-1;
}

#if 0
/* Stores a single bit */
static void store_bit_MSB(cram_block *block, unsigned int bit) {
    if (block->byte >= block->alloc) {
	block->alloc = block->alloc ? block->alloc*2 : 1024;
	block->data = realloc(block->data, block->alloc);
    }

    if (bit)
	block->data[block->byte] |= (1 << block->bit);

    if (--block->bit == -1) {
	block->bit = 7;
	block->byte++;
	block->data[block->byte] = 0;
    }
}
#endif

#if 0
/* Rounds to the next whole byte boundary first */
static void store_bytes_MSB(cram_block *block, char *bytes, int len) {
    if (block->bit != 7) {
	block->bit = 7;
	block->byte++;
    }

    while (block->byte + len >= block->alloc) {
	block->alloc = block->alloc ? block->alloc*2 : 1024;
	block->data = realloc(block->data, block->alloc);
    }

    memcpy(&block->data[block->byte], bytes, len);
    block->byte += len;
}
#endif

/* Local optimised copy for inlining */
static inline unsigned int get_bits_MSB(cram_block *block, int nbits) {
    unsigned int val = 0;
    int i;

#if 0
    // Fits within the current byte */
    if (nbits <= block->bit+1) {
	val = (block->data[block->byte]>>(block->bit-(nbits-1))) & ((1<<nbits)-1);
	if ((block->bit -= nbits) == -1) {
	    block->bit = 7;
	    block->byte++;
	}
	return val;
    }

    // partial first byte
    val = block->data[block->byte] & ((1<<(block->bit+1))-1);
    nbits -= block->bit+1;
    block->bit = 7;
    block->byte++;

    // whole middle bytes
    while (nbits >= 8) {
	val = (val << 8) | block->data[block->byte++];
	nbits -= 8;
    }

    val <<= nbits;
    val |= (block->data[block->byte]>>(block->bit-(nbits-1))) & ((1<<nbits)-1);
    block->bit -= nbits;
    return val;
#endif

#if 0
    /* Inefficient implementation! */
    //printf("{");
    for (i = 0; i < nbits; i++)
	//val = (val << 1) | get_bit_MSB(block);
	GET_BIT_MSB(block, val);
#endif

#if 1
    /* Combination of 1st two methods */
    if (nbits <= block->bit+1) {
	val = (block->data[block->byte]>>(block->bit-(nbits-1))) & ((1<<nbits)-1);
	if ((block->bit -= nbits) == -1) {
	    block->bit = 7;
	    block->byte++;
	}
	return val;
    }

    switch(nbits) {
//    case 15: GET_BIT_MSB(block, val);
//    case 14: GET_BIT_MSB(block, val);
//    case 13: GET_BIT_MSB(block, val);
//    case 12: GET_BIT_MSB(block, val);
//    case 11: GET_BIT_MSB(block, val);
//    case 10: GET_BIT_MSB(block, val);
//    case  9: GET_BIT_MSB(block, val);
    case  8: GET_BIT_MSB(block, val);
    case  7: GET_BIT_MSB(block, val);
    case  6: GET_BIT_MSB(block, val);
    case  5: GET_BIT_MSB(block, val);
    case  4: GET_BIT_MSB(block, val);
    case  3: GET_BIT_MSB(block, val);
    case  2: GET_BIT_MSB(block, val);
    case  1: GET_BIT_MSB(block, val);
	break;

    default:
	for (i = 0; i < nbits; i++)
	    //val = (val << 1) | get_bit_MSB(block);
	    GET_BIT_MSB(block, val);
    }
#endif

    //printf("=0x%x}", val);

    return val;
}

/*
 * Can store up to 24-bits worth of data encoded in an integer value
 * Possibly we'd want to have a less optimal store_bits function when dealing
 * with nbits > 24, but for now we assume the codes generated are never
 * that big. (Given this is only possible with 121392 or more
 * characters with exactly the correct frequency distribution we check
 * for it elsewhere.)
 */
static int store_bits_MSB(cram_block *block, unsigned int val, int nbits) {
    /* fprintf(stderr, " store_bits: %02x %d\n", val, nbits); */

    /*
     * Use slow mode until we tweak the huffman generator to never generate
     * codes longer than 24-bits.
     */
    unsigned int mask;

    if (block->byte+4 >= block->alloc) {
	if (block->byte) {
	    block->alloc *= 2;
	    block->data = realloc(block->data, block->alloc + 4);
	    if (!block->data)
		return -1;
	} else {
	    block->alloc = 1024;
	    block->data = realloc(block->data, block->alloc + 4);
	    if (!block->data)
		return -1;
	    block->data[0] = 0; // initialise first byte of buffer
	}
    }

    
    
    if (nbits <= block->bit+1) {
	block->data[block->byte] |= (val << (block->bit+1-nbits));
	if ((block->bit-=nbits) == -1) {
	    block->bit = 7;
	    block->byte++;
	    block->data[block->byte] = 0;
	}
	return 0;
    }

    block->data[block->byte] |= (val >> (nbits -= block->bit+1));
    block->bit = 7;
    block->byte++;
    block->data[block->byte] = 0;
				 
    mask = 1<<(nbits-1);
    do {
	if (val & mask)
	    block->data[block->byte] |= (1 << block->bit);
	if (--block->bit == -1) {
	    block->bit = 7;
	    block->byte++;
	    block->data[block->byte] = 0;
	}
	mask >>= 1;
    } while(--nbits);

    return 0;
}

/*
 * Returns the next 'size' bytes from a block, or NULL if insufficient
 * data left.This is just a pointer into the block data and not an
 * allocated object, so do not free the result.
 */
static char *cram_extract_block(cram_block *b, int size) {
    char *cp = (char *)b->data + b->idx;
    b->idx += size;
    if (b->idx > b->uncomp_size)
	return NULL;

    return cp;
}

/*
 * ---------------------------------------------------------------------------
 * EXTERNAL
 */
int cram_external_decode_int(cram_slice *slice, cram_codec *c,
			     cram_block *in, char *out, int *out_size) {
    int i;
    char *cp;
    cram_block *b = NULL;

    /* Find the external block */
    if (slice->block_by_id) {
	if (!(b = slice->block_by_id[c->external.content_id]))
	    return -1;
    } else {
	for (i = 0; i < slice->hdr->num_blocks; i++) {
	    b = slice->block[i];
	    if (b->content_type == EXTERNAL &&
		b->content_id == c->external.content_id) {
		break;
	    }
	}
	if (i == slice->hdr->num_blocks || !b)
	    return -1;
    }

    cp = (char *)b->data + b->idx;
    // E_INT and E_LONG are guaranteed single item queries
    b->idx += itf8_get(cp, (int32_t *)out);
    *out_size = 1;

    return 0;
}

int cram_external_decode_char(cram_slice *slice, cram_codec *c,
			      cram_block *in, char *out,
			      int *out_size) {
    int i;
    char *cp;
    cram_block *b = NULL;

    /* Find the external block */
    if (slice->block_by_id) {
	if (!(b = slice->block_by_id[c->external.content_id]))
	    return -1;
    } else {
	for (i = 0; i < slice->hdr->num_blocks; i++) {
	    b = slice->block[i];
	    if (b->content_type == EXTERNAL &&
		b->content_id == c->external.content_id) {
		break;
	    }
	}
	if (i == slice->hdr->num_blocks || !b)
	    return -1;
    }

    cp = cram_extract_block(b, *out_size);
    if (!cp)
	return -1;

    memcpy(out, cp, *out_size);
    return 0;
}

int cram_external_decode_block(cram_slice *slice, cram_codec *c,
			      cram_block *in, char *out_,
			      int *out_size) {
    int i;
    char *cp;
    cram_block *b = NULL;
    cram_block *out = (cram_block *)out_;

    /* Find the external block */
    if (slice->block_by_id) {
	if (!(b = slice->block_by_id[c->external.content_id]))
	    return -1;
    } else {
	for (i = 0; i < slice->hdr->num_blocks; i++) {
	    b = slice->block[i];
	    if (b->content_type == EXTERNAL &&
		b->content_id == c->external.content_id) {
		break;
	    }
	}
	if (i == slice->hdr->num_blocks || !b)
	    return -1;
    }

    cp = cram_extract_block(b, *out_size);
    if (!cp)
	return -1;

    BLOCK_APPEND(out, cp, *out_size);
    return 0;
}

void cram_external_decode_free(cram_codec *c) {
    if (c)
	free(c);
}

cram_codec *cram_external_decode_init(char *data, int size,
				      enum cram_external_type option,
				      int version) {
    cram_codec *c;
    char *cp = data;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_EXTERNAL;
    if (option == E_INT || option == E_LONG)
	c->decode = cram_external_decode_int;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
	c->decode = cram_external_decode_char;
    else
	c->decode = cram_external_decode_block;
    c->free   = cram_external_decode_free;
    
    cp += itf8_get(cp, &c->external.content_id);

    if (cp - data != size) {
	fprintf(stderr, "Malformed external header stream\n");
	free(c);
	return NULL;
    }

    c->external.type = option;

    return c;
}

int cram_external_encode(cram_slice *slice, cram_codec *c,
			cram_block *out, char *in, int in_size) {
    uint32_t *i32 = (uint32_t *)in;

    itf8_put_blk(out, *i32);
    return 0;
}

void cram_external_encode_free(cram_codec *c) {
    if (!c)
	return;
    free(c);
}

int cram_external_encode_store(cram_codec *c, cram_block *b, char *prefix,
			       int version) {
    char tmp[99], *tp = tmp;
    int len = 0;

    if (prefix) {
	size_t l = strlen(prefix);
	BLOCK_APPEND(b, prefix, l);
	len += l;
    }

    tp += itf8_put(tp, c->e_external.content_id);
    len += itf8_put_blk(b, c->codec);
    len += itf8_put_blk(b, tp-tmp);
    BLOCK_APPEND(b, tmp, tp-tmp);
    len += tp-tmp;

    return len;
}

cram_codec *cram_external_encode_init(cram_stats *st,
				      enum cram_external_type option,
				      void *dat,
				      int version) {
    cram_codec *c;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec = E_EXTERNAL;
    c->free = cram_external_encode_free;
    c->encode = cram_external_encode;
    c->store = cram_external_encode_store;

    c->e_external.content_id = (size_t)dat;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * BETA
 */
int cram_beta_decode_int(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n;

    if (c->beta.nbits) {
	for (i = 0, n = *out_size; i < n; i++)
	    out_i[i] = get_bits_MSB(in, c->beta.nbits) - c->beta.offset;
    } else {
	for (i = 0, n = *out_size; i < n; i++)
	    out_i[i] = 0;
    }

    return 0;
}

int cram_beta_decode_char(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int i, n;

    if (c->beta.nbits) {
	for (i = 0, n = *out_size; i < n; i++)
	    out[i] = get_bits_MSB(in, c->beta.nbits) - c->beta.offset;
    } else {
	for (i = 0, n = *out_size; i < n; i++)
	    out[i] = 0;
    }

    return 0;
}

void cram_beta_decode_free(cram_codec *c) {
    if (c)
	free(c);
}

cram_codec *cram_beta_decode_init(char *data, int size,
				  enum cram_external_type option,
				  int version) {
    cram_codec *c;
    char *cp = data;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_BETA;
    if (option == E_INT || option == E_LONG)
	c->decode = cram_beta_decode_int;
    else if (option == E_BYTE_ARRAY || option == E_BYTE)
	c->decode = cram_beta_decode_char;
    else
	abort();
    c->free   = cram_beta_decode_free;
    
    cp += itf8_get(cp, &c->beta.offset);
    cp += itf8_get(cp, &c->beta.nbits);

    if (cp - data != size) {
	fprintf(stderr, "Malformed beta header stream\n");
	free(c);
	return NULL;
    }

    return c;
}

int cram_beta_encode_store(cram_codec *c, cram_block *b,
			   char *prefix, int version) {
    int len = 0;

    if (prefix) {
	size_t l = strlen(prefix);
	BLOCK_APPEND(b, prefix, l);
	len += l;
    }

    len += itf8_put_blk(b, c->codec);
    len += itf8_put_blk(b, itf8_size(c->e_beta.offset)
			+ itf8_size(c->e_beta.nbits)); // codec length
    len += itf8_put_blk(b, c->e_beta.offset);
    len += itf8_put_blk(b, c->e_beta.nbits);

    return len;
}

int cram_beta_encode_int(cram_slice *slice, cram_codec *c,
			 cram_block *out, char *in, int in_size) {
    int *syms = (int *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
	r |= store_bits_MSB(out, syms[i] + c->e_beta.offset, c->e_beta.nbits);

    return r;
}

int cram_beta_encode_char(cram_slice *slice, cram_codec *c,
			 cram_block *out, char *in, int in_size) {
    unsigned char *syms = (unsigned char *)in;
    int i, r = 0;

    for (i = 0; i < in_size; i++)
	r |= store_bits_MSB(out, syms[i] + c->e_beta.offset, c->e_beta.nbits);

    return r;
}

void cram_beta_encode_free(cram_codec *c) {
    if (c) free(c);
}

cram_codec *cram_beta_encode_init(cram_stats *st,
				  enum cram_external_type option,
				  void *dat,
				  int version) {
    cram_codec *c;
    int min_val, max_val, len = 0;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec  = E_BETA;
    c->free   = cram_beta_encode_free;
    if (option == E_INT)
	c->encode = cram_beta_encode_int;
    else
	c->encode = cram_beta_encode_char;
    c->store  = cram_beta_encode_store;

    if (dat) {
	min_val = ((int *)dat)[0];
	max_val = ((int *)dat)[1];
    } else {
	min_val = INT_MAX;
	max_val = INT_MIN;
	int i;
	for (i = 0; i < MAX_STAT_VAL; i++) {
	    if (!st->freqs[i])
		continue;
	    if (min_val > i)
		min_val = i;
	    max_val = i;
	}
	if (st->h) {
	    khint_t k;

	    for (k = kh_begin(st->h); k != kh_end(st->h); k++) {
		if (!kh_exist(st->h, k))
		    continue;
	    
		i = kh_key(st->h, k);
		if (min_val > i)
		    min_val = i;
		if (max_val < i)
		    max_val = i;
	    }
	}
    }

    assert(max_val >= min_val);
    c->e_beta.offset = -min_val;
    max_val -= min_val;
    while (max_val) {
	len++;
	max_val >>= 1;
    }
    c->e_beta.nbits = len;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * SUBEXP
 */
int cram_subexp_decode(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int n, count;
    int k = c->subexp.k;

    for (count = 0, n = *out_size; count < n; count++) {
	int i = 0, tail;
	int val;

	/* Get number of 1s */
	//while (get_bit_MSB(in) == 1) i++;
	i = get_one_bits_MSB(in);

	/*
	 * Val is
	 * i > 0:  2^(k+i-1) + k+i-1 bits
	 * i = 0:  k bits
	 */
	if (i) {
	    tail = i + k-1;
	    val = 0;
	    while (tail) {
		//val = val<<1; val |= get_bit_MSB(in);
		GET_BIT_MSB(in, val);
		tail--;
	    }
	    val += 1 << (i + k-1);
	} else {
	    tail = k;
	    val = 0;
	    while (tail) {
		//val = val<<1; val |= get_bit_MSB(in);
		GET_BIT_MSB(in, val);
		tail--;
	    }
	}

	out_i[count] = val - c->subexp.offset;
    }

    return 0;
}

void cram_subexp_decode_free(cram_codec *c) {
    if (c)
	free(c);
}

cram_codec *cram_subexp_decode_init(char *data, int size,
				    enum cram_external_type option,
				    int version) {
    cram_codec *c;
    char *cp = data;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_SUBEXP;
    c->decode = cram_subexp_decode;
    c->free   = cram_subexp_decode_free;
    
    cp += itf8_get(cp, &c->subexp.offset);
    cp += itf8_get(cp, &c->subexp.k);

    if (cp - data != size) {
	fprintf(stderr, "Malformed subexp header stream\n");
	free(c);
	return NULL;
    }

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * GAMMA
 */
int cram_gamma_decode(cram_slice *slice, cram_codec *c, cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n;

    for (i = 0, n = *out_size; i < n; i++) {
	int nz = 0;
	int val;
	//while (get_bit_MSB(in) == 0) nz++;
	nz = get_zero_bits_MSB(in);
	val = 1;
	while (nz > 0) {
	    //val <<= 1; val |= get_bit_MSB(in);
	    GET_BIT_MSB(in, val);
	    nz--;
	}

	out_i[i] = val - c->gamma.offset;
    }

    return 0;
}

void cram_gamma_decode_free(cram_codec *c) {
    if (c)
	free(c);
}

cram_codec *cram_gamma_decode_init(char *data, int size,
				   enum cram_external_type option,
				   int version) {
    cram_codec *c;
    char *cp = data;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_GAMMA;
    c->decode = cram_gamma_decode;
    c->free   = cram_gamma_decode_free;
    
    cp += itf8_get(cp, &c->gamma.offset);

    if (cp - data != size) {
	fprintf(stderr, "Malformed gamma header stream\n");
	free(c);
	return NULL;
    }

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * HUFFMAN
 */

static int code_sort(const void *vp1, const void *vp2) {
    const cram_huffman_code *c1 = (const cram_huffman_code *)vp1;
    const cram_huffman_code *c2 = (const cram_huffman_code *)vp2;

    if (c1->len != c2->len)
	return c1->len - c2->len;
    else
	return c1->symbol - c2->symbol;
}

void cram_huffman_decode_free(cram_codec *c) {
    if (!c)
	return;

    if (c->huffman.codes)
	free(c->huffman.codes);
    free(c);
}

int cram_huffman_decode_char0(cram_slice *slice, cram_codec *c,
			      cram_block *in, char *out, int *out_size) {
    int i, n;

    /* Special case of 0 length codes */
    for (i = 0, n = *out_size; i < n; i++) {
	out[i] = c->huffman.codes[0].symbol;
    }
    return 0;
}

int cram_huffman_decode_char(cram_slice *slice, cram_codec *c,
			     cram_block *in, char *out, int *out_size) {
    int i, n, ncodes = c->huffman.ncodes;
    const cram_huffman_code * const codes = c->huffman.codes;

    for (i = 0, n = *out_size; i < n; i++) {
	int idx = 0;
	int val = 0, len = 0, last_len = 0;

	for (;;) {
	    int dlen = codes[idx].len - last_len;
	    if (dlen <= 0 || (in->alloc - in->byte)*8 + in->bit + 7 < dlen)
		return -1;

	    //val <<= dlen;
	    //val  |= get_bits_MSB(in, dlen);
	    //last_len = (len  += dlen);

	    last_len = (len  += dlen);
	    for (; dlen; dlen--) GET_BIT_MSB(in, val);

	    idx = val - codes[idx].p;
	    if (idx >= ncodes || idx < 0)
		return -1;

	    if (codes[idx].code == val && codes[idx].len == len) {
		out[i] = codes[idx].symbol;
		break;
	    }
	}
    }

    return 0;
}

int cram_huffman_decode_int0(cram_slice *slice, cram_codec *c,
			     cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n;
    const cram_huffman_code * const codes = c->huffman.codes;

    /* Special case of 0 length codes */
    for (i = 0, n = *out_size; i < n; i++) {
	out_i[i] = codes[0].symbol;
    }
    return 0;
}

int cram_huffman_decode_int(cram_slice *slice, cram_codec *c,
			    cram_block *in, char *out, int *out_size) {
    int32_t *out_i = (int32_t *)out;
    int i, n, ncodes = c->huffman.ncodes;
    const cram_huffman_code * const codes = c->huffman.codes;

    for (i = 0, n = *out_size; i < n; i++) {
	int idx = 0;
	int val = 0, len = 0, last_len = 0;

	// Now one bit at a time for remaining checks
	for (;;) {
	    int dlen = codes[idx].len - last_len;
	    if (dlen <= 0 || (in->alloc - in->byte)*8 + in->bit + 7 < dlen)
		return -1;
	    
	    //val <<= dlen;
	    //val  |= get_bits_MSB(in, dlen);
	    //last_len = (len  += dlen);

	    last_len = (len  += dlen);
	    for (; dlen; dlen--) GET_BIT_MSB(in, val);

	    idx = val - codes[idx].p;
	    if (idx >= ncodes || idx < 0)
		return -1;

	    if (codes[idx].code == val && codes[idx].len == len) {
		out_i[i] = codes[idx].symbol;
		break;
	    }
	}
    }

    return 0;
}

/*
 * Initialises a huffman decoder from an encoding data stream.
 */
cram_codec *cram_huffman_decode_init(char *data, int size,
				     enum cram_external_type option,
				     int version) {
    int32_t ncodes, i, j;
    char *cp = data, *data_end = &data[size];
    cram_codec *h;
    cram_huffman_code *codes;
    int32_t val, last_len, max_len = 0;
    
    cp += itf8_get(cp, &ncodes);
    h = calloc(1, sizeof(*h));
    if (!h)
	return NULL;

    h->free   = cram_huffman_decode_free;

    h->huffman.ncodes = ncodes;
    codes = h->huffman.codes = malloc(ncodes * sizeof(*codes));
    if (!codes) {
	free(h);
	return NULL;
    }

    /* Read symbols and bit-lengths */
    for (i = 0; i < ncodes && cp < data_end; i++) {
	cp += itf8_get(cp, &codes[i].symbol);
    }

    if (cp >= data_end) {
	fprintf(stderr, "Malformed huffman header stream\n");
	free(h);
	return NULL;
    }
    cp += itf8_get(cp, &i);
    if (i != ncodes) {
	fprintf(stderr, "Malformed huffman header stream\n");
	free(h);
	return NULL;
    }

    if (ncodes == 0) {
	/* NULL huffman stream */
	return h;
    }

    for (i = 0; i < ncodes && cp < data_end; i++) {
	cp += itf8_get(cp, &codes[i].len);
	if (max_len < codes[i].len)
	    max_len = codes[i].len;
    }
    if (cp - data != size || max_len >= ncodes) {
	fprintf(stderr, "Malformed huffman header stream\n");
	free(h);
	return NULL;
    }

    /* Sort by bit length and then by symbol value */
    qsort(codes, ncodes, sizeof(*codes), code_sort);

    /* Assign canonical codes */
    val = -1, last_len = 0;
    for (i = 0; i < ncodes; i++) {
	val++;
	if (codes[i].len > last_len) {
	    while (codes[i].len > last_len) {
		val <<= 1;
		last_len++;
	    }
	}
	codes[i].code = val;
    }

    /*
     * Compute the next starting point, offset by the i'th value.
     * For example if codes 10, 11, 12, 13 are 30, 31, 32, 33 then
     * codes[10..13].p = 30 - 10.
     */
    last_len = 0;
    for (i = j = 0; i < ncodes; i++) {
	if (codes[i].len > last_len) {
	    j = codes[i].code - i;
	    last_len = codes[i].len;
	}
	codes[i].p = j;
    }

//    puts("==HUFF LEN==");
//    for (i = 0; i <= last_len+1; i++) {
//	printf("len %d=%d prefix %d\n", i, h->huffman.lengths[i], h->huffman.prefix[i]); 
//    }
//    puts("===HUFFMAN CODES===");
//    for (i = 0; i < ncodes; i++) {
//	int j;
//	printf("%d: %d %d %d ", i, codes[i].symbol, codes[i].len, codes[i].code);
//	j = codes[i].len;
//	while (j) {
//	    putchar(codes[i].code & (1 << --j) ? '1' : '0');
//	}
//	printf(" %d\n", codes[i].code);
//    }

    h->codec  = E_HUFFMAN;
    if (option == E_BYTE || option == E_BYTE_ARRAY) {
	if (h->huffman.codes[0].len == 0)
	    h->decode = cram_huffman_decode_char0;
	else
	    h->decode = cram_huffman_decode_char;
    } else if (option == E_BYTE_ARRAY_BLOCK) {
	abort();
    } else {
	if (h->huffman.codes[0].len == 0)
	    h->decode = cram_huffman_decode_int0;
	else
	    h->decode = cram_huffman_decode_int;
    }

    return (cram_codec *)h;
}

int cram_huffman_encode_char0(cram_slice *slice, cram_codec *c,
			      cram_block *out, char *in, int in_size) {
    return 0;
}

int cram_huffman_encode_char(cram_slice *slice, cram_codec *c,
			     cram_block *out, char *in, int in_size) {
    int i, code, len, r = 0;
    unsigned char *syms = (unsigned char *)in;

    do {
	int sym = *syms++;
	if (sym >= -1 && sym < MAX_HUFF) {
	    i = c->e_huffman.val2code[sym+1];
	    assert(c->e_huffman.codes[i].symbol == sym);
	    code = c->e_huffman.codes[i].code;
	    len  = c->e_huffman.codes[i].len;
	} else {
	    /* Slow - use a lookup table for when sym < MAX_HUFF? */
	    for (i = 0; i < c->e_huffman.nvals; i++) {
		if (c->e_huffman.codes[i].symbol == sym)
		    break;
	    }
	    if (i == c->e_huffman.nvals)
		return -1;
    
	    code = c->e_huffman.codes[i].code;
	    len  = c->e_huffman.codes[i].len;
	}

	r |= store_bits_MSB(out, code, len);
    } while (--in_size);

    return r;
}

int cram_huffman_encode_int0(cram_slice *slice, cram_codec *c,
			     cram_block *out, char *in, int in_size) {
    return 0;
}

int cram_huffman_encode_int(cram_slice *slice, cram_codec *c,
			    cram_block *out, char *in, int in_size) {
    int i, code, len, r = 0;
    int *syms = (int *)in;

    do {
	int sym = *syms++;

	if (sym >= -1 && sym < MAX_HUFF) {
	    i = c->e_huffman.val2code[sym+1];
	    assert(c->e_huffman.codes[i].symbol == sym);
	    code = c->e_huffman.codes[i].code;
	    len  = c->e_huffman.codes[i].len;
	} else {
	    /* Slow - use a lookup table for when sym < MAX_HUFFMAN_SYM? */
	    for (i = 0; i < c->e_huffman.nvals; i++) {
		if (c->e_huffman.codes[i].symbol == sym)
		    break;
	    }
	    if (i == c->e_huffman.nvals)
		return -1;
    
	    code = c->e_huffman.codes[i].code;
	    len  = c->e_huffman.codes[i].len;
	}

	r |= store_bits_MSB(out, code, len);
    } while (--in_size);

    return r;
}

void cram_huffman_encode_free(cram_codec *c) {
    if (!c)
	return;

    if (c->e_huffman.codes)
	free(c->e_huffman.codes);
    free(c);
}

/*
 * Encodes a huffman tree.
 * Returns number of bytes written.
 */
int cram_huffman_encode_store(cram_codec *c, cram_block *b, char *prefix,
			      int version) {
    int i, len = 0;
    cram_huffman_code *codes = c->e_huffman.codes;
    /*
     * Up to code length 127 means 2.5e+26 bytes of data required (worst
     * case huffman tree needs symbols with freqs matching the Fibonacci
     * series). So guaranteed 1 byte per code.
     *
     * Symbols themselves could be 5 bytes (eg -1 is 5 bytes in itf8).
     *
     * Therefore 6*ncodes + 5 + 5 + 1 + 5 is max memory
     */
    char *tmp = malloc(6*c->e_huffman.nvals+16);
    char *tp = tmp;

    if (!tmp)
	return -1;

    if (prefix) {
	size_t l = strlen(prefix);
	BLOCK_APPEND(b, prefix, l);
	len += l;
    }

    tp += itf8_put(tp, c->e_huffman.nvals);
    for (i = 0; i < c->e_huffman.nvals; i++) {
	tp += itf8_put(tp, codes[i].symbol);
    }

    tp += itf8_put(tp, c->e_huffman.nvals);
    for (i = 0; i < c->e_huffman.nvals; i++) {
	tp += itf8_put(tp, codes[i].len);
    }

    len += itf8_put_blk(b, c->codec);
    len += itf8_put_blk(b, tp-tmp);
    BLOCK_APPEND(b, tmp, tp-tmp);
    len += tp-tmp;

    free(tmp);

    return len;
}

cram_codec *cram_huffman_encode_init(cram_stats *st,
				     enum cram_external_type option,
				     void *dat,
				     int version) {
    int *vals = NULL, *freqs = NULL, vals_alloc = 0, *lens, code, len;
    int nvals, i, ntot = 0, max_val = 0, min_val = INT_MAX, k;
    cram_codec *c;
    cram_huffman_code *codes;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec = E_HUFFMAN;

    /* Count number of unique symbols */
    for (nvals = i = 0; i < MAX_STAT_VAL; i++) {
	if (!st->freqs[i])
	    continue;
	if (nvals >= vals_alloc) {
	    vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
	    vals  = realloc(vals,  vals_alloc * sizeof(int));
	    freqs = realloc(freqs, vals_alloc * sizeof(int));
	    if (!vals || !freqs) {
		if (vals)  free(vals);
		if (freqs) free(freqs);
		free(c);
		return NULL;
	    }
	}
	vals[nvals] = i;
	freqs[nvals] = st->freqs[i];
	assert(st->freqs[i] > 0);
	ntot += freqs[nvals];
	if (max_val < i) max_val = i;
	if (min_val > i) min_val = i;
	nvals++;
    }
    if (st->h) {
	khint_t k;

	for (k = kh_begin(st->h); k != kh_end(st->h); k++) {
	    if (!kh_exist(st->h, k))
		continue;
	    if (nvals >= vals_alloc) {
		vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
		vals  = realloc(vals,  vals_alloc * sizeof(int));
		freqs = realloc(freqs, vals_alloc * sizeof(int));
		if (!vals || !freqs)
		    return NULL;
	    }
	    vals[nvals]= kh_key(st->h, k);
	    freqs[nvals] = kh_val(st->h, k);
	    assert(freqs[nvals] > 0);
	    ntot += freqs[nvals];
	    if (max_val < i) max_val = i;
	    if (min_val > i) min_val = i;
	    nvals++;
	}
    }

    assert(nvals > 0);

    freqs = realloc(freqs, 2*nvals*sizeof(*freqs));
    lens = calloc(2*nvals, sizeof(*lens));
    if (!lens || !freqs)
	return NULL;

    /* Inefficient, use pointers to form chain so we can insert and maintain
     * a sorted list? This is currently O(nvals^2) complexity.
     */
    for (;;) {
	int low1 = INT_MAX, low2 = INT_MAX;
	int ind1 = 0, ind2 = 0;
	for (i = 0; i < nvals; i++) {
	    if (freqs[i] < 0)
		continue;
	    if (low1 > freqs[i]) 
		low2 = low1, ind2 = ind1, low1 = freqs[i], ind1 = i;
	    else if (low2 > freqs[i])
		low2 = freqs[i], ind2 = i;
	}
	if (low2 == INT_MAX)
	    break;

	freqs[nvals] = low1 + low2;
	lens[ind1] = nvals;
	lens[ind2] = nvals;
	freqs[ind1] *= -1;
	freqs[ind2] *= -1;
	nvals++;
    }
    nvals = nvals/2+1;

    /* Assign lengths */
    for (i = 0; i < nvals; i++) {
	int code_len = 0;
	for (k = lens[i]; k; k = lens[k])
	    code_len++;
	lens[i] = code_len;
	freqs[i] *= -1;
	//fprintf(stderr, "%d / %d => %d\n", vals[i], freqs[i], lens[i]);
    }


    /* Sort, need in a struct */
    if (!(codes = malloc(nvals * sizeof(*codes))))
	return NULL;
    for (i = 0; i < nvals; i++) {
	codes[i].symbol = vals[i];
	codes[i].len = lens[i];
    }
    qsort(codes, nvals, sizeof(*codes), code_sort);

    /*
     * Generate canonical codes from lengths.
     * Sort by length.
     * Start with 0.
     * Every new code of same length is +1.
     * Every new code of new length is +1 then <<1 per extra length.
     *
     * /\
     * a/\
     * /\/\
     * bcd/\
     *    ef
     * 
     * a 1  0
     * b 3  4 (0+1)<<2
     * c 3  5
     * d 3  6
     * e 4  14  (6+1)<<1
     * f 5  15     
     */
    code = 0; len = codes[0].len;
    for (i = 0; i < nvals; i++) {
	while (len != codes[i].len) {
	    code<<=1;
	    len++;
	}
	codes[i].code = code++;

	if (codes[i].symbol >= -1 && codes[i].symbol < MAX_HUFF)
	    c->e_huffman.val2code[codes[i].symbol+1] = i;

	//fprintf(stderr, "sym %d, code %d, len %d\n",
	//	codes[i].symbol, codes[i].code, codes[i].len);
    }

    free(lens);
    free(vals);
    free(freqs);

    c->e_huffman.codes = codes;
    c->e_huffman.nvals = nvals;

    c->free = cram_huffman_encode_free;
    if (option == E_BYTE || option == E_BYTE_ARRAY) {
	if (c->e_huffman.codes[0].len == 0)
	    c->encode = cram_huffman_encode_char0;
	else
	    c->encode = cram_huffman_encode_char;
    } else {
	if (c->e_huffman.codes[0].len == 0)
	    c->encode = cram_huffman_encode_int0;
	else
	    c->encode = cram_huffman_encode_int;
    }
    c->store = cram_huffman_encode_store;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * BYTE_ARRAY_LEN
 */
int cram_byte_array_len_decode(cram_slice *slice, cram_codec *c,
			       cram_block *in, char *out,
			       int *out_size) {
    /* Fetch length */
    int32_t len, one = 1;

    c->byte_array_len.len_codec->decode(slice, c->byte_array_len.len_codec, in, (char *)&len, &one);
    //printf("ByteArray Len=%d\n", len);

    if (c->byte_array_len.value_codec) {
	c->byte_array_len.value_codec->decode(slice,
					      c->byte_array_len.value_codec,
					      in, out, &len);
    } else {
	return -1;
    }

    *out_size = len;

    return 0;
}

void cram_byte_array_len_decode_free(cram_codec *c) {
    if (!c) return;

    if (c->byte_array_len.len_codec)
	c->byte_array_len.len_codec->free(c->byte_array_len.len_codec);

    if (c->byte_array_len.value_codec)
	c->byte_array_len.value_codec->free(c->byte_array_len.value_codec);

    free(c);
}

cram_codec *cram_byte_array_len_decode_init(char *data, int size,
					    enum cram_external_type option,
					    int version) {
    cram_codec *c;
    char *cp = data;
    int32_t encoding;
    int32_t sub_size;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_BYTE_ARRAY_LEN;
    c->decode = cram_byte_array_len_decode;
    c->free   = cram_byte_array_len_decode_free;
    
    cp += itf8_get(cp, &encoding);
    cp += itf8_get(cp, &sub_size);
    c->byte_array_len.len_codec = cram_decoder_init(encoding, cp, sub_size,
						    E_INT, version);
    cp += sub_size;

    cp += itf8_get(cp, &encoding);
    cp += itf8_get(cp, &sub_size);
    c->byte_array_len.value_codec = cram_decoder_init(encoding, cp, sub_size,
						      option, version);
    cp += sub_size;

    if (cp - data != size) {
	fprintf(stderr, "Malformed byte_array_len header stream\n");
	free(c);
	return NULL;
    }

    return c;
}

int cram_byte_array_len_encode(cram_slice *slice, cram_codec *c,
			cram_block *out, char *in, int in_size) {
    return -1; // not imp.
}

void cram_byte_array_len_encode_free(cram_codec *c) {
    if (!c)
	return;
    free(c);
}

int cram_byte_array_len_encode_store(cram_codec *c, cram_block *b,
				     char *prefix, int version) {
    int len = 0;

    if (prefix) {
	size_t l = strlen(prefix);
	BLOCK_APPEND(b, prefix, l);
	len += l;
    }

    len += itf8_put_blk(b, c->codec);
    len += itf8_put_blk(b, c->e_byte_array_len.len_len +
			   c->e_byte_array_len.val_len);
    BLOCK_APPEND(b, c->e_byte_array_len.len_dat, c->e_byte_array_len.len_len);
    len += c->e_byte_array_len.len_len;

    BLOCK_APPEND(b, c->e_byte_array_len.val_dat, c->e_byte_array_len.val_len);
    len += c->e_byte_array_len.val_len;

    return len;
}

cram_codec *cram_byte_array_len_encode_init(cram_stats *st,
					    enum cram_external_type option,
					    void *dat,
					    int version) {
    cram_codec *c;
    cram_byte_array_len_encoder *e = (cram_byte_array_len_encoder *)dat;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec = E_BYTE_ARRAY_LEN;
    c->free = cram_byte_array_len_encode_free;
    c->encode = cram_byte_array_len_encode;
    c->store = cram_byte_array_len_encode_store;

    c->e_byte_array_len.len_len = e->len_len;
    c->e_byte_array_len.len_dat = e->len_dat;
    c->e_byte_array_len.val_len = e->val_len;
    c->e_byte_array_len.val_dat = e->val_dat;

    return c;
}

/*
 * ---------------------------------------------------------------------------
 * BYTE_ARRAY_STOP
 */
int cram_byte_array_stop_decode_char(cram_slice *slice, cram_codec *c,
				     cram_block *in, char *out,
				     int *out_size) {
    int i;
    cram_block *b = NULL;
    char *cp, ch;

    if (slice->block_by_id) {
	if (!(b = slice->block_by_id[c->byte_array_stop.content_id]))
	    return -1;
    } else {
	for (i = 0; i < slice->hdr->num_blocks; i++) {
	    b = slice->block[i];
	    if (b->content_type == EXTERNAL &&
		b->content_id == c->byte_array_stop.content_id) {
		break;
	    }
	}
	if (i == slice->hdr->num_blocks || !b)
	    return -1;
    }

    if (b->idx >= b->uncomp_size)
	return -1;

    cp = (char *)b->data + b->idx;
    while ((ch = *cp) != (char)c->byte_array_stop.stop) {
	if (cp - (char *)b->data >= b->uncomp_size)
	    return -1;
	*out++ = ch;
	cp++;
    }

    *out_size = cp - (char *)(b->data + b->idx);
    b->idx = cp - (char *)b->data + 1;

    return 0;
}

int cram_byte_array_stop_decode_block(cram_slice *slice, cram_codec *c,
				      cram_block *in, char *out_,
				      int *out_size) {
    int space = 256;
    cram_block *b = NULL;
    cram_block *out = (cram_block *)out_;
    char *cp, ch, *out_cp, *cp_end, *out_end;
    char stop;

    if (slice->block_by_id) {
	if (!(b = slice->block_by_id[c->byte_array_stop.content_id]))
	    return -1;
    } else {
	int i;
	for (i = 0; i < slice->hdr->num_blocks; i++) {
	    b = slice->block[i];
	    if (b->content_type == EXTERNAL &&
		b->content_id == c->byte_array_stop.content_id) {
		break;
	    }
	}
	if (i == slice->hdr->num_blocks || !b)
	    return -1;
    }

    if (b->idx >= b->uncomp_size)
	return -1;
    cp = (char *)b->data + b->idx;
    cp_end = (char *)b->data + b->uncomp_size;
    BLOCK_GROW(out, space);
    out_cp = (char *)BLOCK_END(out);
    out_end = out_cp + space;

    stop = c->byte_array_stop.stop;
    while ((ch = *cp) != stop) {
	if (cp++ == cp_end)
	    return -1;
	*out_cp++ = ch;

	if (out_cp == out_end) {
	    BLOCK_SIZE(out) = out_cp - (char *)BLOCK_DATA(out);
	    space *= 2;
	    BLOCK_GROW(out, space);
	    out_cp = (char *)BLOCK_END(out);
	    out_end = out_cp + space;
	}
    }
    BLOCK_SIZE(out) = out_cp - (char *)BLOCK_DATA(out);

    *out_size = cp - (char *)(b->data + b->idx);
    b->idx = cp - (char *)b->data + 1;

    return 0;
}

void cram_byte_array_stop_decode_free(cram_codec *c) {
    if (!c) return;

    free(c);
}

cram_codec *cram_byte_array_stop_decode_init(char *data, int size,
					     enum cram_external_type option,
					     int version) {
    cram_codec *c;
    unsigned char *cp = (unsigned char *)data;

    if (!(c = malloc(sizeof(*c))))
	return NULL;

    c->codec  = E_BYTE_ARRAY_STOP;
    c->decode = (option == E_BYTE_ARRAY_BLOCK)
	? cram_byte_array_stop_decode_block
	: cram_byte_array_stop_decode_char;
    c->free   = cram_byte_array_stop_decode_free;
    
    c->byte_array_stop.stop = *cp++;
    if (version == CRAM_1_VERS) {
	c->byte_array_stop.content_id = cp[0] + (cp[1]<<8) + (cp[2]<<16)
	    + (cp[3]<<24);
	cp += 4;
    } else {
	cp += itf8_get(cp, &c->byte_array_stop.content_id);
    }

    if ((char *)cp - data != size) {
	fprintf(stderr, "Malformed byte_array_stop header stream\n");
	free(c);
	return NULL;
    }

    return c;
}

int cram_byte_array_stop_encode(cram_slice *slice, cram_codec *c,
				cram_block *out, char *in, int in_size) {
    return -1; // not imp.
}

void cram_byte_array_stop_encode_free(cram_codec *c) {
    if (!c)
	return;
    free(c);
}

int cram_byte_array_stop_encode_store(cram_codec *c, cram_block *b,
				      char *prefix, int version) {
    int len = 0;
    char buf[20], *cp = buf;

    if (prefix) {
	size_t l = strlen(prefix);
	BLOCK_APPEND(b, prefix, l);
	len += l;
    }

    cp += itf8_put(cp, c->codec);

    if (version == CRAM_1_VERS) {
	cp += itf8_put(cp, 5);
	*cp++ = c->e_byte_array_stop.stop;
	*cp++ = (c->e_byte_array_stop.content_id >>  0) & 0xff;
	*cp++ = (c->e_byte_array_stop.content_id >>  8) & 0xff;
	*cp++ = (c->e_byte_array_stop.content_id >> 16) & 0xff;
	*cp++ = (c->e_byte_array_stop.content_id >> 24) & 0xff;
    } else {
	cp += itf8_put(cp, 1 + itf8_size(c->e_byte_array_stop.content_id));
	*cp++ = c->e_byte_array_stop.stop;
	cp += itf8_put(cp, c->e_byte_array_stop.content_id);
    }

    BLOCK_APPEND(b, buf, cp-buf);
    len += cp-buf;

    return len;
}

cram_codec *cram_byte_array_stop_encode_init(cram_stats *st,
					     enum cram_external_type option,
					     void *dat,
					     int version) {
    cram_codec *c;

    c = malloc(sizeof(*c));
    if (!c)
	return NULL;
    c->codec = E_BYTE_ARRAY_STOP;
    c->free = cram_byte_array_stop_encode_free;
    c->encode = cram_byte_array_stop_encode;
    c->store = cram_byte_array_stop_encode_store;

    c->e_byte_array_stop.stop = ((int *)dat)[0];
    c->e_byte_array_stop.content_id = ((int *)dat)[1];

    return c;
}

/*
 * ---------------------------------------------------------------------------
 */

char *cram_encoding2str(enum cram_encoding t) {
    switch (t) {
    case E_NULL:            return "NULL";
    case E_EXTERNAL:        return "EXTERNAL";
    case E_GOLOMB:          return "GOLOMB";
    case E_HUFFMAN:         return "HUFFMAN";
    case E_BYTE_ARRAY_LEN:  return "BYTE_ARRAY_LEN";
    case E_BYTE_ARRAY_STOP: return "BYTE_ARRAY_STOP";
    case E_BETA:            return "BETA";
    case E_SUBEXP:          return "SUBEXP";
    case E_GOLOMB_RICE:     return "GOLOMB_RICE";
    case E_GAMMA:           return "GAMMA";
    }
    return "?";
}

static cram_codec *(*decode_init[])(char *data,
				    int size,
				    enum cram_external_type option,
				    int version) = {
    NULL,
    cram_external_decode_init,
    NULL,
    cram_huffman_decode_init,
    cram_byte_array_len_decode_init,
    cram_byte_array_stop_decode_init,
    cram_beta_decode_init,
    cram_subexp_decode_init,
    NULL,
    cram_gamma_decode_init,
};

cram_codec *cram_decoder_init(enum cram_encoding codec,
			      char *data, int size,
			      enum cram_external_type option,
			      int version) {
    if (decode_init[codec]) {
	return decode_init[codec](data, size, option, version);
    } else {
	fprintf(stderr, "Unimplemented codec of type %s\n", codec2str(codec));
	return NULL;
    }
}

static cram_codec *(*encode_init[])(cram_stats *stx,
				    enum cram_external_type option,
				    void *opt,
				    int version) = {
    NULL,
    cram_external_encode_init,
    NULL,
    cram_huffman_encode_init,
    cram_byte_array_len_encode_init,
    cram_byte_array_stop_encode_init,
    cram_beta_encode_init,
    NULL, //cram_subexp_encode_init,
    NULL,
    NULL, //cram_gamma_encode_init,
};

cram_codec *cram_encoder_init(enum cram_encoding codec,
			      cram_stats *st,
			      enum cram_external_type option,
			      void *dat,
			      int version) {
    if (st && !st->nvals)
	return NULL;

    if (encode_init[codec]) {
	return encode_init[codec](st, option, dat, version);
    } else {
	fprintf(stderr, "Unimplemented codec of type %s\n", codec2str(codec));
	abort();
    }
}
