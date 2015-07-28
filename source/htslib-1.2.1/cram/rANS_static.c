/*
 * Copyright (c) 2014 Genome Research Ltd.
 * Author(s): James Bonfield
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *       Institute nor the names of its contributors may be used to endorse
 *       or promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2014
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#ifdef _WIN32
#else
#include <sys/time.h>
#endif


#include "cram/rANS_static.h"
#include "cram/rANS_byte.h"

#define TF_SHIFT 12
#define TOTFREQ (1<<TF_SHIFT)

#define ABS(a) ((a)>0?(a):-(a))
#ifndef BLK_SIZE
#  define BLK_SIZE 1024*1024
#endif

// Room to allow for expanded BLK_SIZE on worst case compression.
#define BLK_SIZE2 ((int)(1.05*BLK_SIZE))

/*-----------------------------------------------------------------------------
 * Memory to memory compression functions.
 *
 * These are original versions without any manual loop unrolling. They
 * are easier to understand, but can be up to 2x slower.
 */

unsigned char *rans_compress_O0(unsigned char *in, unsigned int in_size,
				unsigned int *out_size) {
	unsigned char *out_buf = (unsigned char*)malloc(1.05*in_size + 257 * 257 * 3 + 9);
    unsigned char *cp, *out_end;
    RansEncSymbol syms[256];
    RansState rans0, rans1, rans2, rans3;
    uint8_t* ptr;
    int F[256] = {0}, i, j, tab_size, rle, x, fsum = 0;
    int m = 0, M = 0;
    uint64_t tr;

    if (!out_buf)
	return NULL;

    ptr = out_end = out_buf + (int)(1.05*in_size) + 257*257*3 + 9;

    // Compute statistics
    for (i = 0; i < in_size; i++) {
	F[in[i]]++;
    }
    tr = ((uint64_t)TOTFREQ<<31)/in_size + (1<<30)/in_size;

    // Normalise so T[i] == TOTFREQ
    for (m = M = j = 0; j < 256; j++) {
	if (!F[j])
	    continue;

	if (m < F[j])
	    m = F[j], M = j;

	if ((F[j] = (F[j]*tr)>>31) == 0)
	    F[j] = 1;
	fsum += F[j];
    }

    fsum++;
    if (fsum < TOTFREQ)
	F[M] += TOTFREQ-fsum;
    else
	F[M] -= fsum-TOTFREQ;

    //printf("F[%d]=%d\n", M, F[M]);
    assert(F[M]>0);

    // Encode statistics.
    cp = out_buf+9;

    for (x = rle = j = 0; j < 256; j++) {
	if (F[j]) {
	    // j
	    if (rle) {
		rle--;
	    } else {
		*cp++ = j;
		if (!rle && j && F[j-1])  {
		    for(rle=j+1; rle<256 && F[rle]; rle++)
			;
		    rle -= j+1;
		    *cp++ = rle;
		}
		//fprintf(stderr, "%d: %d %d\n", j, rle, N[j]);
	    }

	    // F[j]
	    if (F[j]<128) {
		*cp++ = F[j];
	    } else {
		*cp++ = 128 | (F[j]>>8);
		*cp++ = F[j]&0xff;
	    }
	    RansEncSymbolInit(&syms[j], x, F[j], TF_SHIFT);
	    x += F[j];
	}
    }
    *cp++ = 0;

    //write(1, out_buf+4, cp-(out_buf+4));
    tab_size = cp-out_buf;

    RansEncInit(&rans0);
    RansEncInit(&rans1);
    RansEncInit(&rans2);
    RansEncInit(&rans3);

    switch (i=(in_size&3)) {
    case 3: RansEncPutSymbol(&rans2, &ptr, &syms[in[in_size-(i-2)]]);
    case 2: RansEncPutSymbol(&rans1, &ptr, &syms[in[in_size-(i-1)]]);
    case 1: RansEncPutSymbol(&rans0, &ptr, &syms[in[in_size-(i-0)]]);
    case 0:
	break;
    }
    for (i=(in_size &~3); i>0; i-=4) {
	RansEncSymbol *s3 = &syms[in[i-1]];
	RansEncSymbol *s2 = &syms[in[i-2]];
	RansEncSymbol *s1 = &syms[in[i-3]];
	RansEncSymbol *s0 = &syms[in[i-4]];

	RansEncPutSymbol(&rans3, &ptr, s3);
	RansEncPutSymbol(&rans2, &ptr, s2);
	RansEncPutSymbol(&rans1, &ptr, s1);
	RansEncPutSymbol(&rans0, &ptr, s0);
    }

    RansEncFlush(&rans3, &ptr);
    RansEncFlush(&rans2, &ptr);
    RansEncFlush(&rans1, &ptr);
    RansEncFlush(&rans0, &ptr);

    // Finalise block size and return it
    *out_size = (out_end - ptr) + tab_size;

    cp = out_buf;

    *cp++ = 0; // order
    *cp++ = ((*out_size-9)>> 0) & 0xff;
    *cp++ = ((*out_size-9)>> 8) & 0xff;
    *cp++ = ((*out_size-9)>>16) & 0xff;
    *cp++ = ((*out_size-9)>>24) & 0xff;

    *cp++ = (in_size>> 0) & 0xff;
    *cp++ = (in_size>> 8) & 0xff;
    *cp++ = (in_size>>16) & 0xff;
    *cp++ = (in_size>>24) & 0xff;

    memmove(out_buf + tab_size, ptr, out_end-ptr);

    return out_buf;
}

typedef struct {
    struct {
	int F;
	int C;
    } fc[256];
    unsigned char *R;
} ari_decoder;

unsigned char *rans_uncompress_O0(unsigned char *in, unsigned int in_size,
				  unsigned int *out_size) {
    /* Load in the static tables */
    unsigned char *cp = in + 9;
    int i, j, x, out_sz, in_sz, rle;
    char *out_buf;
    ari_decoder D;
    RansDecSymbol syms[256];

    memset(&D, 0, sizeof(D));

    if (*in++ != 0) // Order-0 check
	return NULL;
    
    in_sz  = ((in[0])<<0) | ((in[1])<<8) | ((in[2])<<16) | ((in[3])<<24);
    out_sz = ((in[4])<<0) | ((in[5])<<8) | ((in[6])<<16) | ((in[7])<<24);
    if (in_sz != in_size-9)
	return NULL;

	out_buf = (char*)malloc(out_sz);
    if (!out_buf)
	return NULL;

    //fprintf(stderr, "out_sz=%d\n", out_sz);

    // Precompute reverse lookup of frequency.
    rle = x = 0;
    j = *cp++;
    do {
	if ((D.fc[j].F = *cp++) >= 128) {
	    D.fc[j].F &= ~128;
	    D.fc[j].F = ((D.fc[j].F & 127) << 8) | *cp++;
	}
	D.fc[j].C = x;

	RansDecSymbolInit(&syms[j], D.fc[j].C, D.fc[j].F);

	/* Build reverse lookup table */
	if (!D.R) D.R = (unsigned char *)malloc(TOTFREQ);
	memset(&D.R[x], j, D.fc[j].F);

	x += D.fc[j].F;

	if (!rle && j+1 == *cp) {
	    j = *cp++;
	    rle = *cp++;
	} else if (rle) {
	    rle--;
	    j++;
	} else {
	    j = *cp++;
	}
    } while(j);

    assert(x < TOTFREQ);

    RansState rans0, rans1, rans2, rans3;
    uint8_t *ptr = cp;
    RansDecInit(&rans0, &ptr);
    RansDecInit(&rans1, &ptr);
    RansDecInit(&rans2, &ptr);
    RansDecInit(&rans3, &ptr);

    int out_end = (out_sz&~3);

    RansState R[4];
    R[0] = rans0;
    R[1] = rans1;
    R[2] = rans2;
    R[3] = rans3;
    uint32_t mask = (1u << TF_SHIFT)-1;

    for (i=0; i < out_end; i+=4) {
	uint32_t m[4] = {R[0] & mask,
			 R[1] & mask,
			 R[2] & mask,
			 R[3] & mask};
	uint8_t c[4] = {D.R[m[0]],
			D.R[m[1]],
			D.R[m[2]],
			D.R[m[3]]};
	out_buf[i+0] = c[0];
	out_buf[i+1] = c[1];
	out_buf[i+2] = c[2];
	out_buf[i+3] = c[3];

	// RansDecAdvanceSymbolStep(&R[0], &syms[c[0]], TF_SHIFT);
	// RansDecAdvanceSymbolStep(&R[1], &syms[c[1]], TF_SHIFT);
	// RansDecAdvanceSymbolStep(&R[2], &syms[c[2]], TF_SHIFT);
	// RansDecAdvanceSymbolStep(&R[3], &syms[c[3]], TF_SHIFT);
	R[0] = syms[c[0]].freq * (R[0]>>TF_SHIFT);
	R[1] = syms[c[1]].freq * (R[1]>>TF_SHIFT);
	R[2] = syms[c[2]].freq * (R[2]>>TF_SHIFT);
	R[3] = syms[c[3]].freq * (R[3]>>TF_SHIFT);

	R[0] += m[0] - syms[c[0]].start;
	R[1] += m[1] - syms[c[1]].start;
	R[2] += m[2] - syms[c[2]].start;
	R[3] += m[3] - syms[c[3]].start;

	RansDecRenorm(&R[0], &ptr);
	RansDecRenorm(&R[1], &ptr);
	RansDecRenorm(&R[2], &ptr);
	RansDecRenorm(&R[3], &ptr);
    }

    rans0 = R[0];
    rans1 = R[1];
    rans2 = R[2];
    rans3 = R[3];

    switch(out_sz&3) {
	unsigned char c;
    case 0:
	break;
    case 1:
	c = D.R[RansDecGet(&rans0, TF_SHIFT)];
	RansDecAdvanceSymbol(&rans0, &ptr, &syms[c], TF_SHIFT);
	out_buf[out_end] = c;
	break;

    case 2:
	c = D.R[RansDecGet(&rans0, TF_SHIFT)];
	RansDecAdvanceSymbol(&rans0, &ptr, &syms[c], TF_SHIFT);
	out_buf[out_end] = c;

	c = D.R[RansDecGet(&rans1, TF_SHIFT)];
	RansDecAdvanceSymbol(&rans1, &ptr, &syms[c], TF_SHIFT);
	out_buf[out_end+1] = c;
	break;

    case 3:
	c = D.R[RansDecGet(&rans0, TF_SHIFT)];
	RansDecAdvanceSymbol(&rans0, &ptr, &syms[c], TF_SHIFT);
	out_buf[out_end] = c;

	c = D.R[RansDecGet(&rans1, TF_SHIFT)];
	RansDecAdvanceSymbol(&rans1, &ptr, &syms[c], TF_SHIFT);
	out_buf[out_end+1] = c;

	c = D.R[RansDecGet(&rans2, TF_SHIFT)];
	RansDecAdvanceSymbol(&rans2, &ptr, &syms[c], TF_SHIFT);
	out_buf[out_end+2] = c;
	break;
    }

    *out_size = out_sz;

    if (D.R) free(D.R);

    return (unsigned char *)out_buf;
}

unsigned char *rans_compress_O1(unsigned char *in, unsigned int in_size,
				unsigned int *out_size) {
    unsigned char *out_buf, *out_end, *cp;
    unsigned int last_i, tab_size, rle_i, rle_j;
    RansEncSymbol syms[256][256];

    if (in_size < 4)
	return rans_compress_O0(in, in_size, out_size);

	out_buf = (unsigned char*)malloc(1.05*in_size + 257 * 257 * 3 + 9);
    if (!out_buf)
	return NULL;

    out_end = out_buf + (int)(1.05*in_size) + 257*257*3 + 9;
    cp = out_buf+9;

    int F[256][256], T[256], i, j;
    unsigned char c;

    memset(F, 0, 256*256*sizeof(int));
    memset(T, 0, 256*sizeof(int));
    //for (last = 0, i=in_size-1; i>=0; i--) {
    //	F[last][c = in[i]]++;
    //	T[last]++;
    //	last = c;
    //}

    for (last_i=i=0; i<in_size; i++) {
	F[last_i][c = in[i]]++;
	T[last_i]++;
	last_i = c;
    }
    F[0][in[1*(in_size>>2)]]++;
    F[0][in[2*(in_size>>2)]]++;
    F[0][in[3*(in_size>>2)]]++;
    T[0]+=3;

    // Normalise so T[i] == TOTFREQ
    for (rle_i = i = 0; i < 256; i++) {
	int t2, m, M;
	unsigned int x;

	if (T[i] == 0)
	    continue;

	//uint64_t p = (TOTFREQ * TOTFREQ) / t;
	double p = ((double)TOTFREQ)/T[i];
	for (t2 = m = M = j = 0; j < 256; j++) {
	    if (!F[i][j])
		continue;

	    if (m < F[i][j])
		m = F[i][j], M = j;

	    //if ((F[i][j] = (F[i][j] * p) / TOTFREQ) == 0)
	    if ((F[i][j] *= p) == 0)
		F[i][j] = 1;
	    t2 += F[i][j];
	}

	t2++;
	if (t2 < TOTFREQ)
	    F[i][M] += TOTFREQ-t2;
	else
	    F[i][M] -= t2-TOTFREQ;

	// Store frequency table
	// i
	if (rle_i) {
	    rle_i--;
	} else {
	    *cp++ = i;
	    // FIXME: could use order-0 statistics to observe which alphabet
	    // symbols are present and base RLE on that ordering instead.
	    if (i && T[i-1]) {
		for(rle_i=i+1; rle_i<256 && T[rle_i]; rle_i++)
		    ;
		rle_i -= i+1;
		*cp++ = rle_i;
	    }
	}

	int *F_i_ = F[i];
	x = 0;
	rle_j = 0;
	for (j = 0; j < 256; j++) {
	    if (F_i_[j]) {
		//fprintf(stderr, "F[%d][%d]=%d, x=%d\n", i, j, F_i_[j], x);

		// j
		if (rle_j) {
		    rle_j--;
		} else {
		    *cp++ = j;
		    if (!rle_j && j && F_i_[j-1]) {
			for(rle_j=j+1; rle_j<256 && F_i_[rle_j]; rle_j++)
			    ;
			rle_j -= j+1;
			*cp++ = rle_j;
		    }
		}

		// F_i_[j]
		if (F_i_[j]<128) {
		    *cp++ = F_i_[j];
		} else {
		    *cp++ = 128 | (F_i_[j]>>8);
		    *cp++ = F_i_[j]&0xff;
		}

		RansEncSymbolInit(&syms[i][j], x, F_i_[j], TF_SHIFT);
		x += F_i_[j];
	    }
	}
	*cp++ = 0;
    }
    *cp++ = 0;

    //write(1, out_buf+4, cp-(out_buf+4));
    tab_size = cp - out_buf;
    assert(tab_size < 257*257*3);

    RansState rans0, rans1, rans2, rans3;
    RansEncInit(&rans0);
    RansEncInit(&rans1);
    RansEncInit(&rans2);
    RansEncInit(&rans3);

    uint8_t* ptr = out_end;

    int isz4 = in_size>>2;
    int i0 = 1*isz4-2;
    int i1 = 2*isz4-2;
    int i2 = 3*isz4-2;
    int i3 = 4*isz4-2;

    unsigned char l0 = in[i0+1];
    unsigned char l1 = in[i1+1];
    unsigned char l2 = in[i2+1];
    unsigned char l3 = in[i3+1];

    // Deal with the remainder
    l3 = in[in_size-1];
    for (i3 = in_size-2; i3 > 4*isz4-2; i3--) {
	unsigned char c3 = in[i3];
	RansEncPutSymbol(&rans3, &ptr, &syms[c3][l3]);
	l3 = c3;
    }

    for (; i0 >= 0; i0--, i1--, i2--, i3--) {
	unsigned char c0, c1, c2, c3;
	RansEncSymbol *s3 = &syms[c3 = in[i3]][l3];
	RansEncSymbol *s2 = &syms[c2 = in[i2]][l2];
	RansEncSymbol *s1 = &syms[c1 = in[i1]][l1];
	RansEncSymbol *s0 = &syms[c0 = in[i0]][l0];

	RansEncPutSymbol(&rans3, &ptr, s3);
	RansEncPutSymbol(&rans2, &ptr, s2);
	RansEncPutSymbol(&rans1, &ptr, s1);
	RansEncPutSymbol(&rans0, &ptr, s0);

	l0 = c0;
	l1 = c1;
	l2 = c2;
	l3 = c3;
    }

    RansEncPutSymbol(&rans3, &ptr, &syms[0][l3]);
    RansEncPutSymbol(&rans2, &ptr, &syms[0][l2]);
    RansEncPutSymbol(&rans1, &ptr, &syms[0][l1]);
    RansEncPutSymbol(&rans0, &ptr, &syms[0][l0]);

    RansEncFlush(&rans3, &ptr);
    RansEncFlush(&rans2, &ptr);
    RansEncFlush(&rans1, &ptr);
    RansEncFlush(&rans0, &ptr);

    *out_size = (out_end - ptr) + tab_size;

    cp = out_buf;
    *cp++ = 1; // order

    *cp++ = ((*out_size-9)>> 0) & 0xff;
    *cp++ = ((*out_size-9)>> 8) & 0xff;
    *cp++ = ((*out_size-9)>>16) & 0xff;
    *cp++ = ((*out_size-9)>>24) & 0xff;

    *cp++ = (in_size>> 0) & 0xff;
    *cp++ = (in_size>> 8) & 0xff;
    *cp++ = (in_size>>16) & 0xff;
    *cp++ = (in_size>>24) & 0xff;

    memmove(out_buf + tab_size, ptr, out_end-ptr);

    return out_buf;
}

unsigned char *rans_uncompress_O1(unsigned char *in, unsigned int in_size,
				  unsigned int *out_size) {
    /* Load in the static tables */
    unsigned char *cp = in + 9;
    int i, j = -999, x, out_sz, in_sz, rle_i, rle_j;
    char *out_buf;
    ari_decoder D[256];
    RansDecSymbol syms[256][256];

    memset(D, 0, 256*sizeof(*D));

    if (*in++ != 1) // Order-1 check
	return NULL;

    in_sz  = ((in[0])<<0) | ((in[1])<<8) | ((in[2])<<16) | ((in[3])<<24);
    out_sz = ((in[4])<<0) | ((in[5])<<8) | ((in[6])<<16) | ((in[7])<<24);
    if (in_sz != in_size-9)
	return NULL;

	out_buf = (char*)malloc(out_sz);
    if (!out_buf)
	return NULL;

    //fprintf(stderr, "out_sz=%d\n", out_sz);

    //i = *cp++;
    rle_i = 0;
    i = *cp++;
    do {
	rle_j = x = 0;
	j = *cp++;
	do {
	    if ((D[i].fc[j].F = *cp++) >= 128) {
		D[i].fc[j].F &= ~128;
		D[i].fc[j].F = ((D[i].fc[j].F & 127) << 8) | *cp++;
	    }
	    D[i].fc[j].C = x;

	    //fprintf(stderr, "i=%d j=%d F=%d C=%d\n", i, j, D[i].fc[j].F, D[i].fc[j].C);

	    if (!D[i].fc[j].F)
		D[i].fc[j].F = TOTFREQ;

	    RansDecSymbolInit(&syms[i][j], D[i].fc[j].C, D[i].fc[j].F);

	    /* Build reverse lookup table */
	    if (!D[i].R) D[i].R = (unsigned char *)malloc(TOTFREQ);
	    memset(&D[i].R[x], j, D[i].fc[j].F);

	    x += D[i].fc[j].F;
	    assert(x <= TOTFREQ);

	    if (!rle_j && j+1 == *cp) {
		j = *cp++;
		rle_j = *cp++;
	    } else if (rle_j) {
		rle_j--;
		j++;
	    } else {
		j = *cp++;
	    }
	} while(j);

	if (!rle_i && i+1 == *cp) {
	    i = *cp++;
	    rle_i = *cp++;
	} else if (rle_i) {
	    rle_i--;
	    i++;
	} else {
	    i = *cp++;
	}
    } while (i);

    // Precompute reverse lookup of frequency.

    RansState rans0, rans1, rans2, rans3;
    uint8_t *ptr = cp;
    RansDecInit(&rans0, &ptr);
    RansDecInit(&rans1, &ptr);
    RansDecInit(&rans2, &ptr);
    RansDecInit(&rans3, &ptr);

    int isz4 = out_sz>>2;
    int l0 = 0;
    int l1 = 0;
    int l2 = 0;
    int l3 = 0;
    int i4[] = {0*isz4, 1*isz4, 2*isz4, 3*isz4};

    RansState R[4];
    R[0] = rans0;
    R[1] = rans1;
    R[2] = rans2;
    R[3] = rans3;

    for (; i4[0] < isz4; i4[0]++, i4[1]++, i4[2]++, i4[3]++) {
	uint32_t m[4] = {R[0] & ((1u << TF_SHIFT)-1),
			 R[1] & ((1u << TF_SHIFT)-1),
			 R[2] & ((1u << TF_SHIFT)-1),
			 R[3] & ((1u << TF_SHIFT)-1)};

	uint8_t c[4] = {D[l0].R[m[0]],
			D[l1].R[m[1]],
			D[l2].R[m[2]],
			D[l3].R[m[3]]};

	out_buf[i4[0]] = c[0];
	out_buf[i4[1]] = c[1];
	out_buf[i4[2]] = c[2];
	out_buf[i4[3]] = c[3];

	//RansDecAdvanceSymbolStep(&R[0], &syms[l0][c[0]], TF_SHIFT);
	//RansDecAdvanceSymbolStep(&R[1], &syms[l1][c[1]], TF_SHIFT);
	//RansDecAdvanceSymbolStep(&R[2], &syms[l2][c[2]], TF_SHIFT);
	//RansDecAdvanceSymbolStep(&R[3], &syms[l3][c[3]], TF_SHIFT);

	R[0] = syms[l0][c[0]].freq * (R[0]>>TF_SHIFT);
	R[1] = syms[l1][c[1]].freq * (R[1]>>TF_SHIFT);
	R[2] = syms[l2][c[2]].freq * (R[2]>>TF_SHIFT);
	R[3] = syms[l3][c[3]].freq * (R[3]>>TF_SHIFT);

	R[0] += m[0] - syms[l0][c[0]].start;
	R[1] += m[1] - syms[l1][c[1]].start;
	R[2] += m[2] - syms[l2][c[2]].start;
	R[3] += m[3] - syms[l3][c[3]].start;

	RansDecRenorm(&R[0], &ptr);
	RansDecRenorm(&R[1], &ptr);
	RansDecRenorm(&R[2], &ptr);
	RansDecRenorm(&R[3], &ptr);

	l0 = c[0];
	l1 = c[1];
	l2 = c[2];
	l3 = c[3];
    }

    rans0 = R[0];
    rans1 = R[1];
    rans2 = R[2];
    rans3 = R[3];

    // Remainder
    for (; i4[3] < out_sz; i4[3]++) {
	unsigned char c3 = D[l3].R[RansDecGet(&rans3, TF_SHIFT)];
	out_buf[i4[3]] = c3;
	RansDecAdvanceSymbol(&rans3, &ptr, &syms[l3][c3], TF_SHIFT);
	l3 = c3;
    }

    *out_size = out_sz;

    for (i = 0; i < 256; i++)
	if (D[i].R) free(D[i].R);

    return (unsigned char *)out_buf;
}

/*-----------------------------------------------------------------------------
 * Simple interface to the order-0 vs order-1 encoders and decoders.
 */
unsigned char *rans_compress(unsigned char *in, unsigned int in_size,
			     unsigned int *out_size, int order) {
    return order
	? rans_compress_O1(in, in_size, out_size)
	: rans_compress_O0(in, in_size, out_size);
}

unsigned char *rans_uncompress(unsigned char *in, unsigned int in_size,
			       unsigned int *out_size) {
    return in[0]
	? rans_uncompress_O1(in, in_size, out_size)
	: rans_uncompress_O0(in, in_size, out_size);
}


#ifdef TEST_MAIN
/*-----------------------------------------------------------------------------
 * Main.
 *
 * This is a simple command line tool for testing order-0 and order-1
 * compression using the rANS codec. Simply compile with
 *
 * gcc -DTEST_MAIN -O3 -I. cram/rANS_static.c -o cram/rANS_static
 *
 * Usage: cram/rANS_static -o0 < file    > file.o0
 *        cram/rANS_static -d  < file.o0 > file2
 *
 *        cram/rANS_static -o1 < file    > file.o1
 *        cram/rANS_static -d  < file.o1 > file2
 */
int main(int argc, char **argv) {
    int opt, order = 0;
    unsigned char in_buf[BLK_SIZE2+257*257*3];
    int decode = 0;
    FILE *infp = stdin, *outfp = stdout;
    struct timeval tv1, tv2;
    size_t bytes = 0;

    extern char *optarg;
    extern int optind;

    while ((opt = getopt(argc, argv, "o:d")) != -1) {
	switch (opt) {
	case 'o':
	    order = atoi(optarg);
	    break;

	case 'd':
	    decode = 1;
	    break;
	}
    }

    order = order ? 1 : 0; // Only support O(0) and O(1)

    if (optind < argc) {
	if (!(infp = fopen(argv[optind], "rb"))) {
	    perror(argv[optind]);
	    return 1;
	}
	optind++;
    }

    if (optind < argc) {
	if (!(outfp = fopen(argv[optind], "wb"))) {
	    perror(argv[optind]);
	    return 1;
	}
	optind++;
    }

    gettimeofday(&tv1, NULL);

    if (decode) {
	// Only used in some test implementations of RC_GetFreq()
	//RC_init();
	//RC_init2();

	for (;;) {
	    uint32_t in_size, out_size;
	    unsigned char *out;

	    if (4 != fread(&in_size, 1, 4, infp))
		break;
	    if (in_size != fread(in_buf, 1, in_size, infp)) {
		fprintf(stderr, "Truncated input\n");
		exit(1);
	    }
	    out = rans_uncompress(in_buf, in_size, &out_size);
	    if (!out)
		abort();

	    fwrite(out, 1, out_size, outfp);
	    free(out);

	    bytes += out_size;
	}
    } else {
	for (;;) {
	    uint32_t in_size, out_size;
	    unsigned char *out;

	    in_size = fread(in_buf, 1, BLK_SIZE, infp);
	    if (in_size <= 0)
		break;

	    out = rans_compress(in_buf, in_size, &out_size, order);

	    fwrite(&out_size, 1, 4, outfp);
	    fwrite(out, 1, out_size, outfp);
	    free(out);

	    bytes += in_size;
	}
    }

    gettimeofday(&tv2, NULL);

    fprintf(stderr, "Took %ld microseconds, %5.1f MB/s\n",
	    (long)(tv2.tv_sec - tv1.tv_sec)*1000000 +
	    tv2.tv_usec - tv1.tv_usec,
	    (double)bytes / ((long)(tv2.tv_sec - tv1.tv_sec)*1000000 +
			     tv2.tv_usec - tv1.tv_usec));
    return 0;
}
#endif
