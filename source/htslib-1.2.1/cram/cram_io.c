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
 * CRAM I/O primitives.
 *
 * - ITF8 encoding and decoding.
 * - Block based I/O
 * - Zlib inflating and deflating (memory)
 * - CRAM basic data structure reading and writing
 * - File opening / closing
 * - Reference sequence handling
 */

/*
 * TODO: BLOCK_GROW, BLOCK_RESIZE, BLOCK_APPEND and itf8_put_blk all need
 * a way to return errors for when malloc fails.
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
#ifdef HAVE_LIBBZ2
#include <bzlib.h>
#endif
#ifdef HAVE_LIBLZMA
#include <lzma.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <ctype.h>

#include "cram/cram.h"
#include "cram/os.h"
#include "cram/md5.h"
#include "cram/open_trace_file.h"
#include "cram/rANS_static.h"
#ifdef _WIN32
#include <stdio.h>
#endif

//#define REF_DEBUG

#ifdef REF_DEBUG
#include <sys/syscall.h>
#define gettid() (int)syscall(SYS_gettid)

#define RP(...) fprintf (stderr, __VA_ARGS__)
#else
#define RP(...) 
#endif

#include "htslib/hfile.h"
#include "htslib/bgzf.h"
#include "htslib/faidx.h"

#define TRIAL_SPAN 50
#define NTRIALS 3


/* ----------------------------------------------------------------------
 * ITF8 encoding and decoding.
 *
* Also see the itf8_get and itf8_put macros in cram_io.h
 */

/*
 * Reads an integer in ITF-8 encoding from 'cp' and stores it in
 * *val.
 *
 * Returns the number of bytes read on success
 *        -1 on failure
 */
int itf8_decode(cram_fd *fd, int32_t *val_p) {
    static int nbytes[16] = {
	0,0,0,0, 0,0,0,0,                               // 0000xxxx - 0111xxxx
	1,1,1,1,                                        // 1000xxxx - 1011xxxx
	2,2,                                            // 1100xxxx - 1101xxxx
	3,                                              // 1110xxxx
	4,                                              // 1111xxxx
    };

    static int nbits[16] = {
	0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, // 0000xxxx - 0111xxxx
	0x3f, 0x3f, 0x3f, 0x3f,                         // 1000xxxx - 1011xxxx
	0x1f, 0x1f,                                     // 1100xxxx - 1101xxxx
	0x0f,                                           // 1110xxxx
	0x0f,                                           // 1111xxxx
    };

    int32_t val = hgetc(fd->fp);
    if (val == -1)
	return -1;

    int i = nbytes[val>>4];
    val &= nbits[val>>4];

    switch(i) {
    case 0:
	*val_p = val;
	return 1;

    case 1:
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	*val_p = val;
	return 2;

    case 2:
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	*val_p = val;
	return 3;

    case 3:
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	*val_p = val;
	return 4;

    case 4: // really 3.5 more, why make it different?
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<4) | (((unsigned char)hgetc(fd->fp)) & 0x0f);
	*val_p = val;
    }

    return 5;
}

/*
 * Encodes and writes a single integer in ITF-8 format.
 * Returns 0 on success
 *        -1 on failure
 */
int itf8_encode(cram_fd *fd, int32_t val) {
    char buf[5];
    int len = itf8_put(buf, val);
    return hwrite(fd->fp, buf, len) == len ? 0 : -1;
}

#ifndef ITF8_MACROS
/*
 * As above, but decoding from memory
 */
int itf8_get(char *cp, int32_t *val_p) {
    unsigned char *up = (unsigned char *)cp;
    
    if (up[0] < 0x80) {
	*val_p =   up[0];
	return 1;
    } else if (up[0] < 0xc0) {
	*val_p = ((up[0] <<8) |  up[1])                           & 0x3fff;
	return 2;
    } else if (up[0] < 0xe0) {
	*val_p = ((up[0]<<16) | (up[1]<< 8) |  up[2])             & 0x1fffff;
	return 3;
    } else if (up[0] < 0xf0) {
	*val_p = ((up[0]<<24) | (up[1]<<16) | (up[2]<<8) | up[3]) & 0x0fffffff;
	return 4;
    } else {
	*val_p = ((up[0] & 0x0f)<<28) | (up[1]<<20) | (up[2]<<12) | (up[3]<<4) | (up[4] & 0x0f);
	return 5;
    }
}

/*
 * Stores a value to memory in ITF-8 format.
 *
 * Returns the number of bytes required to store the number.
 * This is a maximum of 5 bytes.
 */
int itf8_put(char *cp, int32_t val) {
    if        (!(val & ~0x00000007f)) { // 1 byte
	*cp = val;
	return 1;
    } else if (!(val & ~0x00003fff)) { // 2 byte
	*cp++ = (val >> 8 ) | 0x80;
	*cp   = val & 0xff;
	return 2;
    } else if (!(val & ~0x01fffff)) { // 3 byte
	*cp++ = (val >> 16) | 0xc0;
	*cp++ = (val >> 8 ) & 0xff;
	*cp   = val & 0xff;
	return 3;
    } else if (!(val & ~0x0fffffff)) { // 4 byte
	*cp++ = (val >> 24) | 0xe0;
	*cp++ = (val >> 16) & 0xff;
	*cp++ = (val >> 8 ) & 0xff;
	*cp   = val & 0xff;
	return 4;
    } else {                           // 5 byte
	*cp++ = 0xf0 | ((val>>28) & 0xff);
	*cp++ = (val >> 20) & 0xff;
	*cp++ = (val >> 12) & 0xff;
	*cp++ = (val >> 4 ) & 0xff;
	*cp = val & 0x0f;
	return 5;
    }
}
#endif

/* 64-bit itf8 variant */
int ltf8_put(char *cp, int64_t val) {
    if        (!(val & ~((1LL<<7)-1))) {
	*cp = val;
	return 1;
    } else if (!(val & ~((1LL<<(6+8))-1))) {
	*cp++ = (val >> 8 ) | 0x80;
	*cp   = val & 0xff;
	return 2;
    } else if (!(val & ~((1LL<<(5+2*8))-1))) {
	*cp++ = (val >> 16) | 0xc0;
	*cp++ = (val >> 8 ) & 0xff;
	*cp   = val & 0xff;
	return 3;
    } else if (!(val & ~((1LL<<(4+3*8))-1))) {
	*cp++ = (val >> 24) | 0xe0;
	*cp++ = (val >> 16) & 0xff;
	*cp++ = (val >> 8 ) & 0xff;
	*cp   = val & 0xff;
	return 4;
    } else if (!(val & ~((1LL<<(3+4*8))-1))) {
	*cp++ = (val >> 32) | 0xf0;
	*cp++ = (val >> 24) & 0xff;
	*cp++ = (val >> 16) & 0xff;
	*cp++ = (val >> 8 ) & 0xff;
	*cp   = val & 0xff;
	return 5;
    } else if (!(val & ~((1LL<<(2+5*8))-1))) {
	*cp++ = (val >> 40) | 0xf8;
	*cp++ = (val >> 32) & 0xff;
	*cp++ = (val >> 24) & 0xff;
	*cp++ = (val >> 16) & 0xff;
	*cp++ = (val >> 8 ) & 0xff;
	*cp   = val & 0xff;
	return 6;
    } else if (!(val & ~((1LL<<(1+6*8))-1))) {
	*cp++ = (val >> 48) | 0xfc;
	*cp++ = (val >> 40) & 0xff;
	*cp++ = (val >> 32) & 0xff;
	*cp++ = (val >> 24) & 0xff;
	*cp++ = (val >> 16) & 0xff;
	*cp++ = (val >> 8 ) & 0xff;
	*cp   = val & 0xff;
	return 7;
    } else if (!(val & ~((1LL<<(7*8))-1))) {
	*cp++ = (val >> 56) | 0xfe;
	*cp++ = (val >> 48) & 0xff;
	*cp++ = (val >> 40) & 0xff;
	*cp++ = (val >> 32) & 0xff;
	*cp++ = (val >> 24) & 0xff;
	*cp++ = (val >> 16) & 0xff;
	*cp++ = (val >> 8 ) & 0xff;
	*cp   = val & 0xff;
	return 8;
    } else {
	*cp++ = 0xff;
	*cp++ = (val >> 56) & 0xff;
	*cp++ = (val >> 48) & 0xff;
	*cp++ = (val >> 40) & 0xff;
	*cp++ = (val >> 32) & 0xff;
	*cp++ = (val >> 24) & 0xff;
	*cp++ = (val >> 16) & 0xff;
	*cp++ = (val >> 8 ) & 0xff;
	*cp   = val & 0xff;
	return 9;
    }
}

int ltf8_get(char *cp, int64_t *val_p) {
    unsigned char *up = (unsigned char *)cp;
    
    if (up[0] < 0x80) {
	*val_p =   up[0];
	return 1;
    } else if (up[0] < 0xc0) {
	*val_p = (((uint64_t)up[0]<< 8) |
		   (uint64_t)up[1]) & (((1LL<<(6+8)))-1);
	return 2;
    } else if (up[0] < 0xe0) {
	*val_p = (((uint64_t)up[0]<<16) |
		  ((uint64_t)up[1]<< 8) |
		   (uint64_t)up[2]) & ((1LL<<(5+2*8))-1);
	return 3;
    } else if (up[0] < 0xf0) {
	*val_p = (((uint64_t)up[0]<<24) |
		  ((uint64_t)up[1]<<16) |
		  ((uint64_t)up[2]<< 8) |
		   (uint64_t)up[3]) & ((1LL<<(4+3*8))-1);
	return 4;
    } else if (up[0] < 0xf8) {
	*val_p = (((uint64_t)up[0]<<32) |
		  ((uint64_t)up[1]<<24) |
		  ((uint64_t)up[2]<<16) |
		  ((uint64_t)up[3]<< 8) |
		   (uint64_t)up[4]) & ((1LL<<(3+4*8))-1);
	return 5;
    } else if (up[0] < 0xfc) {
	*val_p = (((uint64_t)up[0]<<40) |
		  ((uint64_t)up[1]<<32) |
		  ((uint64_t)up[2]<<24) |
		  ((uint64_t)up[3]<<16) |
		  ((uint64_t)up[4]<< 8) |
		   (uint64_t)up[5]) & ((1LL<<(2+5*8))-1);
	return 6;
    } else if (up[0] < 0xfe) {
	*val_p = (((uint64_t)up[0]<<48) |
		  ((uint64_t)up[1]<<40) |
		  ((uint64_t)up[2]<<32) |
		  ((uint64_t)up[3]<<24) |
		  ((uint64_t)up[4]<<16) |
		  ((uint64_t)up[5]<< 8) |
		   (uint64_t)up[6]) & ((1LL<<(1+6*8))-1);
	return 7;
    } else if (up[0] < 0xff) {
	*val_p = (((uint64_t)up[1]<<48) |
		  ((uint64_t)up[2]<<40) |
		  ((uint64_t)up[3]<<32) |
		  ((uint64_t)up[4]<<24) |
		  ((uint64_t)up[5]<<16) |
		  ((uint64_t)up[6]<< 8) |
		   (uint64_t)up[7]) & ((1LL<<(7*8))-1);
	return 8;
    } else {
	*val_p = (((uint64_t)up[1]<<56) |
		  ((uint64_t)up[2]<<48) |
		  ((uint64_t)up[3]<<40) |
		  ((uint64_t)up[4]<<32) |
		  ((uint64_t)up[5]<<24) |
		  ((uint64_t)up[6]<<16) |
		  ((uint64_t)up[7]<< 8) |
		   (uint64_t)up[8]);
	return 9;
    }
}

int ltf8_decode(cram_fd *fd, int64_t *val_p) {
    int c = hgetc(fd->fp);
    int64_t val = (unsigned char)c;
    if (c == -1)
	return -1;

    if (val < 0x80) {
	*val_p =   val;
	return 1;

    } else if (val < 0xc0) {
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	*val_p = val & (((1LL<<(6+8)))-1);
	return 2;

    } else if (val < 0xe0) {
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	*val_p = val & ((1LL<<(5+2*8))-1);
	return 3;

    } else if (val < 0xf0) {
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	*val_p = val & ((1LL<<(4+3*8))-1);
	return 4;

    } else if (val < 0xf8) {
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	*val_p = val & ((1LL<<(3+4*8))-1);
	return 5;

    } else if (val < 0xfc) {
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	*val_p = val & ((1LL<<(2+5*8))-1);
	return 6;

    } else if (val < 0xfe) {
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	*val_p = val & ((1LL<<(1+6*8))-1);
	return 7;

    } else if (val < 0xff) {
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	*val_p = val & ((1LL<<(7*8))-1);
	return 8;

    } else {
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	val = (val<<8) | (unsigned char)hgetc(fd->fp);
	*val_p = val;
    }

    return 9;
}

/*
 * Pushes a value in ITF8 format onto the end of a block.
 * This shouldn't be used for high-volume data as it is not the fastest
 * method.
 *
 * Returns the number of bytes written
 */
int itf8_put_blk(cram_block *blk, int val) {
    char buf[5];
    int sz;

    sz = itf8_put(buf, val);
    BLOCK_APPEND(blk, buf, sz);
    return sz;
}

/*
 * Decodes a 32-bit little endian value from fd and stores in val.
 *
 * Returns the number of bytes read on success
 *         -1 on failure
 */
int int32_decode(cram_fd *fd, int32_t *val) {
    int32_t i;
    if (4 != hread(fd->fp, &i, 4))
	return -1;

    *val = le_int4(i);
    return 4;
}

/*
 * Encodes a 32-bit little endian value 'val' and writes to fd.
 *
 * Returns the number of bytes written on success
 *         -1 on failure
 */
int int32_encode(cram_fd *fd, int32_t val) {
    val = le_int4(val);
    if (4 != hwrite(fd->fp, &val, 4))
	return -1;

    return 4;
}

/* As int32_decoded/encode, but from/to blocks instead of cram_fd */
int int32_get(cram_block *b, int32_t *val) {
    if (b->uncomp_size - BLOCK_SIZE(b) < 4)
	return -1;

    *val =
	 b->data[b->byte  ]        |
	(b->data[b->byte+1] <<  8) |
	(b->data[b->byte+2] << 16) |
	(b->data[b->byte+3] << 24);
    BLOCK_SIZE(b) += 4;
    return 4;
}

/* As int32_decoded/encode, but from/to blocks instead of cram_fd */
int int32_put(cram_block *b, int32_t val) {
    unsigned char cp[4];
    cp[0] = ( val      & 0xff);
    cp[1] = ((val>>8)  & 0xff);
    cp[2] = ((val>>16) & 0xff);
    cp[3] = ((val>>24) & 0xff);

    BLOCK_APPEND(b, cp, 4);
    return b->data ? 0 : -1;
}

/* ----------------------------------------------------------------------
 * zlib compression code - from Gap5's tg_iface_g.c
 * They're static here as they're only used within the cram_compress_block
 * and cram_uncompress_block functions, which are the external interface.
 */
char *zlib_mem_inflate(char *cdata, size_t csize, size_t *size) {
    z_stream s;
    unsigned char *data = NULL; /* Uncompressed output */
    int data_alloc = 0;
    int err;

    /* Starting point at uncompressed size, and scale after that */
	data = (unsigned char*)malloc(data_alloc = csize*1.2 + 100);
    if (!data)
	return NULL;

    /* Initialise zlib stream */
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)cdata;
    s.avail_in = csize;
    s.total_in = 0;
    s.next_out  = data;
    s.avail_out = data_alloc;
    s.total_out = 0;

    //err = inflateInit(&s);
    err = inflateInit2(&s, 15 + 32);
    if (err != Z_OK) {
	fprintf(stderr, "zlib inflateInit error: %s\n", s.msg);
	free(data);
	return NULL;
    }

    /* Decode to 'data' array */
    for (;s.avail_in;) {
	unsigned char *data_tmp;
	int alloc_inc;

	s.next_out = &data[s.total_out];
	err = inflate(&s, Z_NO_FLUSH);
	if (err == Z_STREAM_END)
	    break;

	if (err != Z_OK) {
	    fprintf(stderr, "zlib inflate error: %s\n", s.msg);
	    break;
	}

	/* More to come, so realloc based on growth so far */
	alloc_inc = (double)s.avail_in/s.total_in * s.total_out + 100;
	data = (unsigned char*)realloc((data_tmp = data), data_alloc += alloc_inc);
	if (!data) {
	    free(data_tmp);
	    return NULL;
	}
	s.avail_out += alloc_inc;
    }
    inflateEnd(&s);

    *size = s.total_out;
    return (char *)data;
}

static char *zlib_mem_deflate(char *data, size_t size, size_t *cdata_size,
			      int level, int strat) {
    z_stream s;
    unsigned char *cdata = NULL; /* Compressed output */
    int cdata_alloc = 0;
    int cdata_pos = 0;
    int err;

	cdata = (unsigned char*)malloc(cdata_alloc = size*1.05 + 100);
    if (!cdata)
	return NULL;
    cdata_pos = 0;

    /* Initialise zlib stream */
    s.zalloc = Z_NULL; /* use default allocation functions */
    s.zfree  = Z_NULL;
    s.opaque = Z_NULL;
    s.next_in  = (unsigned char *)data;
    s.avail_in = size;
    s.total_in = 0;
    s.next_out  = cdata;
    s.avail_out = cdata_alloc;
    s.total_out = 0;
    s.data_type = Z_BINARY;

    err = deflateInit2(&s, level, Z_DEFLATED, 15|16, 9, strat);
    if (err != Z_OK) {
	fprintf(stderr, "zlib deflateInit2 error: %s\n", s.msg);
	return NULL;
    }

    /* Encode to 'cdata' array */
    for (;s.avail_in;) {
	s.next_out = &cdata[cdata_pos];
	s.avail_out = cdata_alloc - cdata_pos;
	if (cdata_alloc - cdata_pos <= 0) {
	    fprintf(stderr, "Deflate produced larger output than expected. Abort\n"); 
	    return NULL;
	}
	err = deflate(&s, Z_NO_FLUSH);
	cdata_pos = cdata_alloc - s.avail_out;
	if (err != Z_OK) {
	    fprintf(stderr, "zlib deflate error: %s\n", s.msg);
	    break;
	}
    }
    if (deflate(&s, Z_FINISH) != Z_STREAM_END) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    *cdata_size = s.total_out;

    if (deflateEnd(&s) != Z_OK) {
	fprintf(stderr, "zlib deflate error: %s\n", s.msg);
    }
    return (char *)cdata;
}

#ifdef HAVE_LIBLZMA
/* ------------------------------------------------------------------------ */
/*
 * Data compression routines using liblzma (xz)
 *
 * On a test set this shrunk the main db from 136157104 bytes to 114796168, but
 * caused tg_index to grow from 2m43.707s to 15m3.961s. Exporting as bfastq
 * went from 18.3s to 36.3s. So decompression suffers too, but not as bad
 * as compression times.
 *
 * For now we disable this functionality. If it's to be reenabled make sure you
 * improve the mem_inflate implementation as it's just a test hack at the
 * moment.
 */

static char *lzma_mem_deflate(char *data, size_t size, size_t *cdata_size,
			      int level) {
    char *out;
    size_t out_size = lzma_stream_buffer_bound(size);
    *cdata_size = 0;

    out = malloc(out_size);

    /* Single call compression */
    if (LZMA_OK != lzma_easy_buffer_encode(level, LZMA_CHECK_CRC32, NULL,
					   (uint8_t *)data, size,
					   (uint8_t *)out, cdata_size,
					   out_size))
    	return NULL;

    return out;
}

static char *lzma_mem_inflate(char *cdata, size_t csize, size_t *size) {
    lzma_stream strm = LZMA_STREAM_INIT;
    size_t out_size = 0, out_pos = 0;
    char *out = NULL;
    int r;

    /* Initiate the decoder */
    if (LZMA_OK != lzma_stream_decoder(&strm, 50000000, 0))
	return NULL;

    /* Decode loop */
    strm.avail_in = csize;
    strm.next_in = (uint8_t *)cdata;

    for (;strm.avail_in;) {
	if (strm.avail_in > out_size - out_pos) {
	    out_size += strm.avail_in * 4 + 32768;
	    out = realloc(out, out_size);
	}
	strm.avail_out = out_size - out_pos;
	strm.next_out = (uint8_t *)&out[out_pos];

	r = lzma_code(&strm, LZMA_RUN);
	if (LZMA_OK != r && LZMA_STREAM_END != r) {
	    fprintf(stderr, "r=%d\n", r);
	    fprintf(stderr, "mem=%"PRId64"d\n", (int64_t)lzma_memusage(&strm));
	    return NULL;
	}

	out_pos = strm.total_out;

	if (r == LZMA_STREAM_END)
	    break;
    }

    /* finish up any unflushed data; necessary? */
    r = lzma_code(&strm, LZMA_FINISH);
    if (r != LZMA_OK && r != LZMA_STREAM_END) {
	fprintf(stderr, "r=%d\n", r);
	return NULL;
    }

    out = realloc(out, strm.total_out);
    *size = strm.total_out;

    lzma_end(&strm);

    return out;
}
#endif

/* ----------------------------------------------------------------------
 * CRAM blocks - the dynamically growable data block. We have code to
 * create, update, (un)compress and read/write.
 *
 * These are derived from the deflate_interlaced.c blocks, but with the
 * CRAM extension of content types and IDs.
 */

/*
 * Allocates a new cram_block structure with a specified content_type and
 * id.
 *
 * Returns block pointer on success
 *         NULL on failure
 */
cram_block *cram_new_block(enum cram_content_type content_type,
			   int content_id) {
	cram_block *b = (cram_block*)malloc(sizeof(*b));
    if (!b)
	return NULL;
    b->method = b->orig_method = RAW;
    b->content_type = content_type;
    b->content_id = content_id;
    b->comp_size = 0;
    b->uncomp_size = 0;
    b->data = NULL;
    b->alloc = 0;
    b->byte = 0;
    b->bit = 7; // MSB

    return b;
}

/*
 * Reads a block from a cram file.
 * Returns cram_block pointer on success.
 *         NULL on failure
 */
cram_block *cram_read_block(cram_fd *fd) {
	cram_block *b = (cram_block*)malloc(sizeof(*b));
    if (!b)
	return NULL;

    //fprintf(stderr, "Block at %d\n", (int)ftell(fd->fp));

    if (-1 == (b->method      = (cram_block_method)hgetc(fd->fp))) { free(b); return NULL; }
    if (-1 == (b->content_type= (cram_content_type)hgetc(fd->fp))) { free(b); return NULL; }
    if (-1 == itf8_decode(fd, &b->content_id))  { free(b); return NULL; }
    if (-1 == itf8_decode(fd, &b->comp_size))   { free(b); return NULL; }
    if (-1 == itf8_decode(fd, &b->uncomp_size)) { free(b); return NULL; }

    //    fprintf(stderr, "  method %d, ctype %d, cid %d, csize %d, ucsize %d\n",
    //	    b->method, b->content_type, b->content_id, b->comp_size, b->uncomp_size);

    if (b->method == RAW) {
	b->alloc = b->uncomp_size;
	if (!(b->data = (unsigned char*)malloc(b->uncomp_size))){ free(b); return NULL; }
	if (b->uncomp_size != hread(fd->fp, b->data, b->uncomp_size)) {
	    free(b->data);
	    free(b);
	    return NULL;
	}
    } else {
	b->alloc = b->comp_size;
	if (!(b->data = (unsigned char*)malloc(b->comp_size)))  { free(b); return NULL; }
	if (b->comp_size != hread(fd->fp, b->data, b->comp_size)) {
	    free(b->data);
	    free(b);
	    return NULL;
	}
    }

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
	unsigned char dat[100], *cp = dat;;
	uint32_t crc;

	
	if (-1 == int32_decode(fd, (int32_t *)&b->crc32)) {
	    free(b);
	    return NULL;
	}

	*cp++ = b->method;
	*cp++ = b->content_type;
	cp += itf8_put(cp, b->content_id);
	cp += itf8_put(cp, b->comp_size);
	cp += itf8_put(cp, b->uncomp_size);
	crc = crc32(0L, dat, cp-dat);
	crc = crc32(crc, b->data ? b->data : (uc *)"", b->alloc);

	if (crc != b->crc32) {
	    fprintf(stderr, "Block CRC32 failure\n");
	    free(b->data);
	    free(b);
	    return NULL;
	}
    }

    b->orig_method = b->method;
    b->idx = 0;
    b->byte = 0;
    b->bit = 7; // MSB

    return b;
}

/*
 * Writes a CRAM block.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_block(cram_fd *fd, cram_block *b) {
    assert(b->method != RAW || (b->comp_size == b->uncomp_size));

    if (hputc(b->method,       fd->fp)  == EOF) return -1;
    if (hputc(b->content_type, fd->fp)  == EOF) return -1;
    if (itf8_encode(fd, b->content_id)  ==  -1) return -1;
    if (itf8_encode(fd, b->comp_size)   ==  -1) return -1;
    if (itf8_encode(fd, b->uncomp_size) ==  -1) return -1;

    if (b->method == RAW) {
	if (b->uncomp_size != hwrite(fd->fp, b->data, b->uncomp_size))
	    return -1;
    } else {
	if (b->comp_size != hwrite(fd->fp, b->data, b->comp_size))
	    return -1;
    }

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
	unsigned char dat[100], *cp = dat;;
	uint32_t crc;

	*cp++ = b->method;
	*cp++ = b->content_type;
	cp += itf8_put(cp, b->content_id);
	cp += itf8_put(cp, b->comp_size);
	cp += itf8_put(cp, b->uncomp_size);
	crc = crc32(0L, dat, cp-dat);

	if (b->method == RAW) {
	    b->crc32 = crc32(crc, b->data ? b->data : (uc*)"", b->uncomp_size);
	} else {
	    b->crc32 = crc32(crc, b->data ? b->data : (uc*)"", b->comp_size);
	}

	if (-1 == int32_encode(fd, b->crc32))
	    return -1;
    }

    return 0;
}

/*
 * Frees a CRAM block, deallocating internal data too.
 */
void cram_free_block(cram_block *b) {
    if (!b)
	return;
    if (b->data)
	free(b->data);
    free(b);
}

/*
 * Uncompresses a CRAM block, if compressed.
 */
int cram_uncompress_block(cram_block *b) {
    char *uncomp;
    size_t uncomp_size = 0;

    if (b->uncomp_size == 0) {
	// blank block
	b->method = RAW;
	return 0;
    }

    switch (b->method) {
    case RAW:
	return 0;

    case GZIP:
	uncomp = zlib_mem_inflate((char *)b->data, b->comp_size, &uncomp_size);
	if (!uncomp)
	    return -1;
	if ((int)uncomp_size != b->uncomp_size) {
	    free(uncomp);
	    return -1;
	}
	free(b->data);
	b->data = (unsigned char *)uncomp;
	b->alloc = uncomp_size;
	b->method = RAW;
	break;

#ifdef HAVE_LIBBZ2
    case BZIP2: {
	unsigned int usize = b->uncomp_size;
	if (!(uncomp = malloc(usize)))
	    return -1;
	if (BZ_OK != BZ2_bzBuffToBuffDecompress(uncomp, &usize,
						(char *)b->data, b->comp_size,
						0, 0)) {
	    free(uncomp);
	    return -1;
	}
	free(b->data);
	b->data = (unsigned char *)uncomp;
	b->alloc = usize;
	b->method = RAW;
	b->uncomp_size = usize; // Just incase it differs
	break;
    }
#else
    case BZIP2:
	fprintf(stderr, "Bzip2 compression is not compiled into this "
		"version.\nPlease rebuild and try again.\n");
	return -1;
#endif

#ifdef HAVE_LIBLZMA
    case LZMA:
	uncomp = lzma_mem_inflate((char *)b->data, b->comp_size, &uncomp_size);
	if (!uncomp)
	    return -1;
	if ((int)uncomp_size != b->uncomp_size)
	    return -1;
	free(b->data);
	b->data = (unsigned char *)uncomp;
	b->alloc = uncomp_size;
	b->method = RAW;
	break;
#else
    case LZMA:
	fprintf(stderr, "Lzma compression is not compiled into this "
		"version.\nPlease rebuild and try again.\n");
	return -1;
	break;
#endif

    case RANS: {
	unsigned int usize = b->uncomp_size, usize2;
	uncomp = (char *)rans_uncompress(b->data, b->comp_size, &usize2);
	assert(usize == usize2);
	free(b->data);
	b->data = (unsigned char *)uncomp;
	b->alloc = usize2;
	b->method = RAW;
	b->uncomp_size = usize2; // Just incase it differs
	//fprintf(stderr, "Expanded %d to %d\n", b->comp_size, b->uncomp_size);
	break;
    }

    default:
	return -1;
    }

    return 0;
}

static char *cram_compress_by_method(char *in, size_t in_size,
				     size_t *out_size,
				     enum cram_block_method method,
				     int level, int strat) {
    switch (method) {
    case GZIP:
	return zlib_mem_deflate(in, in_size, out_size, level, strat);

    case BZIP2: {
#ifdef HAVE_LIBBZ2
	unsigned int comp_size = in_size*1.01 + 600;
	char *comp = malloc(comp_size);
	if (!comp)
	    return NULL;

	if (BZ_OK != BZ2_bzBuffToBuffCompress(comp, &comp_size,
					      in, in_size,
					      level, 0, 30)) {
	    free(comp);
	    return NULL;
	}
	*out_size = comp_size;
	return comp;
#else
	return NULL;
#endif
    }

    case LZMA:
#ifdef HAVE_LIBLZMA
	return lzma_mem_deflate(in, in_size, out_size, level);
#else
	return NULL;
#endif

    case RANS0: {
	unsigned int out_size_i;
	unsigned char *cp;
	cp = rans_compress((unsigned char *)in, in_size, &out_size_i, 0);
	*out_size = out_size_i;
	return (char *)cp;
    }

    case RANS1: {
	unsigned int out_size_i;
	unsigned char *cp;
	
	cp = rans_compress((unsigned char *)in, in_size, &out_size_i, 1);
	*out_size = out_size_i;
	return (char *)cp;
    }

    case RAW:
	break;

    default:
	return NULL;
    }

    return NULL;
}

#if !defined(_WIN32) && defined(USE_PTHREAD)

/*
* Compresses a block using one of two different zlib strategies. If we only
* want one choice set strat2 to be -1.
*
* The logic here is that sometimes Z_RLE does a better job than Z_FILTERED
* or Z_DEFAULT_STRATEGY on quality data. If so, we'd rather use it as it is
* significantly faster.
*/
int cram_compress_block(cram_fd *fd, cram_block *b, cram_metrics *metrics,
	int method, int level) {

	char *comp = NULL;
	size_t comp_size = 0;
	int strat;

	//fprintf(stderr, "IN: block %d, sz %d\n", b->content_id, b->uncomp_size);

	if (method == RAW || level == 0 || b->uncomp_size == 0) {
		b->method = RAW;
		b->comp_size = b->uncomp_size;
		//fprintf(stderr, "Skip block id %d\n", b->content_id);
		return 0;
	}

	if (metrics) {
		pthread_mutex_lock(&fd->metrics_lock);

		if (metrics->trial > 0 || --metrics->next_trial <= 0) {
			size_t sz_best = INT_MAX;
			size_t sz_gz_rle = 0;
			size_t sz_gz_def = 0;
			size_t sz_rans0 = 0;
			size_t sz_rans1 = 0;
			size_t sz_bzip2 = 0;
			size_t sz_lzma = 0;
			int method_best = 0;
			char *c_best = NULL, *c = NULL;

			if (metrics->revised_method)
				method = metrics->revised_method;
			else
				metrics->revised_method = method;

			if (metrics->next_trial == 0) {
				metrics->next_trial = TRIAL_SPAN;
				metrics->trial = NTRIALS;
				metrics->sz_gz_rle /= 2;
				metrics->sz_gz_def /= 2;
				metrics->sz_rans0 /= 2;
				metrics->sz_rans1 /= 2;
				metrics->sz_bzip2 /= 2;
				metrics->sz_lzma /= 2;
			}

			pthread_mutex_unlock(&fd->metrics_lock);
			if (method & (1 << GZIP_RLE)) {
				c = cram_compress_by_method((char *)b->data, b->uncomp_size,
					&sz_gz_rle, GZIP, 1, Z_RLE);
				if (c && sz_best > sz_gz_rle) {
					sz_best = sz_gz_rle;
					method_best = GZIP_RLE;
					if (c_best)
						free(c_best);
					c_best = c;
				}
				else if (c) {
					free(c);
				}
				else {
					sz_gz_rle = b->uncomp_size * 2 + 1000;
				}

				//fprintf(stderr, "Block %d; %d->%d\n", b->content_id, b->uncomp_size, sz_gz_rle);
			}

			if (method & (1 << GZIP)) {
				c = cram_compress_by_method((char *)b->data, b->uncomp_size,
					&sz_gz_def, GZIP, level,
					Z_FILTERED);
				if (c && sz_best > sz_gz_def) {
					sz_best = sz_gz_def;
					method_best = GZIP;
					if (c_best)
						free(c_best);
					c_best = c;
				}
				else if (c) {
					free(c);
				}
				else {
					sz_gz_def = b->uncomp_size * 2 + 1000;
				}

				//fprintf(stderr, "Block %d; %d->%d\n", b->content_id, b->uncomp_size, sz_gz_def);
			}

			if (method & (1 << RANS0)) {
				c = cram_compress_by_method((char *)b->data, b->uncomp_size,
					&sz_rans0, RANS0, 0, 0);
				if (c && sz_best > sz_rans0) {
					sz_best = sz_rans0;
					method_best = RANS0;
					if (c_best)
						free(c_best);
					c_best = c;
				}
				else if (c) {
					free(c);
				}
				else {
					sz_rans0 = b->uncomp_size * 2 + 1000;
				}
			}

			if (method & (1 << RANS1)) {
				c = cram_compress_by_method((char *)b->data, b->uncomp_size,
					&sz_rans1, RANS1, 0, 0);
				if (c && sz_best > sz_rans1) {
					sz_best = sz_rans1;
					method_best = RANS1;
					if (c_best)
						free(c_best);
					c_best = c;
				}
				else if (c) {
					free(c);
				}
				else {
					sz_rans1 = b->uncomp_size * 2 + 1000;
				}
			}

			if (method & (1 << BZIP2)) {
				c = cram_compress_by_method((char *)b->data, b->uncomp_size,
					&sz_bzip2, BZIP2, level, 0);
				if (c && sz_best > sz_bzip2) {
					sz_best = sz_bzip2;
					method_best = BZIP2;
					if (c_best)
						free(c_best);
					c_best = c;
				}
				else if (c) {
					free(c);
				}
				else {
					sz_bzip2 = b->uncomp_size * 2 + 1000;
				}
			}

			if (method & (1 << LZMA)) {
				c = cram_compress_by_method((char *)b->data, b->uncomp_size,
					&sz_lzma, LZMA, level, 0);
				if (c && sz_best > sz_lzma) {
					sz_best = sz_lzma;
					method_best = LZMA;
					if (c_best)
						free(c_best);
					c_best = c;
				}
				else if (c) {
					free(c);
				}
				else {
					sz_lzma = b->uncomp_size * 2 + 1000;
				}
			}

			//fprintf(stderr, "sz_best = %d\n", sz_best);

			free(b->data);
			b->data = (unsigned char *)c_best;
			//printf("method_best = %s\n", cram_block_method2str(method_best));
			b->method = method_best == GZIP_RLE ? GZIP : (cram_block_method)method_best;
			b->comp_size = sz_best;

			pthread_mutex_lock(&fd->metrics_lock);

			metrics->sz_gz_rle += sz_gz_rle;
			metrics->sz_gz_def += sz_gz_def;
			metrics->sz_rans0 += sz_rans0;
			metrics->sz_rans1 += sz_rans1;
			metrics->sz_bzip2 += sz_bzip2;
			metrics->sz_lzma += sz_lzma;
			if (--metrics->trial == 0) {
				int best_method = RAW;
				int best_sz = INT_MAX;

				// Scale methods by cost
				if (fd->level <= 3) {
					metrics->sz_rans1 *= 1.02;
					metrics->sz_gz_def *= 1.04;
					metrics->sz_bzip2 *= 1.08;
					metrics->sz_lzma *= 1.10;
				}
				else if (fd->level <= 6) {
					metrics->sz_rans1 *= 1.01;
					metrics->sz_gz_def *= 1.02;
					metrics->sz_bzip2 *= 1.03;
					metrics->sz_lzma *= 1.05;
				}

				if (method & (1 << GZIP_RLE) && best_sz > metrics->sz_gz_rle)
					best_sz = metrics->sz_gz_rle, best_method = GZIP_RLE;

				if (method & (1 << GZIP) && best_sz > metrics->sz_gz_def)
					best_sz = metrics->sz_gz_def, best_method = GZIP;

				if (method & (1 << RANS0) && best_sz > metrics->sz_rans0)
					best_sz = metrics->sz_rans0, best_method = RANS0;

				if (method & (1 << RANS1) && best_sz > metrics->sz_rans1)
					best_sz = metrics->sz_rans1, best_method = RANS1;

				if (method & (1 << BZIP2) && best_sz > metrics->sz_bzip2)
					best_sz = metrics->sz_bzip2, best_method = BZIP2;

				if (method & (1 << LZMA) && best_sz > metrics->sz_lzma)
					best_sz = metrics->sz_lzma, best_method = LZMA;

				if (best_method == GZIP_RLE) {
					metrics->method = GZIP;
					metrics->strat = Z_RLE;
				}
				else {
					metrics->method = best_method;
					metrics->strat = Z_FILTERED;
				}

				// If we see at least MAXFAIL trials in a row for a specific
				// compression method with more than MAXDELTA aggregate
				// size then we drop this from the list of methods used
				// for this block type.
#define MAXDELTA 0.20
#define MAXFAILS 4
				if (best_method == GZIP_RLE) {
					metrics->gz_rle_cnt = 0;
					metrics->gz_rle_extra = 0;
				}
				else if (best_sz < metrics->sz_gz_rle) {
					double r = (double)metrics->sz_gz_rle / best_sz - 1;
					if (++metrics->gz_rle_cnt >= MAXFAILS &&
						(metrics->gz_rle_extra += r) >= MAXDELTA)
						method &= ~(1 << GZIP_RLE);
				}

				if (best_method == GZIP) {
					metrics->gz_def_cnt = 0;
					metrics->gz_def_extra = 0;
				}
				else if (best_sz < metrics->sz_gz_def) {
					double r = (double)metrics->sz_gz_def / best_sz - 1;
					if (++metrics->gz_def_cnt >= MAXFAILS &&
						(metrics->gz_def_extra += r) >= MAXDELTA)
						method &= ~(1 << GZIP);
				}

				if (best_method == RANS0) {
					metrics->rans0_cnt = 0;
					metrics->rans0_extra = 0;
				}
				else if (best_sz < metrics->sz_rans0) {
					double r = (double)metrics->sz_rans0 / best_sz - 1;
					if (++metrics->rans0_cnt >= MAXFAILS &&
						(metrics->rans0_extra += r) >= MAXDELTA)
						method &= ~(1 << RANS0);
				}

				if (best_method == RANS1) {
					metrics->rans1_cnt = 0;
					metrics->rans1_extra = 0;
				}
				else if (best_sz < metrics->sz_rans1) {
					double r = (double)metrics->sz_rans1 / best_sz - 1;
					if (++metrics->rans1_cnt >= MAXFAILS &&
						(metrics->rans1_extra += r) >= MAXDELTA)
						method &= ~(1 << RANS1);
				}

				if (best_method == BZIP2) {
					metrics->bzip2_cnt = 0;
					metrics->bzip2_extra = 0;
				}
				else if (best_sz < metrics->sz_bzip2) {
					double r = (double)metrics->sz_bzip2 / best_sz - 1;
					if (++metrics->bzip2_cnt >= MAXFAILS &&
						(metrics->bzip2_extra += r) >= MAXDELTA)
						method &= ~(1 << BZIP2);
				}

				if (best_method == LZMA) {
					metrics->lzma_cnt = 0;
					metrics->lzma_extra = 0;
				}
				else if (best_sz < metrics->sz_lzma) {
					double r = (double)metrics->sz_lzma / best_sz - 1;
					if (++metrics->lzma_cnt >= MAXFAILS &&
						(metrics->lzma_extra += r) >= MAXDELTA)
						method &= ~(1 << LZMA);
				}

				//if (method != metrics->revised_method)
				//    fprintf(stderr, "%d: method from %x to %x\n",
				//	    b->content_id, metrics->revised_method, method);
				metrics->revised_method = method;
			}
			pthread_mutex_unlock(&fd->metrics_lock);
		}
		else {
			strat = metrics->strat;
			method = metrics->method;

			pthread_mutex_unlock(&fd->metrics_lock);
			comp = cram_compress_by_method((char *)b->data, b->uncomp_size,
				&comp_size, (cram_block_method)method,
				level, strat);
			if (!comp)
				return -1;
			free(b->data);
			b->data = (unsigned char *)comp;
			b->comp_size = comp_size;
			b->method = (cram_block_method)method;
		}

	}
	else {
		// no cached metrics, so just do zlib?
		comp = cram_compress_by_method((char *)b->data, b->uncomp_size,
			&comp_size, GZIP, level, Z_FILTERED);
		if (!comp) {
			fprintf(stderr, "Compression failed!\n");
			return -1;
		}
		free(b->data);
		b->data = (unsigned char *)comp;
		b->comp_size = comp_size;
		b->method = GZIP;
	}

	if (fd->verbose)
		fprintf(stderr, "Compressed block ID %d from %d to %d by method %s\n",
		b->content_id, b->uncomp_size, b->comp_size,
		cram_block_method2str(b->method));

	if (b->method == RANS1)
		b->method = RANS0; // Spec just has RANS (not 0/1) with auto-sensing

	return 0;
}

#else // ~ !defined(_WIN32) && defined(USE_PTHREAD)

/*
* Compresses a block using one of two different zlib strategies. If we only
* want one choice set strat2 to be -1.
*
* The logic here is that sometimes Z_RLE does a better job than Z_FILTERED
* or Z_DEFAULT_STRATEGY on quality data. If so, we'd rather use it as it is
* significantly faster.
*/
int cram_compress_block(cram_fd *fd, cram_block *b, cram_metrics *metrics,
	int method, int level) {

	char *comp = NULL;
	size_t comp_size = 0;
	int strat;

	//fprintf(stderr, "IN: block %d, sz %d\n", b->content_id, b->uncomp_size);

	if (method == RAW || level == 0 || b->uncomp_size == 0) {
		b->method = RAW;
		b->comp_size = b->uncomp_size;
		//fprintf(stderr, "Skip block id %d\n", b->content_id);
		return 0;
	}

	if (metrics) {
		fd->metrics_lock.lock();

		if (metrics->trial > 0 || --metrics->next_trial <= 0) {
			size_t sz_best = INT_MAX;
			size_t sz_gz_rle = 0;
			size_t sz_gz_def = 0;
			size_t sz_rans0 = 0;
			size_t sz_rans1 = 0;
			size_t sz_bzip2 = 0;
			size_t sz_lzma = 0;
			int method_best = 0;
			char *c_best = NULL, *c = NULL;

			if (metrics->revised_method)
				method = metrics->revised_method;
			else
				metrics->revised_method = method;

			if (metrics->next_trial == 0) {
				metrics->next_trial = TRIAL_SPAN;
				metrics->trial = NTRIALS;
				metrics->sz_gz_rle /= 2;
				metrics->sz_gz_def /= 2;
				metrics->sz_rans0 /= 2;
				metrics->sz_rans1 /= 2;
				metrics->sz_bzip2 /= 2;
				metrics->sz_lzma /= 2;
			}

			fd->metrics_lock.unlock();
			if (method & (1 << GZIP_RLE)) {
				c = cram_compress_by_method((char *)b->data, b->uncomp_size,
					&sz_gz_rle, GZIP, 1, Z_RLE);
				if (c && sz_best > sz_gz_rle) {
					sz_best = sz_gz_rle;
					method_best = GZIP_RLE;
					if (c_best)
						free(c_best);
					c_best = c;
				}
				else if (c) {
					free(c);
				}
				else {
					sz_gz_rle = b->uncomp_size * 2 + 1000;
				}

				//fprintf(stderr, "Block %d; %d->%d\n", b->content_id, b->uncomp_size, sz_gz_rle);
			}

			if (method & (1 << GZIP)) {
				c = cram_compress_by_method((char *)b->data, b->uncomp_size,
					&sz_gz_def, GZIP, level,
					Z_FILTERED);
				if (c && sz_best > sz_gz_def) {
					sz_best = sz_gz_def;
					method_best = GZIP;
					if (c_best)
						free(c_best);
					c_best = c;
				}
				else if (c) {
					free(c);
				}
				else {
					sz_gz_def = b->uncomp_size * 2 + 1000;
				}

				//fprintf(stderr, "Block %d; %d->%d\n", b->content_id, b->uncomp_size, sz_gz_def);
			}

			if (method & (1 << RANS0)) {
				c = cram_compress_by_method((char *)b->data, b->uncomp_size,
					&sz_rans0, RANS0, 0, 0);
				if (c && sz_best > sz_rans0) {
					sz_best = sz_rans0;
					method_best = RANS0;
					if (c_best)
						free(c_best);
					c_best = c;
				}
				else if (c) {
					free(c);
				}
				else {
					sz_rans0 = b->uncomp_size * 2 + 1000;
				}
			}

			if (method & (1 << RANS1)) {
				c = cram_compress_by_method((char *)b->data, b->uncomp_size,
					&sz_rans1, RANS1, 0, 0);
				if (c && sz_best > sz_rans1) {
					sz_best = sz_rans1;
					method_best = RANS1;
					if (c_best)
						free(c_best);
					c_best = c;
				}
				else if (c) {
					free(c);
				}
				else {
					sz_rans1 = b->uncomp_size * 2 + 1000;
				}
			}

			if (method & (1 << BZIP2)) {
				c = cram_compress_by_method((char *)b->data, b->uncomp_size,
					&sz_bzip2, BZIP2, level, 0);
				if (c && sz_best > sz_bzip2) {
					sz_best = sz_bzip2;
					method_best = BZIP2;
					if (c_best)
						free(c_best);
					c_best = c;
				}
				else if (c) {
					free(c);
				}
				else {
					sz_bzip2 = b->uncomp_size * 2 + 1000;
				}
			}

			if (method & (1 << LZMA)) {
				c = cram_compress_by_method((char *)b->data, b->uncomp_size,
					&sz_lzma, LZMA, level, 0);
				if (c && sz_best > sz_lzma) {
					sz_best = sz_lzma;
					method_best = LZMA;
					if (c_best)
						free(c_best);
					c_best = c;
				}
				else if (c) {
					free(c);
				}
				else {
					sz_lzma = b->uncomp_size * 2 + 1000;
				}
			}

			//fprintf(stderr, "sz_best = %d\n", sz_best);

			free(b->data);
			b->data = (unsigned char *)c_best;
			//printf("method_best = %s\n", cram_block_method2str(method_best));
			b->method = method_best == GZIP_RLE ? GZIP : (cram_block_method)method_best;
			b->comp_size = sz_best;

			fd->metrics_lock.lock();

			metrics->sz_gz_rle += sz_gz_rle;
			metrics->sz_gz_def += sz_gz_def;
			metrics->sz_rans0 += sz_rans0;
			metrics->sz_rans1 += sz_rans1;
			metrics->sz_bzip2 += sz_bzip2;
			metrics->sz_lzma += sz_lzma;
			if (--metrics->trial == 0) {
				int best_method = RAW;
				int best_sz = INT_MAX;

				// Scale methods by cost
				if (fd->level <= 3) {
					metrics->sz_rans1 *= 1.02;
					metrics->sz_gz_def *= 1.04;
					metrics->sz_bzip2 *= 1.08;
					metrics->sz_lzma *= 1.10;
				}
				else if (fd->level <= 6) {
					metrics->sz_rans1 *= 1.01;
					metrics->sz_gz_def *= 1.02;
					metrics->sz_bzip2 *= 1.03;
					metrics->sz_lzma *= 1.05;
				}

				if (method & (1 << GZIP_RLE) && best_sz > metrics->sz_gz_rle)
					best_sz = metrics->sz_gz_rle, best_method = GZIP_RLE;

				if (method & (1 << GZIP) && best_sz > metrics->sz_gz_def)
					best_sz = metrics->sz_gz_def, best_method = GZIP;

				if (method & (1 << RANS0) && best_sz > metrics->sz_rans0)
					best_sz = metrics->sz_rans0, best_method = RANS0;

				if (method & (1 << RANS1) && best_sz > metrics->sz_rans1)
					best_sz = metrics->sz_rans1, best_method = RANS1;

				if (method & (1 << BZIP2) && best_sz > metrics->sz_bzip2)
					best_sz = metrics->sz_bzip2, best_method = BZIP2;

				if (method & (1 << LZMA) && best_sz > metrics->sz_lzma)
					best_sz = metrics->sz_lzma, best_method = LZMA;

				if (best_method == GZIP_RLE) {
					metrics->method = GZIP;
					metrics->strat = Z_RLE;
				}
				else {
					metrics->method = best_method;
					metrics->strat = Z_FILTERED;
				}

				// If we see at least MAXFAIL trials in a row for a specific
				// compression method with more than MAXDELTA aggregate
				// size then we drop this from the list of methods used
				// for this block type.
#define MAXDELTA 0.20
#define MAXFAILS 4
				if (best_method == GZIP_RLE) {
					metrics->gz_rle_cnt = 0;
					metrics->gz_rle_extra = 0;
				}
				else if (best_sz < metrics->sz_gz_rle) {
					double r = (double)metrics->sz_gz_rle / best_sz - 1;
					if (++metrics->gz_rle_cnt >= MAXFAILS &&
						(metrics->gz_rle_extra += r) >= MAXDELTA)
						method &= ~(1 << GZIP_RLE);
				}

				if (best_method == GZIP) {
					metrics->gz_def_cnt = 0;
					metrics->gz_def_extra = 0;
				}
				else if (best_sz < metrics->sz_gz_def) {
					double r = (double)metrics->sz_gz_def / best_sz - 1;
					if (++metrics->gz_def_cnt >= MAXFAILS &&
						(metrics->gz_def_extra += r) >= MAXDELTA)
						method &= ~(1 << GZIP);
				}

				if (best_method == RANS0) {
					metrics->rans0_cnt = 0;
					metrics->rans0_extra = 0;
				}
				else if (best_sz < metrics->sz_rans0) {
					double r = (double)metrics->sz_rans0 / best_sz - 1;
					if (++metrics->rans0_cnt >= MAXFAILS &&
						(metrics->rans0_extra += r) >= MAXDELTA)
						method &= ~(1 << RANS0);
				}

				if (best_method == RANS1) {
					metrics->rans1_cnt = 0;
					metrics->rans1_extra = 0;
				}
				else if (best_sz < metrics->sz_rans1) {
					double r = (double)metrics->sz_rans1 / best_sz - 1;
					if (++metrics->rans1_cnt >= MAXFAILS &&
						(metrics->rans1_extra += r) >= MAXDELTA)
						method &= ~(1 << RANS1);
				}

				if (best_method == BZIP2) {
					metrics->bzip2_cnt = 0;
					metrics->bzip2_extra = 0;
				}
				else if (best_sz < metrics->sz_bzip2) {
					double r = (double)metrics->sz_bzip2 / best_sz - 1;
					if (++metrics->bzip2_cnt >= MAXFAILS &&
						(metrics->bzip2_extra += r) >= MAXDELTA)
						method &= ~(1 << BZIP2);
				}

				if (best_method == LZMA) {
					metrics->lzma_cnt = 0;
					metrics->lzma_extra = 0;
				}
				else if (best_sz < metrics->sz_lzma) {
					double r = (double)metrics->sz_lzma / best_sz - 1;
					if (++metrics->lzma_cnt >= MAXFAILS &&
						(metrics->lzma_extra += r) >= MAXDELTA)
						method &= ~(1 << LZMA);
				}

				//if (method != metrics->revised_method)
				//    fprintf(stderr, "%d: method from %x to %x\n",
				//	    b->content_id, metrics->revised_method, method);
				metrics->revised_method = method;
			}
			fd->metrics_lock.unlock();
		}
		else {
			strat = metrics->strat;
			method = metrics->method;

			fd->metrics_lock.unlock();
			comp = cram_compress_by_method((char *)b->data, b->uncomp_size,
				&comp_size, (cram_block_method)method,
				level, strat);
			if (!comp)
				return -1;
			free(b->data);
			b->data = (unsigned char *)comp;
			b->comp_size = comp_size;
			b->method = (cram_block_method)method;
		}

	}
	else {
		// no cached metrics, so just do zlib?
		comp = cram_compress_by_method((char *)b->data, b->uncomp_size,
			&comp_size, GZIP, level, Z_FILTERED);
		if (!comp) {
			fprintf(stderr, "Compression failed!\n");
			return -1;
		}
		free(b->data);
		b->data = (unsigned char *)comp;
		b->comp_size = comp_size;
		b->method = GZIP;
	}

	if (fd->verbose)
		fprintf(stderr, "Compressed block ID %d from %d to %d by method %s\n",
		b->content_id, b->uncomp_size, b->comp_size,
		cram_block_method2str(b->method));

	if (b->method == RANS1)
		b->method = RANS0; // Spec just has RANS (not 0/1) with auto-sensing

	return 0;
}

#endif // ~ !defined(_WIN32) && defined(USE_PTHREAD)

cram_metrics *cram_new_metrics(void) {
	cram_metrics *m = (cram_metrics*)calloc(1, sizeof(*m));
	if (!m)
		return NULL;
	m->trial = NTRIALS - 1;
	m->next_trial = TRIAL_SPAN;
	m->method = RAW;
	m->strat = 0;
	m->revised_method = 0;

	return m;
}

char *cram_block_method2str(enum cram_block_method m) {
    switch(m) {
    case RAW:	   return "RAW";
    case GZIP:	   return "GZIP";
    case BZIP2:	   return "BZIP2";
    case LZMA:     return "LZMA";
    case RANS0:    return "RANS0";
    case RANS1:    return "RANS1";
    case GZIP_RLE: return "GZIP_RLE";
    case ERROR:    break;
    }
    return "?";
}

char *cram_content_type2str(enum cram_content_type t) {
    switch (t) {
    case FILE_HEADER:         return "FILE_HEADER";
    case COMPRESSION_HEADER:  return "COMPRESSION_HEADER";
    case MAPPED_SLICE:        return "MAPPED_SLICE";
    case UNMAPPED_SLICE:      return "UNMAPPED_SLICE";
    case EXTERNAL:            return "EXTERNAL";
    case CORE:                return "CORE";
    case CT_ERROR:            break;
    }
    return "?";
}

/*
 * Extra error checking on fclose to really ensure data is written.
 * Care needs to be taken to handle pipes vs real files.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
int paranoid_fclose(FILE *fp) {
	if (-1 == fflush(fp) && errno != EBADF) {
		fclose(fp);
		return -1;
	}
	errno = 0;
	int ret;
	// TODO : revisit this code, is this a perfect replacement for fsync ? 
#ifdef _WIN32
	ret = fflush(fp);
	if (ret != 0)
	{
		fclose(fp);
		return -1;
	}
#else
	if (-1 == fsync(fileno(fp))) {
		if (errno != EINVAL) { // eg pipe
			fclose(fp);
			return -1;
		}
}
#endif
    return fclose(fp);
}

/* ----------------------------------------------------------------------
 * Reference sequence handling
 *
 * These revolve around the refs_t structure, which may potentially be
 * shared between multiple cram_fd.
 *
 * We start with refs_create() to allocate an empty refs_t and then
 * populate it with @SQ line data using refs_from_header(). This is done on
 * cram_open().  Also at start up we can call cram_load_reference() which
 * is used with "scramble -r foo.fa". This replaces the fd->refs with the
 * new one specified. In either case refs2id() is then called which
 * maps ref_entry names to @SQ ids (refs_t->ref_id[]).
 *
 * Later, possibly within a thread, we will want to know the actual ref
 * seq itself, obtained by calling cram_get_ref().  This may use the
 * UR: or M5: fields or the filename specified in the original
 * cram_load_reference() call.
 *
 * Given the potential for multi-threaded reference usage, we have
 * reference counting (sorry for the confusing double use of "ref") to
 * track the number of callers interested in any specific reference.
 */

void refs_free(refs_t *r) {
    RP("refs_free()\n");

    if (--r->count > 0)
	return;

    if (!r)
	return;

    if (r->pool)
	string_pool_destroy(r->pool);

    if (r->h_meta) {
	khint_t k;

	for (k = kh_begin(r->h_meta); k != kh_end(r->h_meta); k++) {
	    ref_entry *e;

	    if (!kh_exist(r->h_meta, k))
		continue;
	    if (!(e = kh_val(r->h_meta, k)))
		continue;
	    if (e->seq)
		free(e->seq);
	    free(e);
	}

	kh_destroy(refs, r->h_meta);
    }
    
    if (r->ref_id)
	free(r->ref_id);

    if (r->fp)
	bgzf_close(r->fp);

#ifndef _WIN32
    pthread_mutex_destroy(&r->lock);
#endif

    free(r);
}

static refs_t *refs_create(void) {
	refs_t *r = (refs_t*)calloc(1, sizeof(*r));

    RP("refs_create()\n");

    if (!r)
	return NULL;

    if (!(r->pool = string_pool_create(8192)))
	goto err;

    r->ref_id = NULL; // see refs2id() to populate.
    r->count = 1;
    r->last = NULL;
    r->last_id = -1;

    if (!(r->h_meta = kh_init(refs)))
	goto err;

#ifndef _WIN32
    pthread_mutex_init(&r->lock, NULL);
#endif
    return r;

 err:
    refs_free(r);
    return NULL;
}

/*
 * Opens a reference fasta file as a BGZF stream, allowing for
 * compressed files.  It automatically builds a .fai file if
 * required and if compressed a .gzi bgzf index too.
 *
 * Returns a BGZF handle on success;
 *         NULL on failure.
 */
static BGZF *bgzf_open_ref(char *fn, char *mode) {
    BGZF *fp;
    char fai_file[PATH_MAX];
#ifdef _WIN32
	_snprintf_s(fai_file, PATH_MAX, "%s.fai", fn);
#else
    snprintf(fai_file, PATH_MAX, "%s.fai", fn);
#endif
    if (access(fai_file, R_OK) != 0)
	if (fai_build(fn) != 0)
	    return NULL;

    if (!(fp = bgzf_open(fn, mode))) {
	perror(fn);
	return NULL;
    }

    if (fp->is_compressed == 1 && bgzf_index_load(fp, fn, ".gzi") < 0) {
	fprintf(stderr, "Unable to load .gzi index '%s.gzi'\n", fn);
	bgzf_close(fp);
	return NULL;
    }

    return fp;
}

/*
 * Loads a FAI file for a reference.fasta.
 * "is_err" indicates whether failure to load is worthy of emitting an
 * error message. In some cases (eg with embedded references) we
 * speculatively load, just incase, and silently ignore errors.
 *
 * Returns the refs_t struct on success (maybe newly allocated);
 *         NULL on failure
 */
static refs_t *refs_load_fai(refs_t *r_orig, char *fn, int is_err) {
    struct stat sb;
    FILE *fp = NULL;
    char fai_fn[PATH_MAX];
    char line[8192];
    refs_t *r = r_orig;
    size_t fn_l = strlen(fn);
    int id = 0, id_alloc = 0;

    RP("refs_load_fai %s\n", fn);

    if (!r)
	if (!(r = refs_create()))
	    goto err;

    /* Open reference, for later use */
    if (stat(fn, &sb) != 0) {
	if (is_err)
	    perror(fn);
	goto err;
    }

    if (r->fp)
	if (bgzf_close(r->fp) != 0)
	    goto err;
    r->fp = NULL;

    if (!(r->fn = string_dup(r->pool, fn)))
	goto err;
	
    if (fn_l > 4 && strcmp(&fn[fn_l-4], ".fai") == 0)
	r->fn[fn_l-4] = 0;

    if (!(r->fp = bgzf_open_ref(r->fn, "r")))
	goto err;

    /* Parse .fai file and load meta-data */
    sprintf(fai_fn, "%.*s.fai", PATH_MAX-5, r->fn);

    if (stat(fai_fn, &sb) != 0) {
	if (is_err)
	    perror(fai_fn);
	goto err;
    }
    if (!(fp = fopen(fai_fn, "r"))) {
	if (is_err)
	    perror(fai_fn);
	goto err;
    }
    while (fgets(line, 8192, fp) != NULL) {
		ref_entry *e = (ref_entry*)malloc(sizeof(*e));
	char *cp;
	int n;
	khint_t k;

	if (!e)
	    return NULL;

	// id
	for (cp = line; *cp && !isspace(*cp); cp++)
	    ;
	*cp++ = 0;
	e->name = string_dup(r->pool, line);
	
	// length
	while (*cp && isspace(*cp))
	    cp++;
	e->length = strtoll(cp, &cp, 10);

	// offset
	while (*cp && isspace(*cp))
	    cp++;
	e->offset = strtoll(cp, &cp, 10);

	// bases per line
	while (*cp && isspace(*cp))
	    cp++;
	e->bases_per_line = strtol(cp, &cp, 10);

	// line length
	while (*cp && isspace(*cp))
	    cp++;
	e->line_length = strtol(cp, &cp, 10);

	// filename
	e->fn = r->fn;

	e->count = 0;
	e->seq = NULL;

	k = kh_put(refs, r->h_meta, e->name, &n);
	if (-1 == n)  {
	    free(e);
	    return NULL;
	}

	if (n) {
	    kh_val(r->h_meta, k) = e;
	} else {
	    ref_entry *re = kh_val(r->h_meta, k);
	    if (re && (re->count != 0 || re->length != 0)) {
		/* Keep old */
		free(e);
	    } else {
		/* Replace old */
		if (re)
		    free(re);
		kh_val(r->h_meta, k) = e;
	    }
	}

	if (id >= id_alloc) {
	    int x;

	    id_alloc = id_alloc ?id_alloc*2 : 16;
	    r->ref_id = (ref_entry**)realloc(r->ref_id, id_alloc * sizeof(*r->ref_id));

	    for (x = id; x < id_alloc; x++)
		r->ref_id[x] = NULL;
	}
	r->ref_id[id] = e;
	r->nref = ++id;
    }

    return r;

 err:
    if (fp)
	fclose(fp);

    if (!r_orig)
	refs_free(r);
    
    return NULL;
}

/*
 * Indexes references by the order they appear in a BAM file. This may not
 * necessarily be the same order they appear in the fasta reference file.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int refs2id(refs_t *r, SAM_hdr *h) {
    int i;

    if (r->ref_id)
	free(r->ref_id);
    if (r->last)
	r->last = NULL;

    r->ref_id = (ref_entry**)calloc(h->nref, sizeof(*r->ref_id));
    if (!r->ref_id)
	return -1;

    r->nref = h->nref;
    for (i = 0; i < h->nref; i++) {
	khint_t k = kh_get(refs, r->h_meta, h->ref[i].name);
	if (k != kh_end(r->h_meta)) {
	    r->ref_id[i] = kh_val(r->h_meta, k);
	} else {
	    fprintf(stderr, "Unable to find ref name '%s'\n",
		    h->ref[i].name);
	}
    }

    return 0;
}

/*
 * Generates refs_t entries based on @SQ lines in the header.
 * Returns 0 on success
 *         -1 on failure
 */
static int refs_from_header(refs_t *r, cram_fd *fd, SAM_hdr *h) {
    int i, j;

    if (!h || h->nref == 0)
	return 0;

    //fprintf(stderr, "refs_from_header for %p mode %c\n", fd, fd->mode);

    /* Existing refs are fine, as long as they're compatible with the hdr. */
    if (!(r->ref_id = (ref_entry**)realloc(r->ref_id, (r->nref + h->nref) * sizeof(*r->ref_id))))
	return -1;

    /* Copy info from h->ref[i] over to r */
    for (i = 0, j = r->nref; i < h->nref; i++) {
	SAM_hdr_type *ty;
	SAM_hdr_tag *tag;
	khint_t k;
	int n;

	k = kh_get(refs, r->h_meta, h->ref[i].name);
	if (k != kh_end(r->h_meta))
	    // Ref already known about
	    continue;

	if (!(r->ref_id[j] = (ref_entry*)calloc(1, sizeof(ref_entry))))
	    return -1;

	if (!h->ref[j].name)
	    return -1;

	r->ref_id[j]->name = string_dup(r->pool, h->ref[i].name);
	r->ref_id[j]->length = 0; // marker for not yet loaded

	/* Initialise likely filename if known */
	if ((ty = sam_hdr_find(h, "SQ", "SN", h->ref[i].name))) {
	    if ((tag = sam_hdr_find_key(h, ty, "M5", NULL))) {
		r->ref_id[j]->fn = string_dup(r->pool, tag->str+3);
		//fprintf(stderr, "Tagging @SQ %s / %s\n", r->ref_id[h]->name, r->ref_id[h]->fn);
	    }
	}

	k = kh_put(refs, r->h_meta, r->ref_id[j]->name, &n);
	if (n <= 0) // already exists or error
	    return -1;
	kh_val(r->h_meta, k) = r->ref_id[j];

	j++;
    }
    r->nref = j;

    return 0;
}

/*
 * Attaches a header to a cram_fd.
 *
 * This should be used when creating a new cram_fd for writing where
 * we have an SAM_hdr already constructed (eg from a file we've read
 * in).
 */
int cram_set_header(cram_fd *fd, SAM_hdr *hdr) {
    if (fd->header)
	sam_hdr_free(fd->header);
    fd->header = hdr;
    return refs_from_header(fd->refs, fd, hdr);
}

/*
 * Converts a directory and a filename into an expanded path, replacing %s
 * in directory with the filename and %[0-9]+s with portions of the filename
 * Any remaining parts of filename are added to the end with /%s.
 */
void expand_cache_path(char *path, char *dir, char *fn) {
    char *cp;

    while ((cp = strchr(dir, '%'))) {
	strncpy(path, dir, cp-dir);
	path += cp-dir;

	if (*++cp == 's') {
	    strcpy(path, fn);
	    path += strlen(fn);
	    fn += strlen(fn);
	    cp++;
	} else if (*cp >= '0' && *cp <= '9') {
	    char *endp;
	    long l;

	    l = strtol(cp, &endp, 10);
	    l = MIN(l, strlen(fn));
	    if (*endp == 's') {
		strncpy(path, fn, l);
		path += l;
		fn += l;
		*path = 0;
		cp = endp+1;
	    } else {
		*path++ = '%';
		*path++ = *cp++;
	    }
	} else {
	    *path++ = '%';
	    *path++ = *cp++;
	}
	dir = cp;
    }
    strcpy(path, dir);
    path += strlen(dir);
    if (*fn && path[-1] != '/')
	*path++ = '/';
    strcpy(path, fn);
}

/*
 * Make the directory containing path and any prefix directories.
 */
void mkdir_prefix(char *path, int mode) {
    char *cp = strrchr(path, '/');
    if (!cp)
	return;

    *cp = 0;
    if (is_directory(path)) {
	*cp = '/';
	return;
    }

#ifdef _WIN32
    if (_mkdir(path) == 0) 
	{
		*cp = '/';
		return;
    }
#else
	if (mkdir(path, mode) == 0)
	{
		chmod(path, mode);
		*cp = '/';
		return;
	}
#endif

    mkdir_prefix(path, mode);
#ifdef _WIN32
	_mkdir(path); 
#else
    mkdir(path, mode);
    chmod(path, mode);
#endif
    *cp = '/';
}

/*
 * Return the cache directory to use, based on the first of these
 * environment variables to be set to a non-empty value.
 */
static const char *get_cache_basedir(const char **extra) {
    char *base;

    *extra = "";

    base = getenv("XDG_CACHE_HOME");
    if (base && *base) return base;

    base = getenv("HOME");
    if (base && *base) { *extra = "/.cache"; return base; }

    base = getenv("TMPDIR");
    if (base && *base) return base;

    base = getenv("TEMP");
    if (base && *base) return base;

    return "/tmp";
}

/*
 * Queries the M5 string from the header and attempts to populate the
 * reference from this using the REF_PATH environment.
 *
 * Returns 0 on sucess
 *        -1 on failure
 */
static int cram_populate_ref(cram_fd *fd, int id, ref_entry *r) {
    char *ref_path = getenv("REF_PATH");
    SAM_hdr_type *ty;
    SAM_hdr_tag *tag;
    char path[PATH_MAX], path_tmp[PATH_MAX], cache[PATH_MAX];
    char *local_cache = getenv("REF_CACHE");
    mFILE *mf;

    if (fd->verbose)
	fprintf(stderr, "cram_populate_ref on fd %p, id %d\n", fd, id);

    if (!ref_path || *ref_path == '\0') {
	/*
	 * If we have no ref path, we use the EBI server.
	 * However to avoid spamming it we require a local ref cache too.
	 */
	ref_path = "http://www.ebi.ac.uk:80/ena/cram/md5/%s";
	if (!local_cache || *local_cache == '\0') {
	    const char *extra;
	    const char *base = get_cache_basedir(&extra);
#ifdef _WIN32
		_snprintf_s(cache,PATH_MAX, "%s%s/hts-ref/%%2s/%%2s/%%s", base, extra);
#else
	    snprintf(cache,PATH_MAX, "%s%s/hts-ref/%%2s/%%2s/%%s", base, extra);
#endif
	    local_cache = cache;
	    if (fd->verbose)
		fprintf(stderr, "Populating local cache: %s\n", local_cache);
	}
    }

    if (!r->name)
	return -1;

    if (!(ty = sam_hdr_find(fd->header, "SQ", "SN", r->name)))
	return -1;

    if (!(tag = sam_hdr_find_key(fd->header, ty, "M5", NULL)))
	goto no_M5;

    if (fd->verbose)
	fprintf(stderr, "Querying ref %s\n", tag->str+3);

    /* Use cache if available */
    if (local_cache && *local_cache) {
	struct stat sb;
	BGZF *fp;

	expand_cache_path(path, local_cache, tag->str+3);

	if (0 == stat(path, &sb) && (fp = bgzf_open(path, "r"))) {
	    r->length = sb.st_size;
	    r->offset = r->line_length = r->bases_per_line = 0;

	    r->fn = string_dup(fd->refs->pool, path);

	    if (fd->refs->fp)
		if (bgzf_close(fd->refs->fp) != 0)
		    return -1;
	    fd->refs->fp = fp;
	    fd->refs->fn = r->fn;

	    // Fall back to cram_get_ref() where it'll do the actual
	    // reading of the file.
	    return 0;
	}
    }

    /* Otherwise search */
    if ((mf = open_path_mfile(tag->str+3, ref_path, NULL))) {
	size_t sz;
	r->seq = (char*)mfsteal(mf, &sz);
	r->length = sz;
    } else {
	refs_t *refs;
	char *fn;

    no_M5:
	/* Failed to find in search path or M5 cache, see if @SQ UR: tag? */
	if (!(tag = sam_hdr_find_key(fd->header, ty, "UR", NULL)))
	    return -1;

	fn = (strncmp(tag->str+3, "file:", 5) == 0)
	    ? tag->str+8
	    : tag->str+3;

	if (fd->refs->fp) {
	    if (bgzf_close(fd->refs->fp) != 0)
		return -1;
	    fd->refs->fp = NULL;
	}
	if (!(refs = refs_load_fai(fd->refs, fn, 0)))
	    return -1;
	fd->refs = refs;
	if (fd->refs->fp) {
	    if (bgzf_close(fd->refs->fp) != 0)
		return -1;
	    fd->refs->fp = NULL;
	}

	if (!fd->refs->fn)
	    return -1;

	if (-1 == refs2id(fd->refs, fd->header))
	    return -1;
	if (!fd->refs->ref_id || !fd->refs->ref_id[id])
	    return -1;

	// Local copy already, so fall back to cram_get_ref().
	return 0;
    }

    /* Populate the local disk cache if required */
    if (local_cache && *local_cache) {
	FILE *fp;
	int i;

	expand_cache_path(path, local_cache, tag->str+3);
	if (fd->verbose)
	    fprintf(stderr, "Path='%s'\n", path);
	mkdir_prefix(path, 01777);

	i = 0;
	do {
	    sprintf(path_tmp, "%s.tmp_%d", path, /*getpid(),*/ i);
	    i++;
	    fp = fopen(path_tmp, "wx");
	} while (fp == NULL && errno == EEXIST);
	if (!fp) {
	    perror(path_tmp);

	    // Not fatal - we have the data already so keep going.
	    return 0;
	}

	if (r->length != fwrite(r->seq, 1, r->length, fp)) {
	    perror(path);
	}
	if (-1 == paranoid_fclose(fp)) {
	    unlink(path_tmp);
	} else {
	    if (0 == chmod(path_tmp, 0444))
		rename(path_tmp, path);
	    else
		unlink(path_tmp);
	}
    }

    return 0;
}

static void cram_ref_incr_locked(refs_t *r, int id) {
    RP("%d INC REF %d, %d %p\n", gettid(), id, (int)(id>=0?r->ref_id[id]->count+1:-999), id>=0?r->ref_id[id]->seq:(char *)1);

    if (id < 0 || !r->ref_id[id]->seq)
	return;

    if (r->last_id == id)
	r->last_id = -1;

    ++r->ref_id[id]->count;
}

static void cram_ref_decr_locked(refs_t *r, int id) {
    RP("%d DEC REF %d, %d %p\n", gettid(), id, (int)(id>=0?r->ref_id[id]->count-1:-999), id>=0?r->ref_id[id]->seq:(char *)1);

    if (id < 0 || !r->ref_id[id]->seq) {
	assert(r->ref_id[id]->count >= 0);
	return;
    }

    if (--r->ref_id[id]->count <= 0) {
	assert(r->ref_id[id]->count == 0);
	if (r->last_id >= 0) {
	    if (r->ref_id[r->last_id]->count <= 0 &&
		r->ref_id[r->last_id]->seq) {
		RP("%d FREE REF %d (%p)\n", gettid(),
		   r->last_id, r->ref_id[r->last_id]->seq);
		free(r->ref_id[r->last_id]->seq);
		r->ref_id[r->last_id]->seq = NULL;
		r->ref_id[r->last_id]->length = 0;
	    }
	}
	r->last_id = id;
    }
}


/*
 * Used by cram_ref_load and cram_ref_get. The file handle will have
 * already been opened, so we can catch it. The ref_entry *e informs us
 * of whether this is a multi-line fasta file or a raw MD5 style file.
 * Either way we create a single contiguous sequence.
 *
 * Returns all or part of a reference sequence on success (malloced);
 *         NULL on failure.
 */
static char *load_ref_portion(BGZF *fp, ref_entry *e, int start, int end) {
    off_t offset, len;
    char *seq;

    if (end < start)
	end = start;

    /*
     * Compute locations in file. This is trivial for the MD5 files, but
     * is still necessary for the fasta variants.
     */
    offset = e->line_length
	? e->offset + (start-1)/e->bases_per_line * e->line_length +
	  (start-1) % e->bases_per_line
	: start-1;

    len = (e->line_length
	   ? e->offset + (end-1)/e->bases_per_line * e->line_length + 
	     (end-1) % e->bases_per_line
	   : end-1) - offset + 1;

    if (bgzf_useek(fp, offset, SEEK_SET) < 0) {
	perror("bgzf_useek() on reference file");
	return NULL;
    }

	if (len == 0 || !(seq = (char*)malloc(len))) {
	return NULL;
    }

    if (len != bgzf_read(fp, seq, len)) {
	perror("bgzf_read() on reference file");
	free(seq);
	return NULL;
    }

    /* Strip white-space if required. */
    if (len != end-start+1) {
	int i, j;
	char *cp = seq;
	char *cp_to;

	for (i = j = 0; i < len; i++) {
	    if (cp[i] >= '!' && cp[i] <= '~')
		cp[j++] = cp[i] & ~0x20;
	}
	cp_to = cp+j;

	if (cp_to - seq != end-start+1) {
	    fprintf(stderr, "Malformed reference file?\n");
	    free(seq);
	    return NULL;
	}
    } else {
	int i;
	for (i = 0; i < len; i++) {
	    seq[i] = seq[i] & ~0x20; // uppercase in ASCII
	}
    }

    return seq;
}

/*
 * Load the entire reference 'id'.
 * This also increments the reference count by 1.
 *
 * Returns ref_entry on success;
 *         NULL on failure
 */
ref_entry *cram_ref_load(refs_t *r, int id) {
    ref_entry *e = r->ref_id[id];
    int start = 1, end = e->length;
    char *seq;

    if (e->seq) {
	return e;
    }

    assert(e->count == 0);

    if (r->last) {
#ifdef REF_DEBUG
	int idx = 0;
	for (idx = 0; idx < r->nref; idx++)
	    if (r->last == r->ref_id[idx])
		break;
	RP("%d cram_ref_load DECR %d\n", gettid(), idx);
#endif
	assert(r->last->count > 0);
	if (--r->last->count <= 0) {
	    RP("%d FREE REF %d (%p)\n", gettid(), id, r->ref_id[id]->seq);
	    if (r->last->seq) {
		free(r->last->seq);
		r->last->seq = NULL;
	    }
	}
    }

    /* Open file if it's not already the current open reference */
    if (strcmp(r->fn, e->fn) || r->fp == NULL) {
	if (r->fp)
	    if (bgzf_close(r->fp) != 0)
		return NULL;
	r->fn = e->fn;
	if (!(r->fp = bgzf_open_ref(r->fn, "r")))
	    return NULL;
    }

    RP("%d Loading ref %d (%d..%d)\n", gettid(), id, start, end);

    if (!(seq = load_ref_portion(r->fp, e, start, end))) {
	return NULL;
    }

    RP("%d Loaded ref %d (%d..%d) = %p\n", gettid(), id, start, end, seq);

    RP("%d INC REF %d, %d\n", gettid(), id, (int)(e->count+1));
    e->seq = seq;
    e->count++;

    /*
     * Also keep track of last used ref so incr/decr loops on the same
     * sequence don't cause load/free loops.
     */
    RP("%d cram_ref_load INCR %d => %d\n", gettid(), id, e->count+1);
    r->last = e;
    e->count++; 

    return e;
}

#if !defined(_WIN32) && defined(USE_PTHREAD)

void cram_ref_incr(refs_t *r, int id) {

	pthread_mutex_lock(&r->lock);
	cram_ref_incr_locked(r, id);
	pthread_mutex_unlock(&r->lock);

}

void cram_ref_decr(refs_t *r, int id) {

	pthread_mutex_lock(&r->lock);
	cram_ref_decr_locked(r, id);
	pthread_mutex_unlock(&r->lock);
}

/*
 * Returns a portion of a reference sequence from start to end inclusive.
 * The returned pointer is owned by either the cram_file fd or by the
 * internal refs_t structure and should not be freed  by the caller.
 *
 * The difference is whether or not this refs_t is in use by just the one
 * cram_fd or by multiples, or whether we have multiple threads accessing
 * references. In either case fd->shared will be true and we start using
 * reference counting to track the number of users of a specific reference
 * sequence.
 *
 * Otherwise the ref seq returned is allocated as part of cram_fd itself
 * and will be freed up on the next call to cram_get_ref or cram_close.
 *
 * To return the entire reference sequence, specify start as 1 and end
 * as 0.
 *
 * To cease using a reference, call cram_ref_decr().
 *
 * Returns reference on success,
 *         NULL on failure
 */
char *cram_get_ref(cram_fd *fd, int id, int start, int end) {
    ref_entry *r;
    char *seq;
    int ostart = start;

    if (id == -1)
	return NULL;

    /* FIXME: axiomatic query of r->seq being true?
     * Or shortcut for unsorted data where we load once and never free?
     */

    //fd->shared_ref = 1; // hard code for now to simplify things
    pthread_mutex_lock(&fd->ref_lock);
    RP("%d cram_get_ref on fd %p, id %d, range %d..%d\n", gettid(), fd, id, start, end);

    /*
     * Unsorted data implies we want to fetch an entire reference at a time.
     * We just deal with this at the moment by claiming we're sharing
     * references instead, which has the same requirement.
     */
    if (fd->unsorted)
	fd->shared_ref = 1;


    /* Sanity checking: does this ID exist? */
    if (id >= fd->refs->nref) {
	fprintf(stderr, "No reference found for id %d\n", id);
	pthread_mutex_unlock(&fd->ref_lock);
	return NULL;
    }

    if (!fd->refs || !fd->refs->ref_id[id]) {
	fprintf(stderr, "No reference found for id %d\n", id);
	pthread_mutex_unlock(&fd->ref_lock);
	return NULL;
    }

    if (!(r = fd->refs->ref_id[id])) {
	fprintf(stderr, "No reference found for id %d\n", id);
	pthread_mutex_unlock(&fd->ref_lock);
	return NULL;
    }


    /*
     * It has an entry, but may not have been populated yet.
     * Any manually loaded .fai files have their lengths known.
     * A ref entry computed from @SQ lines (M5 or UR field) will have
     * r->length == 0 unless it's been loaded once and verified that we have
     * an on-disk filename for it.
     *
     * 19 Sep 2013: Moved the lock here as the cram_populate_ref code calls
     * open_path_mfile and libcurl, which isn't multi-thread safe unless I
     * rewrite my code to have one curl handle per thread.
     */
    pthread_mutex_lock(&fd->refs->lock);
    if (r->length == 0) {
	if (cram_populate_ref(fd, id, r) == -1) {
	    fprintf(stderr, "Failed to populate reference for id %d\n", id);
	    pthread_mutex_unlock(&fd->refs->lock);
	    pthread_mutex_unlock(&fd->ref_lock);
	    return NULL;
	}
	r = fd->refs->ref_id[id];
	if (fd->unsorted)
	    cram_ref_incr_locked(fd->refs, id);
    }


    /*
     * We now know that we the filename containing the reference, so check
     * for limits. If it's over half the reference we'll load all of it in
     * memory as this will speed up subsequent calls.
     */
    if (end < 1)
	end = r->length;
    if (end >= r->length)
	end  = r->length; 
    assert(start >= 1);

    if (end - start >= 0.5*r->length || fd->shared_ref) {
	start = 1;
	end = r->length;
    }
    
    /*
     * Maybe we have it cached already? If so use it.
     *
     * Alternatively if we don't have the sequence but we're sharing
     * references and/or are asking for the entire length of it, then
     * load the full reference into the refs structure and return
     * a pointer to that one instead.
     */
    if (fd->shared_ref || r->seq || (start == 1 && end == r->length)) {
	char *cp;

	if (id >= 0) {
	    if (r->seq) {
		cram_ref_incr_locked(fd->refs, id);
	    } else {
		ref_entry *e;
		if (!(e = cram_ref_load(fd->refs, id))) {
			pthread_mutex_unlock(&fd->refs->lock);
			pthread_mutex_unlock(&fd->ref_lock);
		    return NULL;
		}

		/* unsorted data implies cache ref indefinitely, to avoid
		 * continually loading and unloading.
		 */
		if (fd->unsorted)
		    cram_ref_incr_locked(fd->refs, id);
	    }	    

	    fd->ref = NULL; /* We never access it directly */
	    fd->ref_start = 1;
	    fd->ref_end   = r->length;
	    fd->ref_id    = id;

	    cp = fd->refs->ref_id[id]->seq + ostart-1;
	} else {
	    fd->ref = NULL;
	    cp = NULL;
	}

	RP("%d cram_get_ref returning for id %d, count %d\n", gettid(), id, (int)r->count);

	pthread_mutex_unlock(&fd->refs->lock);
	pthread_mutex_unlock(&fd->ref_lock);
	return cp;
    }

    /*
     * Otherwise we're not sharing, we don't have a copy of it already and
     * we're only asking for a small portion of it.
     *
     * In this case load up just that segment ourselves, freeing any old
     * small segments in the process.
     */

    /* Unmapped ref ID */
    if (id < 0) {
	if (fd->ref_free) {
	    free(fd->ref_free);
	    fd->ref_free = NULL;
	}
	fd->ref = NULL;
	fd->ref_id = id;
	pthread_mutex_unlock(&fd->refs->lock);
	pthread_mutex_unlock(&fd->ref_lock);
	return NULL;
    }

    /* Open file if it's not already the current open reference */
    if (strcmp(fd->refs->fn, r->fn) || fd->refs->fp == NULL) {
	if (fd->refs->fp)
	    if (bgzf_close(fd->refs->fp) != 0)
		return NULL;
	fd->refs->fn = r->fn;
	if (!(fd->refs->fp = bgzf_open_ref(fd->refs->fn, "r"))) {
		pthread_mutex_unlock(&fd->refs->lock);
		pthread_mutex_unlock(&fd->ref_lock);
	    return NULL;
	}
    }

    if (!(fd->ref = load_ref_portion(fd->refs->fp, r, start, end))) {
		pthread_mutex_unlock(&fd->refs->lock);
		pthread_mutex_unlock(&fd->ref_lock);
	return NULL;
    }

    if (fd->ref_free)
	free(fd->ref_free);

    fd->ref_id    = id;
    fd->ref_start = start;
    fd->ref_end   = end;
    fd->ref_free = fd->ref;
    seq = fd->ref;

	pthread_mutex_unlock(&fd->refs->lock);
	pthread_mutex_unlock(&fd->ref_lock);

    return seq + ostart - start;
}

#else // ~ !defined(_WIN32) && defined(USE_PTHREAD)

void cram_ref_incr(refs_t *r, int id) {

	r->lock.lock();
	cram_ref_incr_locked(r, id);
	r->lock.unlock();

}

void cram_ref_decr(refs_t *r, int id) {

	r->lock.lock();
	cram_ref_decr_locked(r, id);
	r->lock.unlock();
}

/*
* Returns a portion of a reference sequence from start to end inclusive.
* The returned pointer is owned by either the cram_file fd or by the
* internal refs_t structure and should not be freed  by the caller.
*
* The difference is whether or not this refs_t is in use by just the one
* cram_fd or by multiples, or whether we have multiple threads accessing
* references. In either case fd->shared will be true and we start using
* reference counting to track the number of users of a specific reference
* sequence.
*
* Otherwise the ref seq returned is allocated as part of cram_fd itself
* and will be freed up on the next call to cram_get_ref or cram_close.
*
* To return the entire reference sequence, specify start as 1 and end
* as 0.
*
* To cease using a reference, call cram_ref_decr().
*
* Returns reference on success,
*         NULL on failure
*/
char *cram_get_ref(cram_fd *fd, int id, int start, int end) {
	ref_entry *r;
	char *seq;
	int ostart = start;

	if (id == -1)
		return NULL;

	/* FIXME: axiomatic query of r->seq being true?
	* Or shortcut for unsorted data where we load once and never free?
	*/

	//fd->shared_ref = 1; // hard code for now to simplify things
	fd->ref_lock.lock();
	RP("%d cram_get_ref on fd %p, id %d, range %d..%d\n", gettid(), fd, id, start, end);

	/*
	* Unsorted data implies we want to fetch an entire reference at a time.
	* We just deal with this at the moment by claiming we're sharing
	* references instead, which has the same requirement.
	*/
	if (fd->unsorted)
		fd->shared_ref = 1;


	/* Sanity checking: does this ID exist? */
	if (id >= fd->refs->nref) {
		fprintf(stderr, "No reference found for id %d\n", id);
		fd->ref_lock.unlock();
		return NULL;
	}

	if (!fd->refs || !fd->refs->ref_id[id]) {
		fprintf(stderr, "No reference found for id %d\n", id);
		fd->ref_lock.unlock();
		return NULL;
	}

	if (!(r = fd->refs->ref_id[id])) {
		fprintf(stderr, "No reference found for id %d\n", id);
		fd->ref_lock.unlock();
		return NULL;
	}


	/*
	* It has an entry, but may not have been populated yet.
	* Any manually loaded .fai files have their lengths known.
	* A ref entry computed from @SQ lines (M5 or UR field) will have
	* r->length == 0 unless it's been loaded once and verified that we have
	* an on-disk filename for it.
	*
	* 19 Sep 2013: Moved the lock here as the cram_populate_ref code calls
	* open_path_mfile and libcurl, which isn't multi-thread safe unless I
	* rewrite my code to have one curl handle per thread.
	*/
	fd->refs->lock.lock();
	if (r->length == 0) {
		if (cram_populate_ref(fd, id, r) == -1) {
			fprintf(stderr, "Failed to populate reference for id %d\n", id);
			fd->refs->lock.unlock();
			fd->ref_lock.unlock();
			return NULL;
		}
		r = fd->refs->ref_id[id];
		if (fd->unsorted)
			cram_ref_incr_locked(fd->refs, id);
	}


	/*
	* We now know that we the filename containing the reference, so check
	* for limits. If it's over half the reference we'll load all of it in
	* memory as this will speed up subsequent calls.
	*/
	if (end < 1)
		end = r->length;
	if (end >= r->length)
		end = r->length;
	assert(start >= 1);

	if (end - start >= 0.5*r->length || fd->shared_ref) {
		start = 1;
		end = r->length;
	}

	/*
	* Maybe we have it cached already? If so use it.
	*
	* Alternatively if we don't have the sequence but we're sharing
	* references and/or are asking for the entire length of it, then
	* load the full reference into the refs structure and return
	* a pointer to that one instead.
	*/
	if (fd->shared_ref || r->seq || (start == 1 && end == r->length)) {
		char *cp;

		if (id >= 0) {
			if (r->seq) {
				cram_ref_incr_locked(fd->refs, id);
			}
			else {
				ref_entry *e;
				if (!(e = cram_ref_load(fd->refs, id))) {
					fd->refs->lock.unlock();
					fd->ref_lock.unlock();
					return NULL;
				}

				/* unsorted data implies cache ref indefinitely, to avoid
				* continually loading and unloading.
				*/
				if (fd->unsorted)
					cram_ref_incr_locked(fd->refs, id);
			}

			fd->ref = NULL; /* We never access it directly */
			fd->ref_start = 1;
			fd->ref_end = r->length;
			fd->ref_id = id;

			cp = fd->refs->ref_id[id]->seq + ostart - 1;
		}
		else {
			fd->ref = NULL;
			cp = NULL;
		}

		RP("%d cram_get_ref returning for id %d, count %d\n", gettid(), id, (int)r->count);

		fd->refs->lock.unlock();
		fd->ref_lock.unlock();
		return cp;
	}

	/*
	* Otherwise we're not sharing, we don't have a copy of it already and
	* we're only asking for a small portion of it.
	*
	* In this case load up just that segment ourselves, freeing any old
	* small segments in the process.
	*/

	/* Unmapped ref ID */
	if (id < 0) {
		if (fd->ref_free) {
			free(fd->ref_free);
			fd->ref_free = NULL;
		}
		fd->ref = NULL;
		fd->ref_id = id;
		fd->refs->lock.unlock();
		fd->ref_lock.unlock();
		return NULL;
	}

	/* Open file if it's not already the current open reference */
	if (strcmp(fd->refs->fn, r->fn) || fd->refs->fp == NULL) {
		if (fd->refs->fp)
			if (bgzf_close(fd->refs->fp) != 0)
				return NULL;
		fd->refs->fn = r->fn;
		if (!(fd->refs->fp = bgzf_open_ref(fd->refs->fn, "r"))) {
			fd->refs->lock.unlock();
			fd->ref_lock.unlock();
			return NULL;
		}
	}

	if (!(fd->ref = load_ref_portion(fd->refs->fp, r, start, end))) {
		fd->refs->lock.unlock();
		fd->ref_lock.unlock();
		return NULL;
	}

	if (fd->ref_free)
		free(fd->ref_free);

	fd->ref_id = id;
	fd->ref_start = start;
	fd->ref_end = end;
	fd->ref_free = fd->ref;
	seq = fd->ref;

	fd->refs->lock.unlock();
	fd->ref_lock.unlock();

	return seq + ostart - start;
}

#endif // ~ !defined(_WIN32) && defined(USE_PTHREAD)
/*
 * If fd has been opened for reading, it may be permitted to specify 'fn'
 * as NULL and let the code auto-detect the reference by parsing the
 * SAM header @SQ lines.
 */
int cram_load_reference(cram_fd *fd, char *fn) {
    if (fn) {
	fd->refs = refs_load_fai(fd->refs, fn,
				 !(fd->embed_ref && fd->mode == 'r'));
	fn = fd->refs ? fd->refs->fn : NULL;
    }
    fd->ref_fn = fn;

    if ((!fd->refs || (fd->refs->nref == 0 && !fn)) && fd->header) {
	if (fd->refs)
	    refs_free(fd->refs);
	if (!(fd->refs = refs_create()))
	    return -1;
	if (-1 == refs_from_header(fd->refs, fd, fd->header))
	    return -1;
    }

    if (fd->header)
	if (-1 == refs2id(fd->refs, fd->header))
	    return -1;

    return fn ? 0 : -1;
}

/* ----------------------------------------------------------------------
 * Containers
 */

/*
 * Creates a new container, specifying the maximum number of slices
 * and records permitted.
 *
 * Returns cram_container ptr on success
 *         NULL on failure
 */
cram_container *cram_new_container(int nrec, int nslice) {
	cram_container *c = (cram_container*)calloc(1, sizeof(*c));
    enum cram_DS_ID id;

    if (!c)
	return NULL;

    c->curr_ref = -2;

    c->max_c_rec = nrec * nslice;
    c->curr_c_rec = 0;

    c->max_rec = nrec;
    c->record_counter = 0;
    c->num_bases = 0;

    c->max_slice = nslice;
    c->curr_slice = 0;

    c->pos_sorted = 1;
    c->max_apos   = 0;
    c->multi_seq  = 0;

    c->bams = NULL;

    if (!(c->slices = (cram_slice **)calloc(nslice, sizeof(cram_slice *))))
	goto err;
    c->slice = NULL;

    if (!(c->comp_hdr = cram_new_compression_header()))
	goto err;
    c->comp_hdr_block = NULL;

    for (int id = DS_RN; id < DS_TN; id++)
	if (!(c->stats[id] = cram_stats_create())) goto err;
    
    //c->aux_B_stats = cram_stats_create();

    if (!(c->tags_used = kh_init(s_i2i)))
	goto err;
    c->refs_used = 0;

    return c;

 err:
    if (c) {
	if (c->slices)
	    free(c->slices);
	free(c);
    }
    return NULL;
}

void cram_free_container(cram_container *c) {
    enum cram_DS_ID id;
    int i;

    if (!c)
	return;

    if (c->refs_used)
	free(c->refs_used);

    if (c->landmark)
	free(c->landmark);

    if (c->comp_hdr)
	cram_free_compression_header(c->comp_hdr);

    if (c->comp_hdr_block)
	cram_free_block(c->comp_hdr_block);

    if (c->slices) {
	for (i = 0; i < c->max_slice; i++)
	    if (c->slices[i])
		cram_free_slice(c->slices[i]);
	free(c->slices);
    }

    for (int id = DS_RN; id < DS_TN; id++)
	if (c->stats[id]) cram_stats_free(c->stats[id]);

    //if (c->aux_B_stats) cram_stats_free(c->aux_B_stats);
    
    if (c->tags_used) kh_destroy(s_i2i, c->tags_used);

    free(c);
}

/*
 * Reads a container header.
 *
 * Returns cram_container on success
 *         NULL on failure or no container left (fd->err == 0).
 */
cram_container *cram_read_container(cram_fd *fd) {
    cram_container c2, *c;
    int i, s;
    size_t rd = 0;
    
    fd->err = 0;
    fd->eof = 0;

    memset(&c2, 0, sizeof(c2));
    if (CRAM_MAJOR_VERS(fd->version) == 1) {
	if ((s = itf8_decode(fd, &c2.length)) == -1) {
	    fd->eof = fd->empty_container ? 1 : 2;
	    return NULL;
	} else {
	    rd+=s;
	}
    } else {
	if ((s = int32_decode(fd, &c2.length)) == -1) {
	    if (CRAM_MAJOR_VERS(fd->version) == 2 &&
		CRAM_MINOR_VERS(fd->version) == 0)
		fd->eof = 1; // EOF blocks arrived in v2.1
	    else
		fd->eof = fd->empty_container ? 1 : 2;
	    return NULL;
	} else {
	    rd+=s;
	}
    }
    if ((s = itf8_decode(fd, &c2.ref_seq_id))   == -1) return NULL; else rd+=s;
    if ((s = itf8_decode(fd, &c2.ref_seq_start))== -1) return NULL; else rd+=s;
    if ((s = itf8_decode(fd, &c2.ref_seq_span)) == -1) return NULL; else rd+=s;
    if ((s = itf8_decode(fd, &c2.num_records))  == -1) return NULL; else rd+=s;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
	c2.record_counter = 0;
	c2.num_bases = 0;
    } else {
	if (CRAM_MAJOR_VERS(fd->version) >= 3) {
	    if ((s = ltf8_decode(fd, &c2.record_counter)) == -1)
		return NULL;
	    else
		rd += s;
	} else {
	    int32_t i32;
	    if ((s = itf8_decode(fd, &i32)) == -1)
		return NULL;
	    else
		rd += s;
	    c2.record_counter = i32;
	}

	if ((s = ltf8_decode(fd, &c2.num_bases))== -1)
	    return NULL;
	else
	    rd += s;
    }
    if ((s = itf8_decode(fd, &c2.num_blocks))   == -1) return NULL; else rd+=s;
    if ((s = itf8_decode(fd, &c2.num_landmarks))== -1) return NULL; else rd+=s;

    if (!(c = (cram_container*)calloc(1, sizeof(*c))))
	return NULL;

    *c = c2;

    if (!(c->landmark = (int32_t*)malloc(c->num_landmarks * sizeof(int32_t))) &&
	c->num_landmarks) {
	fd->err = errno;
	cram_free_container(c);
	return NULL;
    }  
    for (i = 0; i < c->num_landmarks; i++) {
	if ((s = itf8_decode(fd, &c->landmark[i])) == -1) {
	    cram_free_container(c);
	    return NULL;
	} else {
	    rd += s;
	}
    }

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
	uint32_t crc, i;
	unsigned char *dat = (unsigned char*)malloc(50 + 5 * (c->num_landmarks)), *cp = dat;
	if (!dat) {
	    cram_free_container(c);
	    return NULL;
	}
	if (-1 == int32_decode(fd, (int32_t *)&c->crc32))
	    return NULL;
	else
	    rd+=4;

	/* Reencode first as we can't easily access the original byte stream.
	 *
	 * FIXME: Technically this means this may not be fool proof. We could
	 * create a CRAM file using a 2 byte ITF8 value that can fit in a 
	 * 1 byte field, meaning the encoding is different to the original
	 * form and so has a different CRC.
	 *
	 * The correct implementation would be to have an alternative form
	 * of itf8_decode which also squirrels away the raw byte stream
	 * during decoding so we can then CRC that.
	 */
	*(unsigned int *)cp = le_int4(c->length); cp += 4;
	cp += itf8_put(cp, c->ref_seq_id);
	cp += itf8_put(cp, c->ref_seq_start);
	cp += itf8_put(cp, c->ref_seq_span);
	cp += itf8_put(cp, c->num_records);
	cp += ltf8_put((char *)cp, c->record_counter);
	cp += itf8_put(cp, c->num_bases);
	cp += itf8_put(cp, c->num_blocks);
	cp += itf8_put(cp, c->num_landmarks);
	for (i = 0; i < c->num_landmarks; i++) {
	    cp += itf8_put(cp, c->landmark[i]);
	}

	crc = crc32(0L, dat, cp-dat);
	if (crc != c->crc32) {
	    fprintf(stderr, "Container header CRC32 failure\n");
	    cram_free_container(c);
	    return NULL;
	}
    }

    c->offset = rd;
    c->slices = NULL;
    c->curr_slice = 0;
    c->max_slice = c->num_landmarks;
    c->slice_rec = 0;
    c->curr_rec = 0;
    c->max_rec = 0;

    if (c->ref_seq_id == -2) {
	c->multi_seq = 1;
	fd->multi_seq = 1;
    }

    fd->empty_container =
	(c->num_records == 0 &&
	 c->ref_seq_id == -1 &&
	 c->ref_seq_start == 0x454f46 /* EOF */) ? 1 : 0;

    return c;
}

/*
 * Writes a container structure.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_container(cram_fd *fd, cram_container *c) {
    char buf_a[1024], *buf = buf_a, *cp;
    int i;

    if (55 + c->num_landmarks * 5 >= 1024)
		buf = (char*)malloc(55 + c->num_landmarks * 5);
    cp = buf;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
	cp += itf8_put(cp, c->length);
    } else {
	*(int32_t *)cp = le_int4(c->length);
	cp += 4;
    }
    if (c->multi_seq) {
	cp += itf8_put(cp, -2);
	cp += itf8_put(cp, 0);
	cp += itf8_put(cp, 0);
    } else {
	cp += itf8_put(cp, c->ref_seq_id);
	cp += itf8_put(cp, c->ref_seq_start);
	cp += itf8_put(cp, c->ref_seq_span);
    }
    cp += itf8_put(cp, c->num_records);
    if (CRAM_MAJOR_VERS(fd->version) == 2) {
	cp += itf8_put(cp, c->record_counter);
	cp += ltf8_put(cp, c->num_bases);
    } else if (CRAM_MAJOR_VERS(fd->version) >= 3) {
	cp += ltf8_put(cp, c->record_counter);
	cp += ltf8_put(cp, c->num_bases);
    }

    cp += itf8_put(cp, c->num_blocks);
    cp += itf8_put(cp, c->num_landmarks);
    for (i = 0; i < c->num_landmarks; i++)
	cp += itf8_put(cp, c->landmark[i]);

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
	c->crc32 = crc32(0L, (uc *)buf, cp-buf);
	cp[0] =  c->crc32        & 0xff;
	cp[1] = (c->crc32 >>  8) & 0xff;
	cp[2] = (c->crc32 >> 16) & 0xff;
	cp[3] = (c->crc32 >> 24) & 0xff;
	cp += 4;
    }

    if (cp-buf != hwrite(fd->fp, buf, cp-buf)) {
	if (buf != buf_a)
	    free(buf);
	return -1;
    }

    if (buf != buf_a)
	free(buf);

    return 0;
}

// common component shared by cram_flush_container{,_mt}
static int cram_flush_container2(cram_fd *fd, cram_container *c) {
    int i, j;

    //fprintf(stderr, "Writing container %d, sum %u\n", c->record_counter, sum);

    /* Write the container struct itself */
    if (0 != cram_write_container(fd, c))
	return -1;

    /* And the compression header */
    if (0 != cram_write_block(fd, c->comp_hdr_block))
	return -1;

    /* Followed by the slice blocks */
    for (i = 0; i < c->curr_slice; i++) {
	cram_slice *s = c->slices[i];

	if (0 != cram_write_block(fd, s->hdr_block))
	    return -1;

	for (j = 0; j < s->hdr->num_blocks; j++) {
	    if (0 != cram_write_block(fd, s->block[j]))
		return -1;
	}
    }

    return hflush(fd->fp) == 0 ? 0 : -1;
}

/*
 * Flushes a completely or partially full container to disk, writing
 * container structure, header and blocks. This also calls the encoder
 * functions.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_flush_container(cram_fd *fd, cram_container *c) {
    /* Encode the container blocks and generate compression header */
    if (0 != cram_encode_container(fd, c))
	return -1;

    return cram_flush_container2(fd, c);
}

typedef struct {
    cram_fd *fd;
    cram_container *c;
} cram_job;

void *cram_flush_thread(void *arg) {
    cram_job *j = (cram_job *)arg;

    /* Encode the container blocks and generate compression header */
    if (0 != cram_encode_container(j->fd, j->c)) {
	fprintf(stderr, "cram_encode_container failed\n");
	return NULL;
    }

    return arg;
}

static int cram_flush_result(cram_fd *fd) {
    int i, ret = 0;
    t_pool_result *r;

    while ((r = t_pool_next_result(fd->rqueue))) {
	cram_job *j = (cram_job *)r->data;
	cram_container *c;

	if (!j) {
	    t_pool_delete_result(r, 0);
	    return -1;
	}

	fd = j->fd;
	c = j->c;

	if (0 != cram_flush_container2(fd, c))
	    return -1;

	/* Free the container */
	for (i = 0; i < c->max_slice; i++) {
	    cram_free_slice(c->slices[i]);
	    c->slices[i] = NULL;
	}

	c->slice = NULL;
	c->curr_slice = 0;

	cram_free_container(c);

	ret |= hflush(fd->fp) == 0 ? 0 : -1;

	t_pool_delete_result(r, 1);
    }

    return ret;
}

int cram_flush_container_mt(cram_fd *fd, cram_container *c) {
    cram_job *j;

    if (!fd->pool)
	return cram_flush_container(fd, c);

	if (!(j = (cram_job*)malloc(sizeof(*j))))
	return -1;
    j->fd = fd;
    j->c = c;
    
    t_pool_dispatch(fd->pool, fd->rqueue, cram_flush_thread, j);

    return cram_flush_result(fd);
}

/* ----------------------------------------------------------------------
 * Compression headers; the first part of the container
 */

/*
 * Creates a new blank container compression header
 *
 * Returns header ptr on success
 *         NULL on failure
 */
cram_block_compression_hdr *cram_new_compression_header(void) {
	cram_block_compression_hdr *hdr = (cram_block_compression_hdr*)calloc(1, sizeof(*hdr));
    if (!hdr)
	return NULL;

    if (!(hdr->TD_blk = cram_new_block(CORE, 0))) {
	free(hdr);
	return NULL;
    }

    if (!(hdr->TD_hash = kh_init(m_s2i))) {
	cram_free_block(hdr->TD_blk);
	free(hdr);
	return NULL;
    }

    if (!(hdr->TD_keys = string_pool_create(8192))) {
	kh_destroy(m_s2i, hdr->TD_hash);
	cram_free_block(hdr->TD_blk);
	free(hdr);
	return NULL;
    }

    return hdr;
}

void cram_free_compression_header(cram_block_compression_hdr *hdr) {
    int i;

    if (hdr->landmark)
	free(hdr->landmark);

    if (hdr->preservation_map)
	kh_destroy(map, hdr->preservation_map);

    for (i = 0; i < CRAM_MAP_HASH; i++) {
	cram_map *m, *m2;
	for (m = hdr->rec_encoding_map[i]; m; m = m2) {
	    m2 = m->next;
	    if (m->codec)
		m->codec->free(m->codec);
	    free(m);
	}
    }

    for (i = 0; i < CRAM_MAP_HASH; i++) {
	cram_map *m, *m2;
	for (m = hdr->tag_encoding_map[i]; m; m = m2) {
	    m2 = m->next;
	    if (m->codec)
		m->codec->free(m->codec);
	    free(m);
	}
    }

    for (i = 0; i < DS_END; i++) {
	if (hdr->codecs[i])
	    hdr->codecs[i]->free(hdr->codecs[i]);
    }

    if (hdr->TL)
	free(hdr->TL);
    if (hdr->TD_blk)
	cram_free_block(hdr->TD_blk);
    if (hdr->TD_hash)
	kh_destroy(m_s2i, hdr->TD_hash);
    if (hdr->TD_keys)
	string_pool_destroy(hdr->TD_keys);

    free(hdr);
}


/* ----------------------------------------------------------------------
 * Slices and slice headers
 */

void cram_free_slice_header(cram_block_slice_hdr *hdr) {
    if (!hdr)
	return;

    if (hdr->block_content_ids)
	free(hdr->block_content_ids);

    free(hdr);

    return;
}

void cram_free_slice(cram_slice *s) {
    if (!s)
	return;

    if (s->hdr_block)
	cram_free_block(s->hdr_block);

    if (s->block) {
	int i;

	if (s->hdr) {
	    for (i = 0; i < s->hdr->num_blocks; i++) {
		cram_free_block(s->block[i]);
	    }
	}
	free(s->block);
    }

    if (s->block_by_id)
	free(s->block_by_id);

    if (s->hdr)
	cram_free_slice_header(s->hdr);

    if (s->seqs_blk)
	cram_free_block(s->seqs_blk);

    if (s->qual_blk)
	cram_free_block(s->qual_blk);

    if (s->name_blk)
	cram_free_block(s->name_blk);

    if (s->aux_blk)
	cram_free_block(s->aux_blk);

    if (s->aux_OQ_blk)
	cram_free_block(s->aux_OQ_blk);

    if (s->aux_BQ_blk)
	cram_free_block(s->aux_BQ_blk);

    if (s->aux_FZ_blk)
	cram_free_block(s->aux_FZ_blk);

    if (s->aux_oq_blk)
	cram_free_block(s->aux_oq_blk);

    if (s->aux_os_blk)
	cram_free_block(s->aux_os_blk);

    if (s->aux_oz_blk)
	cram_free_block(s->aux_oz_blk);

    if (s->base_blk)
	cram_free_block(s->base_blk);

    if (s->soft_blk)
	cram_free_block(s->soft_blk);

    if (s->cigar)
	free(s->cigar);

    if (s->crecs)
	free(s->crecs);

    if (s->features)
	free(s->features);

    if (s->TN)
	free(s->TN);

    if (s->pair_keys)
	string_pool_destroy(s->pair_keys);

    if (s->pair[0])
	kh_destroy(m_s2i, s->pair[0]);
    if (s->pair[1])
	kh_destroy(m_s2i, s->pair[1]);

    free(s);
}

/*
 * Creates a new empty slice in memory, for subsequent writing to
 * disk.
 *
 * Returns cram_slice ptr on success
 *         NULL on failure
 */
cram_slice *cram_new_slice(enum cram_content_type type, int nrecs) {
	cram_slice *s = (cram_slice*)calloc(1, sizeof(*s));
    if (!s)
	return NULL;

    if (!(s->hdr = (cram_block_slice_hdr *)calloc(1, sizeof(*s->hdr))))
	goto err;
    s->hdr->content_type = type;

    s->hdr_block = NULL;
    s->block = NULL;
    s->block_by_id = NULL;
    s->last_apos = 0;
	if (!(s->crecs = (cram_record*)malloc(nrecs * sizeof(cram_record))))  goto err;
    s->cigar = NULL;
    s->cigar_alloc = 0;
    s->ncigar = 0;

    if (!(s->seqs_blk = cram_new_block(EXTERNAL, 0)))       goto err;
    if (!(s->qual_blk = cram_new_block(EXTERNAL, DS_QS)))   goto err;
    if (!(s->name_blk = cram_new_block(EXTERNAL, DS_RN)))   goto err;
    if (!(s->aux_blk  = cram_new_block(EXTERNAL, DS_aux)))  goto err;
    if (!(s->base_blk = cram_new_block(EXTERNAL, DS_IN)))   goto err;
    if (!(s->soft_blk = cram_new_block(EXTERNAL, DS_SC)))   goto err;

    s->features = NULL;
    s->nfeatures = s->afeatures = 0;

#ifndef TN_external
    s->TN = NULL;
    s->nTN = s->aTN = 0;
#endif

    // Volatile keys as we do realloc in dstring
    if (!(s->pair_keys = string_pool_create(8192))) goto err;
    if (!(s->pair[0] = kh_init(m_s2i)))             goto err;
    if (!(s->pair[1] = kh_init(m_s2i)))             goto err;
    
#ifdef BA_external
    s->BA_len = 0;
#endif

    return s;

 err:
    if (s)
	cram_free_slice(s);

    return NULL;
}

/*
 * Loads an entire slice.
 * FIXME: In 1.0 the native unit of slices within CRAM is broken
 * as slices contain references to objects in other slices.
 * To work around this while keeping the slice oriented outer loop
 * we read all slices and stitch them together into a fake large
 * slice instead.
 *
 * Returns cram_slice ptr on success
 *         NULL on failure
 */
cram_slice *cram_read_slice(cram_fd *fd) {
    cram_block *b = cram_read_block(fd);
	cram_slice *s = (cram_slice*)calloc(1, sizeof(*s));
    int i, n, max_id, min_id;

    if (!b || !s)
	goto err;

    s->hdr_block = b;
    switch (b->content_type) {
    case MAPPED_SLICE:
    case UNMAPPED_SLICE:
	if (!(s->hdr = cram_decode_slice_header(fd, b)))
	    goto err;
	break;

    default:
	fprintf(stderr, "Unexpected block of type %s\n",
		cram_content_type2str(b->content_type));
	goto err;
    }

    s->block = (cram_block**)calloc(n = s->hdr->num_blocks, sizeof(*s->block));
    if (!s->block)
	goto err;

    for (max_id = i = 0, min_id = INT_MAX; i < n; i++) {
	if (!(s->block[i] = cram_read_block(fd)))
	    goto err;

	if (s->block[i]->content_type == EXTERNAL) {
	    if (max_id < s->block[i]->content_id)
		max_id = s->block[i]->content_id;
	    if (min_id > s->block[i]->content_id)
		min_id = s->block[i]->content_id;
	}
    }
    if (min_id >= 0 && max_id < 1024) {
	if (!(s->block_by_id = (cram_block**)calloc(1024, sizeof(s->block[0]))))
	    goto err;

	for (i = 0; i < n; i++) {
	    if (s->block[i]->content_type != EXTERNAL)
		continue;
	    s->block_by_id[s->block[i]->content_id] = s->block[i];
	}
    }

    /* Initialise encoding/decoding tables */
    s->cigar = NULL;
    s->cigar_alloc = 0;
    s->ncigar = 0;

    if (!(s->seqs_blk = cram_new_block(EXTERNAL, 0)))      goto err;
    if (!(s->qual_blk = cram_new_block(EXTERNAL, DS_QS)))  goto err;
    if (!(s->name_blk = cram_new_block(EXTERNAL, DS_RN)))  goto err;
    if (!(s->aux_blk  = cram_new_block(EXTERNAL, DS_aux))) goto err;
    if (!(s->base_blk = cram_new_block(EXTERNAL, DS_IN)))  goto err;
    if (!(s->soft_blk = cram_new_block(EXTERNAL, DS_SC)))  goto err;

    s->crecs = NULL;

    s->last_apos = s->hdr->ref_seq_start;
    
    return s;

 err:
    if (b)
	cram_free_block(b);
    if (s) {
	s->hdr_block = NULL;
	cram_free_slice(s);
    }
    return NULL;
}


/* ----------------------------------------------------------------------
 * CRAM file definition (header)
 */

/*
 * Reads a CRAM file definition structure.
 * Returns file_def ptr on success
 *         NULL on failure
 */
cram_file_def *cram_read_file_def(cram_fd *fd) {
	cram_file_def *def = (cram_file_def*)malloc(sizeof(*def));
    if (!def)
	return NULL;

    if (26 != hread(fd->fp, &def->magic[0], 26)) {
	free(def);
	return NULL;
    }

    if (memcmp(def->magic, "CRAM", 4) != 0) {
	free(def);
	return NULL;
    }

    if (def->major_version > 3) {
	fprintf(stderr, "CRAM version number mismatch\n"
		"Expected 1.x, 2.x or 3.x, got %d.%d\n",
		def->major_version, def->minor_version);
	free(def);
	return NULL;
    }

    fd->first_container += 26;
    fd->last_slice = 0;

    return def;
}

/*
 * Writes a cram_file_def structure to cram_fd.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_file_def(cram_fd *fd, cram_file_def *def) {
    return (hwrite(fd->fp, &def->magic[0], 26) == 26) ? 0 : -1;
}

void cram_free_file_def(cram_file_def *def) {
    if (def) free(def);
}

/* ----------------------------------------------------------------------
 * SAM header I/O
 */


/*
 * Reads the SAM header from the first CRAM data block.
 * Also performs minimal parsing to extract read-group
 * and sample information.

 * Returns SAM hdr ptr on success
 *         NULL on failure
 */
SAM_hdr *cram_read_SAM_hdr(cram_fd *fd) {
    int32_t header_len;
    char *header;
    SAM_hdr *hdr;

    /* 1.1 onwards stores the header in the first block of a container */
    if (CRAM_MAJOR_VERS(fd->version) == 1) {
	/* Length */
	if (-1 == int32_decode(fd, &header_len))
	    return NULL;

	/* Alloc and read */
	if (NULL == (header = (char*)malloc(header_len + 1)))
	    return NULL;

	*header = 0;
	if (header_len != hread(fd->fp, header, header_len))
	    return NULL;

	fd->first_container += 4 + header_len;
    } else {
	cram_container *c = cram_read_container(fd);
	cram_block *b;
	int i, len;

	if (!c)
	    return NULL;

	if (c->num_blocks < 1) {
	    cram_free_container(c);
	    return NULL;
	}

	if (!(b = cram_read_block(fd))) {
	    cram_free_container(c);
	    return NULL;
	}
	cram_uncompress_block(b);

	len = b->comp_size + 2 + 4*(CRAM_MAJOR_VERS(fd->version) >= 3) +
	    itf8_size(b->content_id) + 
	    itf8_size(b->uncomp_size) + 
	    itf8_size(b->comp_size);

	/* Extract header from 1st block */
	if (-1 == int32_get(b, &header_len) ||
	    b->uncomp_size - 4 < header_len) {
	    cram_free_container(c);
	    cram_free_block(b);
	    return NULL;
	}
	if (NULL == (header = (char*)malloc(header_len + 1))) {
	    cram_free_container(c);
	    cram_free_block(b);
	    return NULL;
	}
	memcpy(header, BLOCK_END(b), header_len);
	header[header_len]='\0';
	cram_free_block(b);

	/* Consume any remaining blocks */
	for (i = 1; i < c->num_blocks; i++) {
	    if (!(b = cram_read_block(fd))) {
		cram_free_container(c);
		return NULL;
	    }
	    len += b->comp_size + 2 + 4*(CRAM_MAJOR_VERS(fd->version) >= 3) +
		itf8_size(b->content_id) + 
		itf8_size(b->uncomp_size) + 
		itf8_size(b->comp_size);
	    cram_free_block(b);
	}

	if (c->length && c->length > len) {
	    // Consume padding
	    char *pads = (char*)malloc(c->length - len);
	    if (!pads) {
		cram_free_container(c);
		return NULL;
	    }

	    if (c->length - len != hread(fd->fp, pads, c->length - len)) {
		cram_free_container(c);
		return NULL;
	    }
	    free(pads);
	}

	cram_free_container(c);
    }

    /* Parse */
    hdr = sam_hdr_parse_(header, header_len);
    free(header);

    return hdr;
}

/*
 * Converts 'in' to a full pathname to store in out.
 * Out must be at least PATH_MAX bytes long.
 */
static void full_path(char *out, char *in) {
    if (*in == '/') {
	strncpy(out, in, PATH_MAX);
	out[PATH_MAX-1] = 0;
    } else {
	int len;

	// unable to get dir or out+in is too long
	if (!getcwd(out, PATH_MAX) ||
	    (len = strlen(out))+1+strlen(in) >= PATH_MAX) {
	    strncpy(out, in, PATH_MAX);
	    out[PATH_MAX-1] = 0;
	    return;
	}

	sprintf(out+len, "/%.*s", PATH_MAX - len, in);

	// FIXME: cope with `pwd`/../../../foo.fa ?
    }
}

/*
 * Writes a CRAM SAM header.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_write_SAM_hdr(cram_fd *fd, SAM_hdr *hdr) {
    int header_len;
    int blank_block = (CRAM_MAJOR_VERS(fd->version) >= 3);

    /* Write CRAM MAGIC if not yet written. */
    if (fd->file_def->major_version == 0) {
	fd->file_def->major_version = CRAM_MAJOR_VERS(fd->version);
	fd->file_def->minor_version = CRAM_MINOR_VERS(fd->version);
	if (0 != cram_write_file_def(fd, fd->file_def))
	    return -1;
    }

    /* 1.0 requires and UNKNOWN read-group */
    if (CRAM_MAJOR_VERS(fd->version) == 1) {
	if (!sam_hdr_find_rg(hdr, "UNKNOWN"))
	    if (sam_hdr_add(hdr, "RG",
			    "ID", "UNKNOWN", "SM", "UNKNOWN", NULL))
		return -1;
    }

    /* Fix M5 strings */
    if (fd->refs && !fd->no_ref) {
	int i;
	for (i = 0; i < hdr->nref; i++) {
	    SAM_hdr_type *ty;
	    char *ref;

	    if (!(ty = sam_hdr_find(hdr, "SQ", "SN", hdr->ref[i].name)))
		return -1;

	    if (!sam_hdr_find_key(hdr, ty, "M5", NULL)) {
		char unsigned buf[16], buf2[33];
		int j, rlen;
		MD5_CTX md5;

		if (!fd->refs ||
		    !fd->refs->ref_id ||
		    !fd->refs->ref_id[i]) {
		    return -1;
		}
		rlen = fd->refs->ref_id[i]->length;
		MD5_Init(&md5);
		ref = cram_get_ref(fd, i, 1, rlen);
		if (NULL == ref) return -1;
		rlen = fd->refs->ref_id[i]->length; /* In case it just loaded */
		MD5_Update(&md5, ref, rlen);
		MD5_Final(buf, &md5);
		cram_ref_decr(fd->refs, i);

		for (j = 0; j < 16; j++) {
		    buf2[j*2+0] = "0123456789abcdef"[buf[j]>>4];
		    buf2[j*2+1] = "0123456789abcdef"[buf[j]&15];
		}
		buf2[32] = 0;
		if (sam_hdr_update(hdr, ty, "M5", buf2, NULL))
		    return -1;
	    }

	    if (fd->ref_fn) {
		char ref_fn[PATH_MAX];
		full_path(ref_fn, fd->ref_fn);
		if (sam_hdr_update(hdr, ty, "UR", ref_fn, NULL))
		    return -1;
	    }
	}
    }

    if (sam_hdr_rebuild(hdr))
	return -1;

    /* Length */
    header_len = sam_hdr_length(hdr);
    if (CRAM_MAJOR_VERS(fd->version) == 1) {
	if (-1 == int32_encode(fd, header_len))
	    return -1;

	/* Text data */
	if (header_len != hwrite(fd->fp, sam_hdr_str(hdr), header_len))
	    return -1;
    } else {
	/* Create block(s) inside a container */
	cram_block *b = cram_new_block(FILE_HEADER, 0);
	cram_container *c = cram_new_container(0, 0);
	int padded_length;
	char *pads;
	int is_cram_3 = (CRAM_MAJOR_VERS(fd->version) >= 3);

	if (!b || !c) {
	    if (b) cram_free_block(b);
	    if (c) cram_free_container(c);
	    return -1;
	}

	int32_put(b, header_len);
	BLOCK_APPEND(b, sam_hdr_str(hdr), header_len);
	BLOCK_UPLEN(b);

	// Compress header block if V3.0 and above
	if (CRAM_MAJOR_VERS(fd->version) >= 3 && fd->level > 0) {
	    int method = 1<<GZIP;
	    if (fd->use_bz2)
		method |= 1<<BZIP2;
	    if (fd->use_lzma)
		method |= 1<<LZMA;
	    cram_compress_block(fd, b, NULL, method, fd->level);
	} 

	if (blank_block) {
	    c->length = b->comp_size + 2 + 4*is_cram_3 +
		itf8_size(b->content_id)   + 
		itf8_size(b->uncomp_size)  +
		itf8_size(b->comp_size);

	    c->num_blocks = 2;
	    c->num_landmarks = 2;
	    if (!(c->landmark = (int32_t*)malloc(2*sizeof(*c->landmark)))) {
		cram_free_block(b);
		cram_free_container(c);
		return -1;
	    }
	    c->landmark[0] = 0;
	    c->landmark[1] = c->length;

	    // Plus extra storage for uncompressed secondary blank block
	    padded_length = MIN(c->length*.5, 10000);
	    c->length += padded_length + 2 + 4*is_cram_3 +
		itf8_size(b->content_id) + 
		itf8_size(padded_length)*2;
	} else {
	    // Pad the block instead.
	    c->num_blocks = 1;
	    c->num_landmarks = 1;
	    if (!(c->landmark = (int32_t*)malloc(sizeof(*c->landmark))))
		return -1;
	    c->landmark[0] = 0;

	    padded_length = MAX(c->length*1.5, 10000) - c->length;

	    c->length = b->comp_size + padded_length +
		2 + 4*is_cram_3 +
		itf8_size(b->content_id)   + 
		itf8_size(b->uncomp_size)  +
		itf8_size(b->comp_size);

		if (NULL == (pads = (char*)calloc(1, padded_length))) {
		cram_free_block(b);
		cram_free_container(c);
		return -1;
	    }
	    BLOCK_APPEND(b, pads, padded_length);
	    BLOCK_UPLEN(b);
	    free(pads);
	}

	if (-1 == cram_write_container(fd, c)) {
	    cram_free_block(b);
	    cram_free_container(c);
	    return -1;
	}

	if (-1 == cram_write_block(fd, b)) {
	    cram_free_block(b);
	    cram_free_container(c);
	    return -1;
	}

	if (blank_block) {
	    BLOCK_RESIZE(b, padded_length);
	    memset(BLOCK_DATA(b), 0, padded_length);
	    BLOCK_SIZE(b) = padded_length;
	    BLOCK_UPLEN(b);
	    b->method = RAW;
	    if (-1 == cram_write_block(fd, b)) {
		cram_free_block(b);
		cram_free_container(c);
		return -1;
	    }
	}

	cram_free_block(b);
	cram_free_container(c);
    }

    if (-1 == refs_from_header(fd->refs, fd, fd->header))
	return -1;
    if (-1 == refs2id(fd->refs, fd->header))
	return -1;

    if (0 != hflush(fd->fp))
	return -1;

    RP("=== Finishing saving header ===\n");

    return 0;
}

/* ----------------------------------------------------------------------
 * The top-level cram opening, closing and option handling
 */

/*
 * Initialises the lookup tables. These could be global statics, but they're
 * clumsy to setup in a multi-threaded environment unless we generate
 * verbatim code and include that.
 */
static void cram_init_tables(cram_fd *fd) {
    int i;

    memset(fd->L1, 4, 256);
    fd->L1['A'] = 0; fd->L1['a'] = 0;
    fd->L1['C'] = 1; fd->L1['c'] = 1;
    fd->L1['G'] = 2; fd->L1['g'] = 2;
    fd->L1['T'] = 3; fd->L1['t'] = 3;

    memset(fd->L2, 5, 256);
    fd->L2['A'] = 0; fd->L2['a'] = 0;
    fd->L2['C'] = 1; fd->L2['c'] = 1;
    fd->L2['G'] = 2; fd->L2['g'] = 2;
    fd->L2['T'] = 3; fd->L2['t'] = 3;
    fd->L2['N'] = 4; fd->L2['n'] = 4;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
	for (i = 0; i < 0x200; i++) {
	    int f = 0;

	    if (i & CRAM_FPAIRED)      f |= BAM_FPAIRED;
	    if (i & CRAM_FPROPER_PAIR) f |= BAM_FPROPER_PAIR;
	    if (i & CRAM_FUNMAP)       f |= BAM_FUNMAP;
	    if (i & CRAM_FREVERSE)     f |= BAM_FREVERSE;
	    if (i & CRAM_FREAD1)       f |= BAM_FREAD1;
	    if (i & CRAM_FREAD2)       f |= BAM_FREAD2;
	    if (i & CRAM_FSECONDARY)   f |= BAM_FSECONDARY;
	    if (i & CRAM_FQCFAIL)      f |= BAM_FQCFAIL;
	    if (i & CRAM_FDUP)         f |= BAM_FDUP;

	    fd->bam_flag_swap[i]  = f;
	}
    
	for (i = 0; i < 0x1000; i++) {
	    int g = 0;

	    if (i & BAM_FPAIRED)	   g |= CRAM_FPAIRED;
	    if (i & BAM_FPROPER_PAIR)  g |= CRAM_FPROPER_PAIR;
	    if (i & BAM_FUNMAP)        g |= CRAM_FUNMAP;
	    if (i & BAM_FREVERSE)      g |= CRAM_FREVERSE;
	    if (i & BAM_FREAD1)        g |= CRAM_FREAD1;
	    if (i & BAM_FREAD2)        g |= CRAM_FREAD2;
	    if (i & BAM_FSECONDARY)    g |= CRAM_FSECONDARY;
	    if (i & BAM_FQCFAIL)       g |= CRAM_FQCFAIL;
	    if (i & BAM_FDUP)          g |= CRAM_FDUP;

	    fd->cram_flag_swap[i] = g;
	}
    } else {
	/* NOP */
	for (i = 0; i < 0x1000; i++)
	    fd->bam_flag_swap[i] = i;
	for (i = 0; i < 0x1000; i++)
	    fd->cram_flag_swap[i] = i;
    }

    memset(fd->cram_sub_matrix, 4, 32*32);
    for (i = 0; i < 32; i++) {
	fd->cram_sub_matrix[i]['A'&0x1f]=0;
	fd->cram_sub_matrix[i]['C'&0x1f]=1;
	fd->cram_sub_matrix[i]['G'&0x1f]=2;
	fd->cram_sub_matrix[i]['T'&0x1f]=3;
	fd->cram_sub_matrix[i]['N'&0x1f]=4;
    }
    for (i = 0; i < 20; i+=4) {
	int j;
	for (j = 0; j < 20; j++) {
	    fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
	    fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
	    fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
	    fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][j]=3;
	}
	fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+0]&0x1f]=0;
	fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+1]&0x1f]=1;
	fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+2]&0x1f]=2;
	fd->cram_sub_matrix["ACGTN"[i>>2]&0x1f][CRAM_SUBST_MATRIX[i+3]&0x1f]=3;
    }
}

// Default version numbers for CRAM
static int major_version = 2;
static int minor_version = 1;

/*
 * Opens a CRAM file for read (mode "rb") or write ("wb").
 * The filename may be "-" to indicate stdin or stdout.
 *
 * Returns file handle on success
 *         NULL on failure.
 */
cram_fd *cram_open(const char *filename, const char *mode) {
    hFILE *fp;
    cram_fd *fd;
    char fmode[3]= { mode[0], '\0', '\0' };

    if (strlen(mode) > 1 && (mode[1] == 'b' || mode[1] == 'c')) {
	fmode[1] = 'b';
    }

    fp = hopen(filename, fmode);
    if (!fp)
	return NULL;

    fd = cram_dopen(fp, filename, mode);
    if (!fd)
	hclose_abruptly(fp);

    return fd;
}

/* Opens an existing stream for reading or writing.
 *
 * Returns file handle on success;
 *         NULL on failure.
 */
cram_fd *cram_dopen(hFILE *fp, const char *filename, const char *mode) {
    int i;
    char *cp;
	cram_fd *fd = (cram_fd*)calloc(1, sizeof(*fd));
    if (!fd)
	return NULL;

    fd->level = 5;
    for (i = 0; mode[i]; i++) {
	if (mode[i] >= '0' && mode[i] <= '9') {
	    fd->level = mode[i] - '0';
	    break;
	}
    }

    fd->fp = fp;
    fd->mode = *mode;
    fd->first_container = 0;

    if (fd->mode == 'r') {
	/* Reader */

	if (!(fd->file_def = cram_read_file_def(fd)))
	    goto err;

	fd->version = fd->file_def->major_version * 256 +
	    fd->file_def->minor_version;

	if (!(fd->header = cram_read_SAM_hdr(fd)))
	    goto err;

    } else {
	/* Writer */
		cram_file_def *def = (cram_file_def*)calloc(1, sizeof(*def));
	if (!def)
	    return NULL;

	fd->file_def = def;

	def->magic[0] = 'C';
	def->magic[1] = 'R';
	def->magic[2] = 'A';
	def->magic[3] = 'M';
	def->major_version = 0; // Indicator to write file def later.
	def->minor_version = 0;
	memset(def->file_id, 0, 20);
	strncpy(def->file_id, filename, 20);

	fd->version = major_version * 256 + minor_version;

	/* SAM header written later along with this file_def */
    }

    cram_init_tables(fd);

    fd->prefix = strdup((cp = (char*)strrchr(filename, '/')) ? cp+1 : filename);
    if (!fd->prefix)
	goto err;
    fd->first_base = fd->last_base = -1;
    fd->record_counter = 0;

    fd->ctr = NULL;
    fd->refs  = refs_create();
    if (!fd->refs)
	goto err;
    fd->ref_id = -2;
    fd->ref = NULL;

    fd->decode_md = 0;
    fd->verbose = 0;
    fd->seqs_per_slice = SEQS_PER_SLICE;
    fd->slices_per_container = SLICE_PER_CNT;
    fd->embed_ref = 0;
    fd->no_ref = 0;
    fd->ignore_md5 = 0;
    fd->use_bz2 = 0;
    fd->use_rans = (CRAM_MAJOR_VERS(fd->version) >= 3);
    fd->use_lzma = 0;
    fd->multi_seq = -1;
    fd->unsorted   = 0;
    fd->shared_ref = 0;

    fd->index       = NULL;
    fd->own_pool    = 0;
    fd->pool        = NULL;
    fd->rqueue      = NULL;
    fd->job_pending = NULL;
    fd->ooc         = 0;
    fd->required_fields = INT_MAX;

    for (i = 0; i < DS_END; i++)
	fd->m[i] = cram_new_metrics();

    fd->range.refid = -2; // no ref.
    fd->eof = 1;          // See samtools issue #150
    fd->ref_fn = NULL;

    fd->bl = NULL;

    /* Initialise dummy refs from the @SQ headers */
    if (-1 == refs_from_header(fd->refs, fd, fd->header))
	goto err;

    return fd;

 err:
    if (fd)
	free(fd);

    return NULL;
}

/*
 * Seek within a CRAM file.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek(cram_fd *fd, off_t offset, int whence) {
    char buf[65536];

    fd->ooc = 0;

    if (hseek(fd->fp, offset, whence) >= 0)
	return 0;

    if (!(whence == SEEK_CUR && offset >= 0))
	return -1;

    /* Couldn't fseek, but we're in SEEK_CUR mode so read instead */
    while (offset > 0) {
	int len = MIN(65536, offset);
	if (len != hread(fd->fp, buf, len))
	    return -1;
	offset -= len;
    }

    return 0;
}

/*
 * Flushes a CRAM file.
 * Useful for when writing to stdout without wishing to close the stream.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_flush(cram_fd *fd) {
    if (!fd)
	return -1;

    if (fd->mode == 'w' && fd->ctr) {
	if(fd->ctr->slice)
	    fd->ctr->curr_slice++;
	if (-1 == cram_flush_container_mt(fd, fd->ctr))
	    return -1;
    }

    return 0;
}

/*
 * Closes a CRAM file.
 * Returns 0 on success
 *        -1 on failure
 */
int cram_close(cram_fd *fd) {
    spare_bams *bl, *next;
    int i;
	
    if (!fd)
	return -1;

    if (fd->mode == 'w' && fd->ctr) {
	if(fd->ctr->slice)
	    fd->ctr->curr_slice++;
	if (-1 == cram_flush_container_mt(fd, fd->ctr))
	    return -1;
    }

    if (fd->pool) {
	t_pool_flush(fd->pool);

	if (0 != cram_flush_result(fd))
	    return -1;

#if !defined(_WIN32) && defined(USE_PTHREAD)
	pthread_mutex_destroy(&fd->metrics_lock);
	pthread_mutex_destroy(&fd->ref_lock);
	pthread_mutex_destroy(&fd->bam_list_lock);
#endif

	fd->ctr = NULL; // prevent double freeing

	//fprintf(stderr, "CRAM: destroy queue %p\n", fd->rqueue);

	t_results_queue_destroy(fd->rqueue);
    }

    if (fd->mode == 'w') {
	/* Write EOF block */
	if (CRAM_MAJOR_VERS(fd->version) == 3) {
	    if (38 != hwrite(fd->fp,
			     "\x0f\x00\x00\x00\xff\xff\xff\xff" // Cont HDR
			     "\x0f\xe0\x45\x4f\x46\x00\x00\x00" // Cont HDR
			     "\x00\x01\x00"                     // Cont HDR
			     "\x05\xbd\xd9\x4f"                 // CRC32
			     "\x00\x01\x00\x06\x06"             // Comp.HDR blk
			     "\x01\x00\x01\x00\x01\x00"         // Comp.HDR blk
			     "\xee\x63\x01\x4b",                // CRC32
			     38))
		return -1;
	} else {
	    if (30 != hwrite(fd->fp,
			     "\x0b\x00\x00\x00\xff\xff\xff\xff"
			     "\x0f\xe0\x45\x4f\x46\x00\x00\x00"
			     "\x00\x01\x00\x00\x01\x00\x06\x06"
			     "\x01\x00\x01\x00\x01\x00", 30))
		return -1;
	}
    }

    for (bl = fd->bl; bl; bl = next) {
	int i, max_rec = fd->seqs_per_slice * fd->slices_per_container;

	next = bl->next;
	for (i = 0; i < max_rec; i++) {
	    if (bl->bams[i])
		bam_free(bl->bams[i]);
	}
	free(bl->bams);
	free(bl);
    }

    if (hclose(fd->fp) != 0)
	return -1;

    if (fd->file_def)
	cram_free_file_def(fd->file_def);

    if (fd->header)
	sam_hdr_free(fd->header);

    free(fd->prefix);

    if (fd->ctr)
	cram_free_container(fd->ctr);

    if (fd->refs)
	refs_free(fd->refs);
    if (fd->ref_free)
        free(fd->ref_free);

    for (i = 0; i < DS_END; i++)
	if (fd->m[i])
	    free(fd->m[i]);

    if (fd->index)
	cram_index_free(fd);

    if (fd->own_pool && fd->pool)
	t_pool_destroy(fd->pool, 0);

    free(fd);
    return 0;
}

/*
 * Returns 1 if we hit an EOF while reading.
 */
int cram_eof(cram_fd *fd) {
    return fd->eof;
}


/* 
 * Sets options on the cram_fd. See CRAM_OPT_* definitions in cram_structs.h.
 * Use this immediately after opening.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_set_option(cram_fd *fd, enum cram_option opt, ...) {
    int r;
    va_list args;

    va_start(args, opt);
    r = cram_set_voption(fd, opt, args);
    va_end(args);

    return r;
}

/*
 * Sets options on the cram_fd. See CRAM_OPT_* definitions in cram_structs.h.
 * Use this immediately after opening.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_set_voption(cram_fd *fd, enum cram_option opt, va_list args) {
    refs_t *refs;

    if (!fd)
	return -1;

    switch (opt) {
    case CRAM_OPT_DECODE_MD:
	fd->decode_md = va_arg(args, int);
	break;

    case CRAM_OPT_PREFIX:
	if (fd->prefix)
	    free(fd->prefix);
	if (!(fd->prefix = strdup(va_arg(args, char *))))
	    return -1;
	break;

    case CRAM_OPT_VERBOSITY:
	fd->verbose = va_arg(args, int);
	break;

    case CRAM_OPT_SEQS_PER_SLICE:
	fd->seqs_per_slice = va_arg(args, int);
	break;

    case CRAM_OPT_SLICES_PER_CONTAINER:
	fd->slices_per_container = va_arg(args, int);
	break;

    case CRAM_OPT_EMBED_REF:
	fd->embed_ref = va_arg(args, int);
	break;

    case CRAM_OPT_NO_REF:
	fd->no_ref = va_arg(args, int);
	break;

    case CRAM_OPT_IGNORE_MD5:
	fd->ignore_md5 = va_arg(args, int);
	break;

    case CRAM_OPT_USE_BZIP2:
	fd->use_bz2 = va_arg(args, int);
	break;

    case CRAM_OPT_USE_RANS:
	fd->use_rans = va_arg(args, int);
	break;

    case CRAM_OPT_USE_LZMA:
	fd->use_lzma = va_arg(args, int);
	break;

    case CRAM_OPT_SHARED_REF:
	fd->shared_ref = 1;
	refs = va_arg(args, refs_t *);
	if (refs != fd->refs) {
	    if (fd->refs)
		refs_free(fd->refs);
	    fd->refs = refs;
	    fd->refs->count++;
	}
	break;

    case CRAM_OPT_RANGE:
	fd->range = *va_arg(args, cram_range *);
	return cram_seek_to_refpos(fd, &fd->range);

    case CRAM_OPT_REFERENCE:
	return cram_load_reference(fd, va_arg(args, char *));

    case CRAM_OPT_VERSION: {
	int major, minor;
	char *s = va_arg(args, char *);
	if (2 != sscanf(s, "%d.%d", &major, &minor)) {
	    fprintf(stderr, "Malformed version string %s\n", s);
	    return -1;
	}
	if (!((major == 1 &&  minor == 0) ||
	      (major == 2 && (minor == 0 || minor == 1)) ||
	      (major == 3 &&  minor == 0))) {
	    fprintf(stderr, "Unknown version string; "
		    "use 1.0, 2.0, 2.1 or 3.0\n");
	    return -1;
	}
	fd->version = major*256 + minor;

	if (CRAM_MAJOR_VERS(fd->version) >= 3)
	    fd->use_rans = 1;
	break;
    }

    case CRAM_OPT_MULTI_SEQ_PER_SLICE:
	fd->multi_seq = va_arg(args, int);
	break;

    case CRAM_OPT_NTHREADS: {
	int nthreads =  va_arg(args, int);
        if (nthreads > 1) {
            if (!(fd->pool = t_pool_init(nthreads*2, nthreads)))
                return -1;

	    fd->rqueue = t_results_queue_init();

#if !defined(_WIN32) && defined(USE_PTHREAD)
	    pthread_mutex_init(&fd->metrics_lock, NULL);
	    pthread_mutex_init(&fd->ref_lock, NULL);
	    pthread_mutex_init(&fd->bam_list_lock, NULL);
#endif
	    fd->shared_ref = 1;
	    fd->own_pool = 1;
        }
	break;
    }

    case CRAM_OPT_THREAD_POOL:
	fd->pool = va_arg(args, t_pool *);
	if (fd->pool) {
	    fd->rqueue = t_results_queue_init();

#if !defined(_WIN32) && defined(USE_PTHREAD)
	    pthread_mutex_init(&fd->metrics_lock, NULL);
	    pthread_mutex_init(&fd->ref_lock, NULL);
	    pthread_mutex_init(&fd->bam_list_lock, NULL);
#endif
	}
	fd->shared_ref = 1; // Needed to avoid clobbering ref between threads
	fd->own_pool = 0;

	//fd->qsize = 1;
	//fd->decoded = calloc(fd->qsize, sizeof(cram_container *));
	//t_pool_dispatch(fd->pool, cram_decoder_thread, fd);
	break;

    case CRAM_OPT_REQUIRED_FIELDS:
	fd->required_fields = va_arg(args, int);
	break;

    default:
	fprintf(stderr, "Unknown CRAM option code %d\n", opt);
	return -1;
    }

    return 0;
}
