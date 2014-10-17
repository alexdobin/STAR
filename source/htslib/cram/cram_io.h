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
 * Implements the low level CRAM I/O primitives.
 * This includes basic data types such as byte, int, ITF-8,
 * maps, bitwise I/O, etc.
 */

#ifndef _CRAM_IO_H_
#define _CRAM_IO_H_

#ifdef __cplusplus
extern "C" {
#endif

#define ITF8_MACROS

#include <stdint.h>
#include <cram/misc.h>

/**@{ ----------------------------------------------------------------------
 * ITF8 encoding and decoding.
 *
 * Also see the itf8_get and itf8_put macros.
 */

/*! INTERNAL: Converts two characters into an integer for use in switch{} */
#define CRAM_KEY(a,b) (((a)<<8)|((b)))

/*! Reads an integer in ITF-8 encoding from 'fd' and stores it in
 * *val.
 *
 * @return
 * Returns the number of bytes read on success;
 *        -1 on failure
 */
int itf8_decode(cram_fd *fd, int32_t *val);

#ifndef ITF8_MACROS
/*! Reads an integer in ITF-8 encoding from 'cp' and stores it in
 * *val.
 *
 * @return
 * Returns the number of bytes read on success;
 *        -1 on failure
 */
int itf8_get(char *cp, int32_t *val_p);

/*! Stores a value to memory in ITF-8 format.
 *
 * @return
 * Returns the number of bytes required to store the number.
 * This is a maximum of 5 bytes.
 */
int itf8_put(char *cp, int32_t val);

#else

/*
 * Macro implementations of the above
 */
#define itf8_get(c,v) (((uc)(c)[0]<0x80)?(*(v)=(uc)(c)[0],1):(((uc)(c)[0]<0xc0)?(*(v)=(((uc)(c)[0]<<8)|(uc)(c)[1])&0x3fff,2):(((uc)(c)[0]<0xe0)?(*(v)=(((uc)(c)[0]<<16)|((uc)(c)[1]<<8)|(uc)(c)[2])&0x1fffff,3):(((uc)(c)[0]<0xf0)?(*(v)=(((uc)(c)[0]<<24)|((uc)(c)[1]<<16)|((uc)(c)[2]<<8)|(uc)(c)[3])&0x0fffffff,4):(*(v)=(((uc)(c)[0]&0x0f)<<28)|((uc)(c)[1]<<20)|((uc)(c)[2]<<12)|((uc)(c)[3]<<4)|((uc)(c)[4]&0x0f),5)))))

#define itf8_put(c,v) ((!((v)&~0x7f))?((c)[0]=(v),1):(!((v)&~0x3fff))?((c)[0]=((v)>>8)|0x80,(c)[1]=(v)&0xff,2):(!((v)&~0x1fffff))?((c)[0]=((v)>>16)|0xc0,(c)[1]=((v)>>8)&0xff,(c)[2]=(v)&0xff,3):(!((v)&~0xfffffff))?((c)[0]=((v)>>24)|0xe0,(c)[1]=((v)>>16)&0xff,(c)[2]=((v)>>8)&0xff,(c)[3]=(v)&0xff,4):((c)[0]=0xf0|(((v)>>28)&0xff),(c)[1]=((v)>>20)&0xff,(c)[2]=((v)>>12)&0xff,(c)[3]=((v)>>4)&0xff,(c)[4]=(v)&0xf,5))

#define itf8_size(v) ((!((v)&~0x7f))?1:(!((v)&~0x3fff))?2:(!((v)&~0x1fffff))?3:(!((v)&~0xfffffff))?4:5)

#endif

/*! Pushes a value in ITF8 format onto the end of a block.
 *
 * This shouldn't be used for high-volume data as it is not the fastest
 * method.
 *
 * @return
 * Returns the number of bytes written
 */
int itf8_put_blk(cram_block *blk, int val);

/**@}*/
/**@{ ----------------------------------------------------------------------
 * CRAM blocks - the dynamically growable data block. We have code to
 * create, update, (un)compress and read/write.
 *
 * These are derived from the deflate_interlaced.c blocks, but with the
 * CRAM extension of content types and IDs.
 */

/*! Allocates a new cram_block structure with a specified content_type and
 * id.
 *
 * @return
 * Returns block pointer on success;
 *         NULL on failure
 */
cram_block *cram_new_block(enum cram_content_type content_type,
			   int content_id);

/*! Reads a block from a cram file.
 *
 * @return
 * Returns cram_block pointer on success;
 *         NULL on failure
 */
cram_block *cram_read_block(cram_fd *fd);

/*! Writes a CRAM block.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_write_block(cram_fd *fd, cram_block *b);

/*! Frees a CRAM block, deallocating internal data too.
 */
void cram_free_block(cram_block *b);

/*! Uncompress a memory block using Zlib.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
char *zlib_mem_inflate(char *cdata, size_t csize, size_t *size);

/*! Uncompresses a CRAM block, if compressed.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_uncompress_block(cram_block *b);

/*! Compresses a block.
 *
 * Compresses a block using one of two different zlib strategies. If we only
 * want one choice set strat2 to be -1.
 *
 * The logic here is that sometimes Z_RLE does a better job than Z_FILTERED
 * or Z_DEFAULT_STRATEGY on quality data. If so, we'd rather use it as it is
 * significantly faster.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_compress_block(cram_fd *fd, cram_block *b, cram_metrics *metrics,
			int level,  int strat,
			int level2, int strat2);

cram_metrics *cram_new_metrics(void);
char *cram_block_method2str(enum cram_block_method m);
char *cram_content_type2str(enum cram_content_type t);

/* --- Accessor macros for manipulating blocks on a byte by byte basis --- */

/* Block size and data pointer. */
#define BLOCK_SIZE(b) ((b)->byte)
#define BLOCK_DATA(b) ((b)->data)

/* Returns the address one past the end of the block */
#define BLOCK_END(b) (&(b)->data[(b)->byte])

/* Request block to be at least 'l' bytes long */
#define BLOCK_RESIZE(b,l)					\
    do {							\
	while((b)->alloc <= (l)) {				\
	    (b)->alloc = (b)->alloc ? (b)->alloc*1.5 : 1024;	\
	    (b)->data = realloc((b)->data, (b)->alloc);		\
	}							\
     } while(0)

/* Ensure the block can hold at least another 'l' bytes */
#define BLOCK_GROW(b,l) BLOCK_RESIZE((b), BLOCK_SIZE((b)) + (l))

/* Append string 's' of length 'l' */
#define BLOCK_APPEND(b,s,l)		  \
    do {				  \
        BLOCK_GROW((b),(l));		  \
        memcpy(BLOCK_END((b)), (s), (l)); \
	BLOCK_SIZE((b)) += (l);		  \
    } while (0)

/* Append as single character 'c' */
#define BLOCK_APPEND_CHAR(b,c)		  \
    do {				  \
        BLOCK_GROW((b),1);		  \
	(b)->data[(b)->byte++] = (c);	  \
    } while (0)

/* Append via sprintf with 1 arg */
#define BLOCK_APPENDF_1(b,buf,fmt, a1)			\
    do {						\
	int l = sprintf((buf), (fmt), (a1));		\
	BLOCK_APPEND((b), (buf), l);			\
    } while (0)

/* Append via sprintf with 2 args */
#define BLOCK_APPENDF_2(b,buf,fmt, a1,a2)		\
    do {						\
	int l = sprintf((buf), (fmt), (a1), (a2));	\
	BLOCK_APPEND((b), (buf), l);			\
    } while (0)

#define BLOCK_UPLEN(b) \
    (b)->comp_size = (b)->uncomp_size = BLOCK_SIZE((b))

/**@}*/
/**@{ ----------------------------------------------------------------------
 * Reference sequence handling
 */

/*! Loads a reference set from fn and stores in the cram_fd.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_load_reference(cram_fd *fd, char *fn);

/*! Generates a lookup table in refs based on the SQ headers in SAM_hdr.
 *
 * Indexes references by the order they appear in a BAM file. This may not
 * necessarily be the same order they appear in the fasta reference file.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int refs2id(refs_t *r, SAM_hdr *bfd);

void refs_free(refs_t *r);

/*! Returns a portion of a reference sequence from start to end inclusive.
 *
 * The returned pointer is owned by the cram_file fd and should not be freed
 * by the caller. It is valid only until the next cram_get_ref is called
 * with the same fd parameter (so is thread-safe if given multiple files).
 *
 * To return the entire reference sequence, specify start as 1 and end
 * as 0.
 *
 * @return
 * Returns reference on success;
 *         NULL on failure
 */
char *cram_get_ref(cram_fd *fd, int id, int start, int end);
void cram_ref_incr(refs_t *r, int id);
void cram_ref_decr(refs_t *r, int id);
/**@}*/
/**@{ ----------------------------------------------------------------------
 * Containers
 */

/*! Creates a new container, specifying the maximum number of slices
 * and records permitted.
 *
 * @return
 * Returns cram_container ptr on success;
 *         NULL on failure
 */
cram_container *cram_new_container(int nrec, int nslice);
void cram_free_container(cram_container *c);

/*! Reads a container header.
 *
 * @return
 * Returns cram_container on success;
 *         NULL on failure or no container left (fd->err == 0).
 */
cram_container *cram_read_container(cram_fd *fd);

/*! Writes a container structure.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_write_container(cram_fd *fd, cram_container *h);

/*! Flushes a container to disk.
 *
 * Flushes a completely or partially full container to disk, writing
 * container structure, header and blocks. This also calls the encoder
 * functions.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_flush_container(cram_fd *fd, cram_container *c);
int cram_flush_container_mt(cram_fd *fd, cram_container *c);


/**@}*/
/**@{ ----------------------------------------------------------------------
 * Compression headers; the first part of the container
 */

/*! Creates a new blank container compression header
 *
 * @return
 * Returns header ptr on success;
 *         NULL on failure
 */
cram_block_compression_hdr *cram_new_compression_header(void);

/*! Frees a cram_block_compression_hdr */
void cram_free_compression_header(cram_block_compression_hdr *hdr);


/**@}*/
/**@{ ----------------------------------------------------------------------
 * Slices and slice headers
 */

/*! Frees a slice header */
void cram_free_slice_header(cram_block_slice_hdr *hdr);

/*! Frees a slice */
void cram_free_slice(cram_slice *s);

/*! Creates a new empty slice in memory, for subsequent writing to
 * disk.
 *
 * @return
 * Returns cram_slice ptr on success;
 *         NULL on failure
 */
cram_slice *cram_new_slice(enum cram_content_type type, int nrecs);

/*! Loads an entire slice.
 *
 * FIXME: In 1.0 the native unit of slices within CRAM is broken
 * as slices contain references to objects in other slices.
 * To work around this while keeping the slice oriented outer loop
 * we read all slices and stitch them together into a fake large
 * slice instead.
 *
 * @return
 * Returns cram_slice ptr on success;
 *         NULL on failure
 */
cram_slice *cram_read_slice(cram_fd *fd);



/**@}*/
/**@{ ----------------------------------------------------------------------
 * CRAM file definition (header)
 */

/*! Reads a CRAM file definition structure.
 *
 * @return
 * Returns file_def ptr on success;
 *         NULL on failure
 */
cram_file_def *cram_read_file_def(cram_fd *fd);

/*! Writes a cram_file_def structure to cram_fd.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_write_file_def(cram_fd *fd, cram_file_def *def);

/*! Frees a cram_file_def structure. */
void cram_free_file_def(cram_file_def *def);


/**@}*/
/**@{ ----------------------------------------------------------------------
 * SAM header I/O
 */

/*! Reads the SAM header from the first CRAM data block.
 *
 * Also performs minimal parsing to extract read-group
 * and sample information.
 *
 * @return
 * Returns SAM hdr ptr on success;
 *         NULL on failure
 */
SAM_hdr *cram_read_SAM_hdr(cram_fd *fd);

/*! Writes a CRAM SAM header.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_write_SAM_hdr(cram_fd *fd, SAM_hdr *hdr);


/**@}*/
/**@{ ----------------------------------------------------------------------
 * The top-level cram opening, closing and option handling
 */

/*! Opens a CRAM file for read (mode "rb") or write ("wb").
 *
 * The filename may be "-" to indicate stdin or stdout.
 *
 * @return
 * Returns file handle on success;
 *         NULL on failure.
 */
cram_fd *cram_open(const char *filename, const char *mode);

/*! Opens an existing stream for reading or writing.
 *
 * @return
 * Returns file handle on success;
 *         NULL on failure.
 *
 * cram_FILE is either htslib's hFILE or stdio's FILE, depending on how
 * cram_structs.h has been configured.
 */
cram_fd *cram_dopen(cram_FILE *fp, const char *filename, const char *mode);

/*! Closes a CRAM file.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_close(cram_fd *fd);

/*
 * Seek within a CRAM file.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek(cram_fd *fd, off_t offset, int whence);

/*
 * Flushes a CRAM file.
 * Useful for when writing to stdout without wishing to close the stream.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_flush(cram_fd *fd);

/*! Checks for end of file on a cram_fd stream.
 *
 * @return
 * Returns 0 if not at end of file
 *         1 if we hit an expected EOF (end of range or EOF block)
 *         2 for other EOF (end of stream without EOF block)
 */
int cram_eof(cram_fd *fd);

/*! Sets options on the cram_fd.
 *
 * See CRAM_OPT_* definitions in cram_structs.h.
 * Use this immediately after opening.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_set_option(cram_fd *fd, enum cram_option opt, ...);

/*! Sets options on the cram_fd.
 *
 * See CRAM_OPT_* definitions in cram_structs.h.
 * Use this immediately after opening.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_set_voption(cram_fd *fd, enum cram_option opt, va_list args);

/*!
 * Attaches a header to a cram_fd.
 *
 * This should be used when creating a new cram_fd for writing where
 * we have an SAM_hdr already constructed (eg from a file we've read
 * in).
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_set_header(cram_fd *fd, SAM_hdr *hdr);


#ifdef __cplusplus
}
#endif

#endif /* _CRAM_IO_H_ */
