/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>

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
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

/* The BGZF library was originally written by Bob Handsaker from the Broad
 * Institute. It was later improved by the SAMtools developers. */

#ifndef __BGZF_H
#define __BGZF_H

#include <stdint.h>
#include <stdio.h>
#include <zlib.h>
#include <sys/types.h>

#define BGZF_BLOCK_SIZE     0xff00 // make sure compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE
#define BGZF_MAX_BLOCK_SIZE 0x10000

#define BGZF_ERR_ZLIB   1
#define BGZF_ERR_HEADER 2
#define BGZF_ERR_IO     4
#define BGZF_ERR_MISUSE 8

struct hFILE;
struct bgzf_mtaux_t;
typedef struct __bgzidx_t bgzidx_t;

struct BGZF {
	int errcode:16, is_write:2, is_be:2, compress_level:9, is_compressed:2, is_gzip:1;
	int cache_size;
    int block_length, block_offset;
    int64_t block_address, uncompressed_address;
    void *uncompressed_block, *compressed_block;
	void *cache; // a pointer to a hash table
    struct hFILE *fp; // actual file handle
    struct bgzf_mtaux_t *mt; // only used for multi-threading
    bgzidx_t *idx;      // BGZF index
    int idx_build_otf;  // build index on the fly, set by bgzf_index_build_init()
    z_stream *gz_stream;// for gzip-compressed files
};
#ifndef HTS_BGZF_TYPEDEF
typedef struct BGZF BGZF;
#define HTS_BGZF_TYPEDEF
#endif

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

	/******************
	 * Basic routines *
	 ******************/

	/**
	 * Open an existing file descriptor for reading or writing.
	 *
	 * @param fd    file descriptor
	 * @param mode  mode matching /[rwa][u0-9]+/: 'r' for reading, 'w' for
	 *              writing, or 'a' for appending, while a digit specifies
	 *              the zlib compression level.
     *              Note that there is a distinction between 'u' and '0': the
     *              first yields plain uncompressed output whereas the latter
     *              outputs uncompressed data wrapped in the zlib format.
	 * @return      BGZF file handler; 0 on error
	 */
	BGZF* bgzf_dopen(int fd, const char *mode);

	#define bgzf_fdopen(fd, mode) bgzf_dopen((fd), (mode)) // for backward compatibility

	/**
	 * Open the specified file for reading or writing.
	 */
	BGZF* bgzf_open(const char* path, const char *mode);

	/**
	 * Open an existing hFILE stream for reading or writing.
	 */
	BGZF* bgzf_hopen(struct hFILE *fp, const char *mode);

	/**
	 * Close the BGZF and free all associated resources.
	 *
	 * @param fp    BGZF file handler
	 * @return      0 on success and -1 on error
	 */
	int bgzf_close(BGZF *fp);

	/**
	 * Read up to _length_ bytes from the file storing into _data_.
	 *
	 * @param fp     BGZF file handler
	 * @param data   data array to read into
	 * @param length size of data to read
	 * @return       number of bytes actually read; 0 on end-of-file and -1 on error
	 */
	ssize_t bgzf_read(BGZF *fp, void *data, size_t length);

	/**
	 * Write _length_ bytes from _data_ to the file.  If no I/O errors occur,
	 * the complete _length_ bytes will be written (or queued for writing).
	 *
	 * @param fp     BGZF file handler
	 * @param data   data array to write
	 * @param length size of data to write
	 * @return       number of bytes written (i.e., _length_); negative on error
	 */
	ssize_t bgzf_write(BGZF *fp, const void *data, size_t length);

	/**
	 * Read up to _length_ bytes directly from the underlying stream without
	 * decompressing.  Bypasses BGZF blocking, so must be used with care in
	 * specialised circumstances only.
	 *
	 * @param fp     BGZF file handler
	 * @param data   data array to read into
	 * @param length number of raw bytes to read
	 * @return       number of bytes actually read; 0 on end-of-file and -1 on error
	 */
	ssize_t bgzf_raw_read(BGZF *fp, void *data, size_t length);

	/**
	 * Write _length_ bytes directly to the underlying stream without
	 * compressing.  Bypasses BGZF blocking, so must be used with care
	 * in specialised circumstances only.
	 *
	 * @param fp     BGZF file handler
	 * @param data   data array to write
	 * @param length number of raw bytes to write
	 * @return       number of bytes actually written; -1 on error
	 */
	ssize_t bgzf_raw_write(BGZF *fp, const void *data, size_t length);

	/**
	 * Write the data in the buffer to the file.
	 */
	int bgzf_flush(BGZF *fp);

	/**
	 * Return a virtual file pointer to the current location in the file.
	 * No interpetation of the value should be made, other than a subsequent
	 * call to bgzf_seek can be used to position the file at the same point.
	 * Return value is non-negative on success.
	 */
	#define bgzf_tell(fp) (((fp)->block_address << 16) | ((fp)->block_offset & 0xFFFF))

	/**
	 * Set the file to read from the location specified by _pos_.
	 *
	 * @param fp     BGZF file handler
	 * @param pos    virtual file offset returned by bgzf_tell()
	 * @param whence must be SEEK_SET
	 * @return       0 on success and -1 on error
	 */
	int64_t bgzf_seek(BGZF *fp, int64_t pos, int whence);

	/**
	 * Check if the BGZF end-of-file (EOF) marker is present
	 *
	 * @param fp    BGZF file handler opened for reading
	 * @return      1 if the EOF marker is present and correct;
	 *              2 if it can't be checked, e.g., because fp isn't seekable;
	 *              0 if the EOF marker is absent;
	 *              -1 (with errno set) on error
	 */
	int bgzf_check_EOF(BGZF *fp);

	/**
	 * Check if a file is in the BGZF format
	 *
	 * @param fn    file name
	 * @return      1 if _fn_ is BGZF; 0 if not or on I/O error
	 */
	 int bgzf_is_bgzf(const char *fn);

	/*********************
	 * Advanced routines *
	 *********************/

	/**
	 * Set the cache size. Only effective when compiled with -DBGZF_CACHE.
	 *
	 * @param fp    BGZF file handler
	 * @param size  size of cache in bytes; 0 to disable caching (default)
	 */
	void bgzf_set_cache_size(BGZF *fp, int size);

	/**
	 * Flush the file if the remaining buffer size is smaller than _size_ 
	 * @return      0 if flushing succeeded or was not needed; negative on error
	 */
	int bgzf_flush_try(BGZF *fp, ssize_t size);

	/**
	 * Read one byte from a BGZF file. It is faster than bgzf_read()
	 * @param fp     BGZF file handler
	 * @return       byte read; -1 on end-of-file or error
	 */
	int bgzf_getc(BGZF *fp);

	/**
	 * Read one line from a BGZF file. It is faster than bgzf_getc()
	 *
	 * @param fp     BGZF file handler
	 * @param delim  delimitor
	 * @param str    string to write to; must be initialized
	 * @return       length of the string; 0 on end-of-file; negative on error
	 */
	int bgzf_getline(BGZF *fp, int delim, kstring_t *str);

	/**
	 * Read the next BGZF block.
	 */
	int bgzf_read_block(BGZF *fp);

	/**
	 * Enable multi-threading (only effective on writing and when the
	 * library was compiled with -DBGZF_MT)
	 *
	 * @param fp          BGZF file handler; must be opened for writing
	 * @param n_threads   #threads used for writing
	 * @param n_sub_blks  #blocks processed by each thread; a value 64-256 is recommended
	 */
	int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks);


	/*******************
	 * bgzidx routines *
	 *******************/

	/**
     *  Position BGZF at the uncompressed offset 
     *
     *  @param fp           BGZF file handler; must be opened for reading
     *  @param uoffset      file offset in the uncompressed data
     *  @param where        SEEK_SET supported atm
     *
     *  Returns 0 on success and -1 on error.
	 */
    int bgzf_useek(BGZF *fp, long uoffset, int where);

	/**
     *  Position in uncompressed BGZF
     *
     *  @param fp           BGZF file handler; must be opened for reading
     *
     *  Returns the current offset on success and -1 on error.
	 */
    long bgzf_utell(BGZF *fp);

	/**
	 * Tell BGZF to build index while compressing.
     *
	 * @param fp          BGZF file handler; can be opened for reading or writing.
     *
     * Returns 0 on success and -1 on error.
	 */
	int bgzf_index_build_init(BGZF *fp);

   	/**
	 * Load BGZF index
	 *
	 * @param fp          BGZF file handler
	 * @param bname       base name
     * @param suffix      suffix to add to bname (can be NULL)
     *
     * Returns 0 on success and -1 on error.
	 */
    int bgzf_index_load(BGZF *fp, const char *bname, const char *suffix);

   	/**
	 * Save BGZF index
	 *
	 * @param fp          BGZF file handler
	 * @param bname       base name
     * @param suffix      suffix to add to bname (can be NULL)
     *
     * Returns 0 on success and -1 on error.
	 */
    int bgzf_index_dump(BGZF *fp, const char *bname, const char *suffix);

#ifdef __cplusplus
}
#endif

#endif
