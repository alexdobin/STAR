/*  hts.h -- format-neutral I/O, indexing, and iterator API functions.

    Copyright (C) 2012-2014 Genome Research Ltd.
    Copyright (C) 2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef HTSLIB_HTS_H
#define HTSLIB_HTS_H

#include <stddef.h>
#include <stdint.h>
#include "common\cross_platform.h"

#ifndef HTS_BGZF_TYPEDEF
typedef struct BGZF BGZF;
#define HTS_BGZF_TYPEDEF
#endif
struct cram_fd;
struct hFILE;

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
    size_t l, m;
    char *s;
} kstring_t;
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/**
 * hts_expand()  - expands memory block pointed to by $ptr;
 * hts_expand0()   the latter sets the newly allocated part to 0.
 *
 * @param n     requested number of elements of type type_t
 * @param m     size of memory allocated
 */
#define hts_expand(type_t, n, m, ptr) if ((n) > (m)) { \
        (m) = (n); kroundup32(m); \
        (ptr) = (type_t*)realloc((ptr), (m) * sizeof(type_t)); \
    }
#define hts_expand0(type_t, n, m, ptr) if ((n) > (m)) { \
        int t = (m); (m) = (n); kroundup32(m); \
        (ptr) = (type_t*)realloc((ptr), (m) * sizeof(type_t)); \
        memset(((type_t*)ptr)+t,0,sizeof(type_t)*((m)-t)); \
    }

/************
 * File I/O *
 ************/

// Add new entries only at the end (but before the *_maximum entry)
// of these enums, as their numbering is part of the htslib ABI.

enum htsFormatCategory {
    unknown_category,
    sequence_data,    // Sequence data -- SAM, BAM, CRAM, etc
    variant_data,     // Variant calling data -- VCF, BCF, etc
    index_file,       // Index file associated with some data file
    region_list,      // Coordinate intervals or regions -- BED, etc
    category_maximum = 32767
};

enum htsExactFormat {
    unknown_format,
    binary_format, text_format,
    sam, bam, bai, cram, crai, vcf, bcf, csi, gzi, tbi, bed,
    format_maximum = 32767
};

enum htsCompression {
    no_compression, gzip, bgzf, custom,
    compression_maximum = 32767
};

typedef struct htsFormat {
    enum htsFormatCategory category;
    enum htsExactFormat format;
    struct { short major, minor; } version;
    enum htsCompression compression;
    short compression_level;  // currently unused
    void *specific;  // currently unused
} htsFormat;

// Maintainers note htsFile cannot be an opaque structure because some of its
// fields are part of libhts.so's ABI (hence these fields must not be moved):
//  - fp is used in the public sam_itr_next()/etc macros
//  - is_bin is used directly in samtools <= 1.1 and bcftools <= 1.1
//  - is_write and is_cram are used directly in samtools <= 1.1
//  - fp is used directly in samtools (up to and including current develop)
//  - line is used directly in bcftools (up to and including current develop)
typedef struct {
    uint32_t is_bin:1, is_write:1, is_be:1, is_cram:1, dummy:28;
    int64_t lineno;
    kstring_t line;
    char *fn, *fn_aux;
    union {
        BGZF *bgzf;
        struct cram_fd *cram;
        struct hFILE *hfile;
        void *voidp;
    } fp;
    htsFormat format;
} htsFile;

// REQUIRED_FIELDS
enum sam_fields {
    SAM_QNAME = 0x00000001,
    SAM_FLAG  = 0x00000002,
    SAM_RNAME = 0x00000004,
    SAM_POS   = 0x00000008,
    SAM_MAPQ  = 0x00000010,
    SAM_CIGAR = 0x00000020,
    SAM_RNEXT = 0x00000040,
    SAM_PNEXT = 0x00000080,
    SAM_TLEN  = 0x00000100,
    SAM_SEQ   = 0x00000200,
    SAM_QUAL  = 0x00000400,
    SAM_AUX   = 0x00000800,
    SAM_RGAUX = 0x00001000,
};

enum cram_option {
    CRAM_OPT_DECODE_MD,
    CRAM_OPT_PREFIX,
    CRAM_OPT_VERBOSITY,
    CRAM_OPT_SEQS_PER_SLICE,
    CRAM_OPT_SLICES_PER_CONTAINER,
    CRAM_OPT_RANGE,
    CRAM_OPT_VERSION,
    CRAM_OPT_EMBED_REF,
    CRAM_OPT_IGNORE_MD5,
    CRAM_OPT_REFERENCE,
    CRAM_OPT_MULTI_SEQ_PER_SLICE,
    CRAM_OPT_NO_REF,
    CRAM_OPT_USE_BZIP2,
    CRAM_OPT_SHARED_REF,
    CRAM_OPT_NTHREADS,
    CRAM_OPT_THREAD_POOL,
    CRAM_OPT_USE_LZMA,
    CRAM_OPT_USE_RANS,
    CRAM_OPT_REQUIRED_FIELDS,
};

/**********************
 * Exported functions *
 **********************/

extern int hts_verbose;

/*! @abstract Table for converting a nucleotide character to 4-bit encoding.
The input character may be either an IUPAC ambiguity code, '=' for 0, or
'0'/'1'/'2'/'3' for a result of 1/2/4/8.  The result is encoded as 1/2/4/8
for A/C/G/T or combinations of these bits for ambiguous bases.
*/
extern const unsigned char seq_nt16_table[256];

/*! @abstract Table for converting a 4-bit encoded nucleotide to an IUPAC
ambiguity code letter (or '=' when given 0).
*/
extern const char seq_nt16_str[];

/*! @abstract Table for converting a 4-bit encoded nucleotide to about 2 bits.
Returns 0/1/2/3 for 1/2/4/8 (i.e., A/C/G/T), or 4 otherwise (0 or ambiguous).
*/
extern const int seq_nt16_int[];

#ifdef __cplusplus
extern "C" {
#endif

/*!
  @abstract  Get the htslib version number
  @return    For released versions, a string like "N.N[.N]"; or git describe
  output if using a library built within a Git repository.
*/
const char *hts_version(void);

/*!
  @abstract    Determine format by peeking at the start of a file
  @param fp    File opened for reading, positioned at the beginning
  @param fmt   Format structure that will be filled out on return
  @return      0 for success, or negative if an error occurred.
*/
int hts_detect_format(struct hFILE *fp, htsFormat *fmt);

/*!
  @abstract    Get a human-readable description of the file format
  @return      Description string, to be freed by the caller after use.
*/
char *hts_format_description(const htsFormat *format);

/*!
  @abstract       Open a SAM/BAM/CRAM/VCF/BCF/etc file
  @param fn       The file name or "-" for stdin/stdout
  @param mode     Mode matching /[rwa][bcuz0-9]+/
  @discussion
      With 'r' opens for reading; any further format mode letters are ignored
      as the format is detected by checking the first few bytes or BGZF blocks
      of the file.  With 'w' or 'a' opens for writing or appending, with format
      specifier letters:
        b  binary format (BAM, BCF, etc) rather than text (SAM, VCF, etc)
        c  CRAM format
        g  gzip compressed
        u  uncompressed
        z  bgzf compressed
        [0-9]  zlib compression level
      Note that there is a distinction between 'u' and '0': the first yields
      plain uncompressed output whereas the latter outputs uncompressed data
      wrapped in the zlib format.
  @example
      [rw]b .. compressed BCF, BAM, FAI
      [rw]u .. uncompressed BCF
      [rw]z .. compressed VCF
      [rw]  .. uncompressed VCF
*/
htsFile *hts_open(const char *fn, const char *mode);

/*!
  @abstract       Open an existing stream as a SAM/BAM/CRAM/VCF/BCF/etc file
  @param fn       The already-open file handle
  @param mode     Open mode, as per hts_open()
*/
htsFile *hts_hopen(struct hFILE *fp, const char *fn, const char *mode);

/*!
  @abstract  Close a file handle, flushing buffered data for output streams
  @param fp  The file handle to be closed
  @return    0 for success, or negative if an error occurred.
*/
int hts_close(htsFile *fp);

/*!
  @abstract  Returns the file's format information
  @param fp  The file handle
  @return    Read-only pointer to the file's htsFormat.
*/
const htsFormat *hts_get_format(htsFile *fp);

/*!
  @abstract  Sets a specified CRAM option on the open file handle.
  @param fp  The file handle open the open file.
  @param opt The CRAM_OPT_* option.
  @param ... Optional arguments, dependent on the option used.
  @return    0 for success, or negative if an error occurred.
*/
int hts_set_opt(htsFile *fp, enum cram_option opt, ...);

int hts_getline(htsFile *fp, int delimiter, kstring_t *str);
char **hts_readlines(const char *fn, int *_n);
/*!
    @abstract       Parse comma-separated list or read list from a file
    @param list     File name or comma-separated list
    @param is_file
    @param _n       Size of the output array (number of items read)
    @return         NULL on failure or pointer to newly allocated array of
                    strings
*/
char **hts_readlist(const char *fn, int is_file, int *_n);

/*!
  @abstract  Create extra threads to aid compress/decompression for this file
  @param fp  The file handle
  @param n   The number of worker threads to create
  @return    0 for success, or negative if an error occurred.
  @notes     THIS THREADING API IS LIKELY TO CHANGE IN FUTURE.
*/
int hts_set_threads(htsFile *fp, int n);

/*!
  @abstract  Set .fai filename for a file opened for reading
  @return    0 for success, negative on failure
  @discussion
      Called before *_hdr_read(), this provides the name of a .fai file
      used to provide a reference list if the htsFile contains no @SQ headers.
*/
int hts_set_fai_filename(htsFile *fp, const char *fn_aux);

#ifdef __cplusplus
}
#endif

/************
 * Indexing *
 ************/

/*!
These HTS_IDX_* macros are used as special tid values for hts_itr_query()/etc,
producing iterators operating as follows:
 - HTS_IDX_NOCOOR iterates over unmapped reads sorted at the end of the file
 - HTS_IDX_START  iterates over the entire file
 - HTS_IDX_REST   iterates from the current position to the end of the file
 - HTS_IDX_NONE   always returns "no more alignment records"
When one of these special tid values is used, beg and end are ignored.
When REST or NONE is used, idx is also ignored and may be NULL.
*/
#define HTS_IDX_NOCOOR (-2)
#define HTS_IDX_START  (-3)
#define HTS_IDX_REST   (-4)
#define HTS_IDX_NONE   (-5)

#define HTS_FMT_CSI 0
#define HTS_FMT_BAI 1
#define HTS_FMT_TBI 2
#define HTS_FMT_CRAI 3

struct __hts_idx_t;
typedef struct __hts_idx_t hts_idx_t;

typedef struct {
    uint64_t u, v;
} hts_pair64_t;

typedef int hts_readrec_func(BGZF *fp, void *data, void *r, int *tid, int *beg, int *end);

typedef struct {
    uint32_t read_rest:1, finished:1, dummy:29;
    int tid, beg, end, n_off, i;
    int curr_tid, curr_beg, curr_end;
    uint64_t curr_off;
    hts_pair64_t *off;
    hts_readrec_func *readrec;
    struct {
        int n, m;
        int *a;
    } bins;
} hts_itr_t;

#ifdef __cplusplus
extern "C" {
#endif

    #define hts_bin_first(l) (((1<<(((l)<<1) + (l))) - 1) / 7)
    #define hts_bin_parent(l) (((l) - 1) >> 3)

    hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls);
    void hts_idx_destroy(hts_idx_t *idx);
    int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int is_mapped);
    void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset);

    void hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt);
    hts_idx_t *hts_idx_load(const char *fn, int fmt);

    uint8_t *hts_idx_get_meta(hts_idx_t *idx, int *l_meta);
    void hts_idx_set_meta(hts_idx_t *idx, int l_meta, uint8_t *meta, int is_copy);

    int hts_idx_get_stat(const hts_idx_t* idx, int tid, uint64_t* mapped, uint64_t* unmapped);
    uint64_t hts_idx_get_n_no_coor(const hts_idx_t* idx);

    const char *hts_parse_reg(const char *s, int *beg, int *end);
    hts_itr_t *hts_itr_query(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec);
    void hts_itr_destroy(hts_itr_t *iter);

    typedef int (*hts_name2id_f)(void*, const char*);
    typedef const char *(*hts_id2name_f)(void*, int);
    typedef hts_itr_t *hts_itr_query_func(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec);

    hts_itr_t *hts_itr_querys(const hts_idx_t *idx, const char *reg, hts_name2id_f getid, void *hdr, hts_itr_query_func *itr_query, hts_readrec_func *readrec);
    int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data);
    const char **hts_idx_seqnames(const hts_idx_t *idx, int *n, hts_id2name_f getid, void *hdr); // free only the array, not the values

    /**
     * hts_file_type() - Convenience function to determine file type
     * DEPRECATED:  This function has been replaced by hts_detect_format().
     * It and these FT_* macros will be removed in a future HTSlib release.
     */
    #define FT_UNKN   0
    #define FT_GZ     1
    #define FT_VCF    2
    #define FT_VCF_GZ (FT_GZ|FT_VCF)
    #define FT_BCF    (1<<2)
    #define FT_BCF_GZ (FT_GZ|FT_BCF)
    #define FT_STDIN  (1<<3)
    int hts_file_type(const char *fname);


#ifdef __cplusplus
}
#endif

static inline int hts_reg2bin(int64_t beg, int64_t end, int min_shift, int n_lvls)
{
    int l, s = min_shift, t = ((1<<((n_lvls<<1) + n_lvls)) - 1) / 7;
    for (--end, l = n_lvls; l > 0; --l, s += 3, t -= 1<<((l<<1)+l))
        if (beg>>s == end>>s) return t + (beg>>s);
    return 0;
}

static inline int hts_bin_bot(int bin, int n_lvls)
{
    int l, b;
    for (l = 0, b = bin; b; ++l, b = hts_bin_parent(b)); // compute the level of bin
    return (bin - hts_bin_first(l)) << (n_lvls - l) * 3;
}

/**************
 * Endianness *
 **************/

static inline int ed_is_big(void)
{
    long one= 1;
    return !(*((char *)(&one)));
}
static inline uint16_t ed_swap_2(uint16_t v)
{
    return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}
static inline void *ed_swap_2p(void *x)
{
    *(uint16_t*)x = ed_swap_2(*(uint16_t*)x);
    return x;
}
static inline uint32_t ed_swap_4(uint32_t v)
{
    v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
    return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}
static inline void *ed_swap_4p(void *x)
{
    *(uint32_t*)x = ed_swap_4(*(uint32_t*)x);
    return x;
}
static inline uint64_t ed_swap_8(uint64_t v)
{
    v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
    v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
    return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}
static inline void *ed_swap_8p(void *x)
{
    *(uint64_t*)x = ed_swap_8(*(uint64_t*)x);
    return x;
}

#endif
