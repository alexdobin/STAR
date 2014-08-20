#ifndef HTS_H
#define HTS_H

#include <stddef.h>
#include <stdint.h>

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

typedef struct {
	uint32_t is_bin:1, is_write:1, is_be:1, is_cram:1, is_compressed:2, is_kstream:1, dummy:25;
	int64_t lineno;
	kstring_t line;
	char *fn, *fn_aux;
	union {
		BGZF *bgzf;
		struct cram_fd *cram;
		struct hFILE *hfile;
		void *voidp;
	} fp;
} htsFile;

/**********************
 * Exported functions *
 **********************/

extern int hts_verbose;

/*! @abstract Table for converting a nucleotide character to the 4-bit encoding. */
extern const unsigned char seq_nt16_table[256];

/*! @abstract Table for converting a 4-bit encoded nucleotide to a letter. */
extern const char seq_nt16_str[];

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
        u  uncompressed
        z  compressed
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
  @abstract  Close a file handle, flushing buffered data for output streams
  @param fp  The file handle to be closed
  @return    0 for success, or negative if an error occurred.
*/
int hts_close(htsFile *fp);

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
     * @fname: the file name
     *
     * Returns one of the FT_* defines.
     *
     * This function was added in order to avoid the need for excessive command
     * line switches.
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
