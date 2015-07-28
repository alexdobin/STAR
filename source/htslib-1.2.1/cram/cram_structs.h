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

#ifndef _CRAM_STRUCTS_H_
#define _CRAM_STRUCTS_H_

/*
 * Defines in-memory structs for the basic file-format objects in the
 * CRAM format.
 *
 * The basic file format is:
 *     File-def SAM-hdr Container Container ...
 *
 * Container:
 *     Service-block data-block data-block ...
 *
 * Multiple blocks in a container are grouped together as slices,
 * also sometimes referred to as landmarks in the spec.
 */


#include <stdint.h>

#include "cram/thread_pool.h"
#include "cram/string_alloc.h"
#include "htslib/khash.h"
#ifdef _WIN32
#include <mutex>
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Generic hash-map integer -> integer
KHASH_MAP_INIT_INT(m_i2i, int)

// Generic hash-set integer -> (existance)
KHASH_SET_INIT_INT(s_i2i)

// For brevity
typedef unsigned char uc;

/*
 * A union for the preservation map. Required for khash.
 */
typedef union {
    int i;
    char *p;
} pmap_t;

// Generates static functions here which isn't ideal, but we have no way
// currently to declare the kh_map_t structure here without also declaring a
// duplicate in the .c files due to the nature of the KHASH macros.
KHASH_MAP_INIT_STR(map, pmap_t)

struct hFILE;

#define SEQS_PER_SLICE 10000
#define SLICE_PER_CNT  1

#define CRAM_SUBST_MATRIX "CGTNAGTNACTNACGNACGT"

#define MAX_STAT_VAL 1024
//#define MAX_STAT_VAL 16
typedef struct {
    int freqs[MAX_STAT_VAL];
    khash_t(m_i2i) *h;
    int nsamp; // total number of values added
    int nvals; // total number of unique values added
} cram_stats;

/* NB: matches java impl, not the spec */
enum cram_encoding {
    E_NULL               = 0,
    E_EXTERNAL           = 1,
    E_GOLOMB             = 2,
    E_HUFFMAN            = 3,
    E_BYTE_ARRAY_LEN     = 4,
    E_BYTE_ARRAY_STOP    = 5,
    E_BETA               = 6,
    E_SUBEXP             = 7,
    E_GOLOMB_RICE        = 8,
    E_GAMMA              = 9
};

enum cram_external_type {
    E_INT                = 1,
    E_LONG               = 2,
    E_BYTE               = 3,
    E_BYTE_ARRAY         = 4,
    E_BYTE_ARRAY_BLOCK   = 5,
};

/* External IDs used by this implementation (only assumed during writing) */
enum cram_DS_ID {
    DS_CORE   = 0,
    DS_aux    = 1, // aux_blk
    DS_aux_OQ = 2,
    DS_aux_BQ = 3,
    DS_aux_BD = 4,
    DS_aux_BI = 5,
    DS_aux_FZ = 6, // also ZM:B
    DS_aux_oq = 7, // other qualities
    DS_aux_os = 8, // other sequences
    DS_aux_oz = 9, // other strings
    DS_ref,
    DS_RN, // name_blk
    DS_QS, // qual_blk
    DS_IN, // base_blk
    DS_SC, // soft_blk

    DS_BF, // start loop
    DS_CF,
    DS_AP,
    DS_RG,
    DS_MQ,
    DS_NS,
    DS_MF,
    DS_TS,
    DS_NP,
    DS_NF,
    DS_RL,
    DS_FN,
    DS_FC,
    DS_FP,
    DS_DL,
    DS_BA,
    DS_BS,
    DS_TL,
    DS_RI,
    DS_RS,
    DS_PD,
    DS_HC,
    DS_BB,
    DS_QQ,

    DS_TN, // end loop

    DS_RN_len,
    DS_SC_len,
    DS_BB_len,
    DS_QQ_len,

    DS_TC, // CRAM v1.0 tags
    DS_TM, // test
    DS_TV, // test
    
    DS_END,
};

/* "File Definition Structure" */
typedef struct {
    char    magic[4];
    uint8_t major_version;
    uint8_t minor_version;
    char    file_id[20];      // Filename or SHA1 checksum
} cram_file_def;

#define CRAM_MAJOR_VERS(v) ((v) >> 8)
#define CRAM_MINOR_VERS(v) ((v) & 0xff)

struct cram_slice;

enum cram_block_method {
    ERROR    = -1,
    RAW      = 0,
    GZIP     = 1,
    BZIP2    = 2,
    LZMA     = 3,
    RANS     = 4,  // Generic; either order
    RANS0    = 4,
    RANS1    = 10, // Not externalised; stored as RANS (generic)
    GZIP_RLE = 11, // NB: not externalised in CRAM
};

enum cram_content_type {
    CT_ERROR           = -1,
    FILE_HEADER        = 0,
    COMPRESSION_HEADER = 1,
    MAPPED_SLICE       = 2,
    UNMAPPED_SLICE     = 3, // CRAM V1.0 only
    EXTERNAL           = 4,
    CORE               = 5,
};

/* Compression metrics */
typedef struct {
    // number of trials and time to next trial
    int trial;
    int next_trial;

    // aggregate sizes during trials
    int sz_gz_rle;
    int sz_gz_def;
    int sz_rans0;
    int sz_rans1;
    int sz_bzip2;
    int sz_lzma;

    // resultant method from trials
    int method;
    int strat;

    // Revisions of method, to allow culling of continually failing ones.
    int gz_rle_cnt;
    int gz_def_cnt;
    int rans0_cnt;
    int rans1_cnt;
    int bzip2_cnt;
    int lzma_cnt;
    int revised_method;

    double gz_rle_extra;
    double gz_def_extra;
    double rans0_extra;
    double rans1_extra;
    double bzip2_extra;
    double lzma_extra;
} cram_metrics;

/* Block */
typedef struct {
    enum cram_block_method  method, orig_method;
    enum cram_content_type  content_type;
    int32_t  content_id;
    int32_t  comp_size;
    int32_t  uncomp_size;
    uint32_t crc32;
    int32_t  idx; /* offset into data */
    unsigned char    *data;

    // For bit I/O
    size_t alloc;
    size_t byte;
    int bit;
} cram_block;

struct cram_codec; /* defined in cram_codecs.h */
struct cram_map;

#define CRAM_MAP_HASH 32
#define CRAM_MAP(a,b) (((a)*3+(b))&(CRAM_MAP_HASH-1))

/* Compression header block */
typedef struct {
    int32_t ref_seq_id;
    int32_t ref_seq_start;
    int32_t ref_seq_span;
    int32_t num_records;
    int32_t num_landmarks;
    int32_t *landmark;

    /* Flags from preservation map */
    int mapped_qs_included;
    int unmapped_qs_included;
    int unmapped_placed;
    int qs_included;
    int read_names_included;
    int AP_delta;
    // indexed by ref-base and subst. code
    char substitution_matrix[5][4];

    // TD Dictionary as a concatenated block
    cram_block *TD_blk;          // Tag Dictionary
    int nTL;		         // number of TL entries in TD
    unsigned char **TL;          // array of size nTL, pointer into TD_blk.
    khash_t(m_s2i) *TD_hash;     // Keyed on TD strings, map to TL[] indices
    string_alloc_t *TD_keys;     // Pooled keys for TD hash.
    
    khash_t(map) *preservation_map;
    struct cram_map *rec_encoding_map[CRAM_MAP_HASH];
    struct cram_map *tag_encoding_map[CRAM_MAP_HASH];

    struct cram_codec *codecs[DS_END];

    char *uncomp; // A single block of uncompressed data
    size_t uncomp_size, uncomp_alloc;

    unsigned int data_series; // See cram_fields enum below
} cram_block_compression_hdr;

typedef struct cram_map {
    int key;    /* 0xe0 + 3 bytes */
    enum cram_encoding encoding;
    int offset; /* Offset into a single block of memory */
    int size;   /* Size */
    struct cram_codec *codec;
    struct cram_map *next; // for noddy internal hash
} cram_map;

/* Mapped or unmapped slice header block */
typedef struct {
    enum cram_content_type content_type;
    int32_t ref_seq_id;     /* if content_type == MAPPED_SLICE */
    int32_t ref_seq_start;  /* if content_type == MAPPED_SLICE */
    int32_t ref_seq_span;   /* if content_type == MAPPED_SLICE */
    int32_t num_records;
    int64_t record_counter;
    int32_t num_blocks;
    int32_t num_content_ids;
    int32_t *block_content_ids;
    int32_t ref_base_id;    /* if content_type == MAPPED_SLICE */
    unsigned char md5[16];
} cram_block_slice_hdr;

struct ref_entry;

/*
 * Container.
 *
 * Conceptually a container is split into slices, and slices into blocks.
 * However on disk it's just a list of blocks and we need to query the
 * block types to identify the start/end points of the slices.
 *
 * OR... are landmarks the start/end points of slices?
 */
typedef struct {
    int32_t  length;
    int32_t  ref_seq_id;
    int32_t  ref_seq_start;
    int32_t  ref_seq_span;
    int64_t  record_counter;
    int64_t  num_bases;
    int32_t  num_records;
    int32_t  num_blocks;
    int32_t  num_landmarks;
    int32_t *landmark;

    /* Size of container header above */
    size_t   offset;
    
    /* Compression header is always the first block? */
    cram_block_compression_hdr *comp_hdr;
    cram_block *comp_hdr_block;

    /* For construction purposes */
    int max_slice, curr_slice;   // maximum number of slices
    int max_rec, curr_rec;       // current and max recs per slice
    int max_c_rec, curr_c_rec;   // current and max recs per container
    int slice_rec;               // rec no. for start of this slice
    int curr_ref;                // current ref ID. -2 for no previous
    int last_pos;                // last record position
    struct cram_slice **slices, *slice;
    int pos_sorted;              // boolean, 1=>position sorted data
    int max_apos;                // maximum position, used if pos_sorted==0
    int last_slice;              // number of reads in last slice (0 for 1st)
    int multi_seq;               // true if packing multi seqs per cont/slice
    int unsorted;		 // true is AP_delta is 0.

    /* Copied from fd before encoding, to allow multi-threading */
    int ref_start, first_base, last_base, ref_id, ref_end;
    char *ref;
    //struct ref_entry *ref;

    /* For multi-threading */
    bam_seq_t **bams;

    /* Statistics for encoding */
    cram_stats *stats[DS_END];

    khash_t(s_i2i) *tags_used; // set of tag types in use, for tag encoding map
    int *refs_used;       // array of frequency of ref seq IDs

    uint32_t crc32;       // CRC32
} cram_container;

/*
 * A single cram record
 */
typedef struct {
    struct cram_slice *s; // Filled out by cram_decode only

    int32_t ref_id;       // fixed for all recs in slice?
    int32_t flags;        // BF
    int32_t cram_flags;   // CF
    int32_t len;          // RL
    int32_t apos;         // AP
    int32_t rg;           // RG
    int32_t name;         // RN; idx to s->names_blk
    int32_t name_len;
    int32_t mate_line;    // index to another cram_record
    int32_t mate_ref_id;
    int32_t mate_pos;     // NP
    int32_t tlen;         // TS

    // Auxiliary data
    int32_t ntags;        // TC
    int32_t aux;          // idx to s->aux_blk
    int32_t aux_size;     // total size of packed ntags in aux_blk
#ifndef TN_external
    int32_t TN_idx;       // TN; idx to s->TN;
#else
    int32_t tn;           // idx to s->tn_blk
#endif
    int     TL;

    int32_t seq;          // idx to s->seqs_blk
    int32_t qual;         // idx to s->qual_blk
    int32_t cigar;        // idx to s->cigar
    int32_t ncigar;
    int32_t aend;         // alignment end
    int32_t mqual;        // MQ

    int32_t feature;      // idx to s->feature
    int32_t nfeature;     // number of features
    int32_t mate_flags;   // MF
} cram_record;

// Accessor macros as an analogue of the bam ones
#define cram_qname(c)    (&(c)->s->name_blk->data[(c)->name])
#define cram_seq(c)      (&(c)->s->seqs_blk->data[(c)->seq])
#define cram_qual(c)     (&(c)->s->qual_blk->data[(c)->qual])
#define cram_aux(c)      (&(c)->s->aux_blk->data[(c)->aux])
#define cram_seqi(c,i)   (cram_seq((c))[(i)])
#define cram_name_len(c) ((c)->name_len)
#define cram_strand(c)   (((c)->flags & BAM_FREVERSE) != 0)
#define cram_mstrand(c)  (((c)->flags & BAM_FMREVERSE) != 0)
#define cram_cigar(c)    (&((cr)->s->cigar)[(c)->cigar])

/*
 * A feature is a base difference, used for the sequence reference encoding.
 * (We generate these internally when writing CRAM.)
 */
typedef struct {
    union {
	struct {
	    int pos;
	    int code;
	    int base;    // substitution code
	} X;
	struct {
	    int pos;
	    int code;
	    int base;    // actual base & qual
	    int qual;
	} B;
	struct {
	    int pos;
	    int code;
	    int seq_idx; // index to s->seqs_blk
	    int len;
	} b;
	struct {
	    int pos;
	    int code;
	    int qual;
	} Q;
	struct {
	    int pos;
	    int code;
	    int len;
	    int seq_idx; // soft-clip multiple bases
	} S;
	struct {
	    int pos;
	    int code;
	    int len;
	    int seq_idx; // insertion multiple bases
	} I;
	struct {
	    int pos;
	    int code;
	    int base; // insertion single base
	} i;
	struct {
	    int pos;
	    int code;
	    int len;
	} D;
	struct {
	    int pos;
	    int code;
	    int len;
	} N;
	struct {
	    int pos;
	    int code;
	    int len;
	} P;
	struct {
	    int pos;
	    int code;
	    int len;
	} H;
    };
} cram_feature;

/*
 * A slice is really just a set of blocks, but it
 * is the logical unit for decoding a number of
 * sequences.
 */
typedef struct cram_slice {
    cram_block_slice_hdr *hdr;
    cram_block *hdr_block;
    cram_block **block;
    cram_block **block_by_id;

    /* State used during encoding/decoding */
    int last_apos, max_apos;

    /* Array of decoded cram records */
    cram_record *crecs;

    /* An dynamically growing buffers for data pointed
     * to by crecs[] array.
     */
    uint32_t  *cigar;
    uint32_t   cigar_alloc;
    uint32_t   ncigar;

    cram_feature *features;
    int           nfeatures;
    int           afeatures; // allocated size of features

#ifndef TN_external
    // TN field (Tag Name)
    uint32_t      *TN;
    int           nTN, aTN;  // used and allocated size for TN[]
#else
    cram_block *tn_blk;
    int tn_id;
#endif

    // For variable sized elements which are always external blocks.
    cram_block *name_blk;
    cram_block *seqs_blk;
    cram_block *qual_blk;
    cram_block *base_blk;
    cram_block *soft_blk;
    cram_block *aux_blk;
    cram_block *aux_OQ_blk;
    cram_block *aux_BQ_blk;
    cram_block *aux_BD_blk;
    cram_block *aux_BI_blk;
    cram_block *aux_FZ_blk;
    cram_block *aux_oq_blk;
    cram_block *aux_os_blk;
    cram_block *aux_oz_blk;

    string_alloc_t *pair_keys; // Pooled keys for pair hash.
    khash_t(m_s2i) *pair[2];   // for identifying read-pairs in this slice.

    char *ref;               // slice of current reference
    int ref_start;           // start position of current reference;
    int ref_end;             // end position of current reference;
    int ref_id;
} cram_slice;

/*-----------------------------------------------------------------------------
 * Consider moving reference handling to cram_refs.[ch]
 */
// from fa.fai / samtools faidx files
typedef struct ref_entry {
    char *name;
    char *fn;
    int64_t length;
    int64_t offset;
    int bases_per_line;
    int line_length;
    int64_t count;	   // for shared references so we know to dealloc seq
    char *seq;
} ref_entry;

KHASH_MAP_INIT_STR(refs, ref_entry*)

// References structure.
typedef struct {
    string_alloc_t *pool;  // String pool for holding filenames and SN vals

    khash_t(refs) *h_meta; // ref_entry*, index by name
    ref_entry **ref_id;    // ref_entry*, index by ID
    int nref;              // number of ref_entry

    char *fn;              // current file opened
    BGZF *fp;              // and the hFILE* to go with it.

    int count;             // how many cram_fd sharing this refs struct

#ifndef _WIN32
    pthread_mutex_t lock;  // Mutex for multi-threaded updating
#else
	std::mutex lock;
#endif
	
    ref_entry *last;       // Last queried sequence
    int last_id;           // Used in cram_ref_decr_locked to delay free
} refs_t;

/*-----------------------------------------------------------------------------
 * CRAM index
 *
 * Detect format by number of entries per line.
 * 5 => 1.0 (refid, start, nseq, C offset, slice)
 * 6 => 1.1 (refid, start, span, C offset, S offset, S size)
 *
 * Indices are stored in a nested containment list, which is trivial to set
 * up as the indices are on sorted data so we're appending to the nclist
 * in sorted order. Basically if a slice entirely fits within a previous
 * slice then we append to that slices list. This is done recursively.
 *
 * Lists are sorted on two dimensions: ref id + slice coords.
 */
typedef struct cram_index {
    int nslice, nalloc;   // total number of slices
    struct cram_index *e; // array of size nslice

    int     refid;  // 1.0                 1.1
    int     start;  // 1.0                 1.1
    int     end;    //                     1.1
    int     nseq;   // 1.0 - undocumented
    int     slice;  // 1.0 landmark index, 1.1 landmark value
    int     len;    //                     1.1 - size of slice in bytes
    int64_t offset; // 1.0                 1.1
} cram_index;

typedef struct {
    int refid;
    int start;
    int end;
} cram_range;

/*-----------------------------------------------------------------------------
 */
/* CRAM File handle */

typedef struct spare_bams {
    bam_seq_t **bams;
    struct spare_bams *next;
} spare_bams;

typedef struct cram_fd {
    struct hFILE  *fp;
    int            mode;     // 'r' or 'w'
    int            version;
    cram_file_def *file_def;
    SAM_hdr       *header;

    char          *prefix;
    int64_t        record_counter;
    int            err;

    // Most recent compression header decoded
    //cram_block_compression_hdr *comp_hdr;
    //cram_block_slice_hdr       *slice_hdr;

    // Current container being processed.
    cram_container *ctr;

    // positions for encoding or decoding
    int first_base, last_base;

    // cached reference portion
    refs_t *refs;              // ref meta-data structure
    char *ref, *ref_free;      // current portion held in memory
    int   ref_id;
    int   ref_start;
    int   ref_end;
    char *ref_fn;   // reference fasta filename

    // compression level and metrics
    int level;
    cram_metrics *m[DS_END];

    // options
    int decode_md; // Whether to export MD and NM tags
    int verbose;
    int seqs_per_slice;
    int slices_per_container;
    int embed_ref;
    int no_ref;
    int ignore_md5;
    int use_bz2;
    int use_rans;
    int use_lzma;
    int shared_ref;
    unsigned int required_fields;
    cram_range range;

    // lookup tables, stored here so we can be trivially multi-threaded
    unsigned int bam_flag_swap[0x1000]; // cram -> bam flags
    unsigned int cram_flag_swap[0x1000];// bam -> cram flags
    unsigned char L1[256];              // ACGT{*} ->0123{4}
    unsigned char L2[256];              // ACGTN{*}->01234{5}
    char cram_sub_matrix[32][32];	// base substituion codes

    int         index_sz;
    cram_index *index;                  // array, sizeof index_sz
    off_t first_container;
    int eof;
    int last_slice;                     // number of recs encoded in last slice
    int multi_seq;
    int unsorted;
    int empty_container; 		// Marker for EOF block
    
    // thread pool
    int own_pool;
    t_pool *pool;
    t_results_queue *rqueue;

#if !defined(_WIN32) && defined(USE_PTHREAD)
    pthread_mutex_t metrics_lock;
    pthread_mutex_t ref_lock;
    pthread_mutex_t bam_list_lock;
#else // ~ #if !defined(_WIN32) && defined(USE_PTHREAD)
	std::mutex metrics_lock;
	std::mutex ref_lock;
	std::mutex bam_list_lock;
#endif // ~ #if !defined(_WIN32) && defined(USE_PTHREAD)

	spare_bams *bl;
    void *job_pending;
    int ooc;                            // out of containers.
} cram_fd;

// Translation of required fields to cram data series
enum cram_fields {
    CRAM_BF = 0x00000001,
    CRAM_AP = 0x00000002,
    CRAM_FP = 0x00000004,
    CRAM_RL = 0x00000008,
    CRAM_DL = 0x00000010,
    CRAM_NF = 0x00000020,
    CRAM_BA = 0x00000040,
    CRAM_QS = 0x00000080,
    CRAM_FC = 0x00000100,
    CRAM_FN = 0x00000200,
    CRAM_BS = 0x00000400,
    CRAM_IN = 0x00000800,
    CRAM_RG = 0x00001000,
    CRAM_MQ = 0x00002000,
    CRAM_TL = 0x00004000,
    CRAM_RN = 0x00008000,
    CRAM_NS = 0x00010000,
    CRAM_NP = 0x00020000,
    CRAM_TS = 0x00040000,
    CRAM_MF = 0x00080000,
    CRAM_CF = 0x00100000,
    CRAM_RI = 0x00200000,
    CRAM_RS = 0x00400000,
    CRAM_PD = 0x00800000,
    CRAM_HC = 0x01000000,
    CRAM_SC = 0x02000000,
    CRAM_BB = 0x04000000,
    CRAM_BB_len = 0x08000000,
    CRAM_QQ = 0x10000000,
    CRAM_QQ_len = 0x20000000,
    CRAM_aux= 0x40000000,
    CRAM_ALL= 0x7fffffff,
};

// A CIGAR opcode, but not necessarily the implications of it. Eg FC/FP may
// encode a base difference, but we don't need to know what it is for CIGAR.
// If we have a soft-clip or insertion, we do need SC/IN though to know how
// long that array is.
#define CRAM_CIGAR (CRAM_FN | CRAM_FP | CRAM_FC | CRAM_DL | CRAM_IN | \
		    CRAM_SC | CRAM_HC | CRAM_PD | CRAM_RS | CRAM_RL | CRAM_BF)

#define CRAM_SEQ (CRAM_CIGAR | CRAM_BA | CRAM_QS | CRAM_BS | \
		  CRAM_RL    | CRAM_AP | CRAM_BB | CRAM_QQ)

/* BF bitfields */
/* Corrected in 1.1. Use bam_flag_swap[bf] and BAM_* macros for 1.0 & 1.1 */
#define CRAM_FPAIRED      256
#define CRAM_FPROPER_PAIR 128
#define CRAM_FUNMAP        64
#define CRAM_FREVERSE      32
#define CRAM_FREAD1        16
#define CRAM_FREAD2         8
#define CRAM_FSECONDARY     4
#define CRAM_FQCFAIL        2
#define CRAM_FDUP           1

#define DS_aux_S "\001"
#define DS_aux_OQ_S "\002"
#define DS_aux_BQ_S "\003"
#define DS_aux_BD_S "\004"
#define DS_aux_BI_S "\005"
#define DS_aux_FZ_S "\006"
#define DS_aux_oq_S "\007"
#define DS_aux_os_S "\010"
#define DS_aux_oz_S "\011"

#define CRAM_M_REVERSE  1
#define CRAM_M_UNMAP    2


/* CF bitfields */
#define CRAM_FLAG_PRESERVE_QUAL_SCORES (1<<0)
#define CRAM_FLAG_DETACHED             (1<<1)
#define CRAM_FLAG_MATE_DOWNSTREAM      (1<<2)

#ifdef __cplusplus
}
#endif

#endif /* _CRAM_STRUCTS_H_ */
