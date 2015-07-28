/*  synced_bcf_reader.h -- stream through multiple VCF files.

    Copyright (C) 2012-2014 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

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

/*
    The synced_bcf_reader allows to keep multiple VCFs open and stream them
    using the next_line iterator in a seamless matter without worrying about
    chromosomes and synchronizing the sites. This is used by vcfcheck to
    compare multiple VCFs simultaneously and is used also for merging,
    creating intersections, etc.

    The synced_bcf_reader also provides API for reading indexed BCF/VCF,
    hiding differences in BCF/VCF opening, indexing and reading.


    Example of usage:

        bcf_srs_t *sr = bcf_sr_init();
        for (i=0; i<nfiles; i++)
            bcf_sr_add_reader(sr,files[i]);
        while ( bcf_sr_next_line(sr) )
        {
            for (i=0; i<nfiles; i++)
            {
                bcf1_t *line = bcf_sr_get_line(sr,i);
                ...
            }
        }
        bcf_sr_destroy(sr);
*/

#ifndef HTSLIB_SYNCED_BCF_READER_H
#define HTSLIB_SYNCED_BCF_READER_H

#include "hts.h"
#include "vcf.h"
#include "tbx.h"

// How should be treated sites with the same position but different alleles
#define COLLAPSE_NONE   0   // require the exact same set of alleles in all files
#define COLLAPSE_SNPS   1   // allow different alleles, as long as they all are SNPs
#define COLLAPSE_INDELS 2   // the same as above, but with indels
#define COLLAPSE_ANY    4   // any combination of alleles can be returned by bcf_sr_next_line()
#define COLLAPSE_SOME   8   // at least some of the ALTs must match
#define COLLAPSE_BOTH  (COLLAPSE_SNPS|COLLAPSE_INDELS)

typedef struct _bcf_sr_regions_t
{
    // for reading from tabix-indexed file (big data)
    tbx_t *tbx;             // tabix index
    hts_itr_t *itr;         // tabix iterator
    kstring_t line;         // holder of the current line, set only when reading from tabix-indexed files
    htsFile *file;
    char *fname;
    int is_bin;             // is open in binary mode (tabix access)
    char **als;             // parsed alleles if targets_als set and _regions_match_alleles called
    kstring_t als_str;      // block of parsed alleles
    int nals, mals;         // number of set alleles and the size of allocated array
    int als_type;           // alleles type, currently VCF_SNP or VCF_INDEL

    // user handler to deal with skipped regions without a counterpart in VCFs
    void (*missed_reg_handler)(struct _bcf_sr_regions_t *, void *);
    void *missed_reg_data;

    // for in-memory regions (small data)
    struct _region_t *regs; // the regions

    // shared by both tabix-index and in-memory regions
    void *seq_hash;         // keys: sequence names, values: index to seqs
    char **seq_names;       // sequence names
    int nseqs;              // number of sequences (chromosomes) in the file
    int iseq;               // current position: chr name, index to snames
    int start, end;         // current position: start, end of the region (0-based)
    int prev_seq, prev_start;
}
bcf_sr_regions_t;

typedef struct
{
    htsFile *file;
    tbx_t *tbx_idx;
    hts_idx_t *bcf_idx;
    bcf_hdr_t *header;
    hts_itr_t *itr;
    const char *fname;
    bcf1_t **buffer;                // cached VCF records. First is the current record synced across the reader
    int nbuffer, mbuffer;           // number of cached records (including the current record); number of allocated records
    int nfilter_ids, *filter_ids;   // -1 for ".", otherwise filter id as returned by bcf_id2int
    int *samples, n_smpl;   // list of columns in the order consistent with bcf_srs_t.samples
}
bcf_sr_t;

typedef enum
{
    open_failed, not_bgzf, idx_load_failed, file_type_error, api_usage_error,
    header_error
}
bcf_sr_error;

typedef struct
{
    // Parameters controlling the logic
    int collapse;       // How should the duplicate sites be treated. One of the COLLAPSE_* types above.
    char *apply_filters;    // If set, sites where none of the FILTER strings is listed
                            // will be skipped. Active only at the time of
                            // initialization, that is during the add_reader()
                            // calls. Therefore, each reader can be initialized with different
                            // filters.
    int require_index;  // Some tools do not need random access
    int max_unpack;     // When reading VCFs and knowing some fields will not be needed, boost performance of vcf_parse1
    int *has_line;      // Corresponds to return value of bcf_sr_next_line but is not limited by sizeof(int). Use bcf_sr_has_line macro to query.
    bcf_sr_error errnum;

    // Auxiliary data
    bcf_sr_t *readers;
    int nreaders;
    int streaming;      // reading mode: index-jumping or streaming
    int explicit_regs;  // was the list of regions se by bcf_sr_set_regions or guessed from tabix index?
    char **samples; // List of samples
    bcf_sr_regions_t *regions, *targets;    // see bcf_sr_set_[targets|regions] for description
    int targets_als;    // subset to targets not only by position but also by alleles?
    int targets_exclude;
    kstring_t tmps;
    int n_smpl;
}
bcf_srs_t;

#ifdef __cplusplus
extern "C" {
#endif

/** Init bcf_srs_t struct */
bcf_srs_t *bcf_sr_init(void);

/** Destroy  bcf_srs_t struct */
void bcf_sr_destroy(bcf_srs_t *readers);

char *bcf_sr_strerror(int errnum);


/**
 *  bcf_sr_add_reader() - open new reader
 *  @readers: holder of the open readers
 *  @fname:   the VCF file
 *
 *  Returns 1 if the call succeeded, or 0 on error.
 *
 *  See also the bcf_srs_t data structure for parameters controlling
 *  the reader's logic.
 */
int bcf_sr_add_reader(bcf_srs_t *readers, const char *fname);
void bcf_sr_remove_reader(bcf_srs_t *files, int i);

/**
 * bcf_sr_next_line() - the iterator
 * @readers:    holder of the open readers
 *
 * Returns the number of readers which have the current line
 * (bcf_sr_t.buffer[0]) set at this position. Use the bcf_sr_has_line macro to
 * determine which of the readers are set.
 */
int bcf_sr_next_line(bcf_srs_t *readers);
#define bcf_sr_has_line(readers, i) (readers)->has_line[i]
#define bcf_sr_get_line(_readers, i) ((_readers)->has_line[i] ? ((_readers)->readers[i].buffer[0]) : NULL)
#define bcf_sr_region_done(_readers,i) (!(_readers)->has_line[i] && !(_readers)->readers[i].nbuffer ? 1 : 0)
#define bcf_sr_get_header(_readers, i) (_readers)->readers[i].header
#define bcf_sr_get_reader(_readers, i) &((_readers)->readers[i])

/**
 *  bcf_sr_seek() - set all readers to selected position
 *  @seq:  sequence name; NULL to seek to start
 *  @pos:  0-based coordinate
 */
int bcf_sr_seek(bcf_srs_t *readers, const char *seq, int pos);

/**
 * bcf_sr_set_samples() - sets active samples
 * @readers: holder of the open readers
 * @samples: this can be one of: file name with one sample per line;
 *           or column-separated list of samples; or '-' for a list of
 *           samples shared by all files. If first character is the
 *           exclamation mark, all but the listed samples are included.
 * @is_file: 0: list of samples; 1: file with sample names
 *
 * Returns 1 if the call succeeded, or 0 on error.
 */
int bcf_sr_set_samples(bcf_srs_t *readers, const char *samples, int is_file);

/**
 *  bcf_sr_set_targets(), bcf_sr_set_regions() - init targets/regions
 *  @readers:   holder of the open readers
 *  @targets:   list of regions, one-based and inclusive.
 *  @is_fname:  0: targets is a comma-separated list of regions (chr,chr:from-to)
 *              1: targets is a tabix indexed file with a list of regions
 *              (<chr,pos> or <chr,from,to>)
 *
 *  Returns 0 if the call succeeded, or -1 on error.
 *
 *  Both functions behave the same way, unlisted positions will be skipped by
 *  bcf_sr_next_line(). However, there is an important difference: regions use
 *  index to jump to desired positions while targets streams the whole files
 *  and merely skip unlisted positions.
 *
 *  Moreover, bcf_sr_set_targets() accepts an optional parameter $alleles which
 *  is intepreted as a 1-based column index in the tab-delimited file where
 *  alleles are listed. This in principle enables to perform the COLLAPSE_*
 *  logic also with tab-delimited files. However, the current implementation
 *  considers the alleles merely as a suggestion for prioritizing one of possibly
 *  duplicate VCF lines. It is up to the caller to examine targets->als if
 *  perfect match is sought after. Note that the duplicate positions in targets
 *  file are currently not supported.
 *  Targets (but not regions) can be prefixed with "^" to request logical complement,
 *  for example "^X,Y,MT" indicates that sequences X, Y and MT should be skipped.
 */
int bcf_sr_set_targets(bcf_srs_t *readers, const char *targets, int is_file, int alleles);
int bcf_sr_set_regions(bcf_srs_t *readers, const char *regions, int is_file);



/*
 *  bcf_sr_regions_init()
 *  @regions:   regions can be either a comma-separated list of regions
 *              (chr|chr:pos|chr:from-to|chr:from-) or VCF, BED, or
 *              tab-delimited file (the default). Uncompressed files
 *              are stored in memory while bgzip-compressed and tabix-indexed
 *              region files are streamed.
 *  @is_file:   0: regions is a comma-separated list of regions
 *                  (chr|chr:pos|chr:from-to|chr:from-)
 *              1: VCF, BED or tab-delimited file
 *  @chr, from, to:
 *              Column indexes of chromosome, start position and end position
 *              in the tab-delimited file. The positions are 1-based and
 *              inclusive.
 *              These parameters are ignored when reading from VCF, BED or
 *              tabix-indexed files. When end position column is not present,
 *              supply 'from' in place of 'to'. When 'to' is negative, first
 *              abs(to) will be attempted and if that fails, 'from' will be used
 *              instead.
 */
bcf_sr_regions_t *bcf_sr_regions_init(const char *regions, int is_file, int chr, int from, int to);
void bcf_sr_regions_destroy(bcf_sr_regions_t *regions);

/*
 *  bcf_sr_regions_seek() - seek to the chromosome block
 *
 *  Returns 0 on success or -1 on failure. Sets reg->seq appropriately and
 *  reg->start,reg->end to -1.
 */
int bcf_sr_regions_seek(bcf_sr_regions_t *regions, const char *chr);

/*
 *  bcf_sr_regions_next() - retrieves next region. Returns 0 on success and -1
 *  when all regions have been read. The fields reg->seq, reg->start and
 *  reg->end are filled with the genomic coordinates on succes or with
 *  NULL,-1,-1 when no region is available. The coordinates are 0-based,
 *  inclusive.
 */
int bcf_sr_regions_next(bcf_sr_regions_t *reg);

/*
 *  bcf_sr_regions_overlap() - checks if the interval <start,end> overlaps any of
 *  the regions, the coordinates are 0-based, inclusive. The coordinate queries
 *  must come in ascending order.
 *
 *  Returns 0 if the position is in regions; -1 if the position is not in the
 *  regions and more regions exist; -2 if not in the regions and there are no more
 *  regions left.
 */
int bcf_sr_regions_overlap(bcf_sr_regions_t *reg, const char *seq, int start, int end);

/*
 *  bcf_sr_regions_flush() - calls repeatedly regs->missed_reg_handler() until
 *  all remaining records are processed.
 */
void bcf_sr_regions_flush(bcf_sr_regions_t *regs);

#ifdef __cplusplus
}
#endif

#endif
