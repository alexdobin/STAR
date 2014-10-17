/*
Copyright (c) 2013-2014 Genome Research Ltd.
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
 * SAM header parsing.
 *
 * These functions can be shared between SAM, BAM and CRAM file
 * formats as all three internally use the same string encoding for
 * header fields.
 *
 * Consider using the scram() generic API and calling
 * scram_get_header() to obtain the format-specific pointer to the
 * SAM_hdr struct.
 */ 

/*
 * TODO.
 *
 * - Sort order (parse to struct, enum type, updating funcs)
 * - Removal of lines.
 * - Updating of lines
 */

#ifndef _SAM_HDR_H_
#define _SAM_HDR_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdarg.h>

#include "cram/string_alloc.h"
#include "cram/pooled_alloc.h"

#include "htslib/khash.h"
#include "htslib/kstring.h"

// For structure assignment. Eg kstring_t s = KS_INITIALIZER;
#define KS_INITIALIZER {0,0,0}

// For initialisation elsewhere. Eg KS_INIT(x->str);
#define KS_INIT(ks) ((ks)->l = 0, (ks)->m = 0, (ks)->s = NULL)

// Frees the string subfield only. Assumes 's' itself is static.
#define KS_FREE(ks) do { if ((ks)->s) free((ks)->s); } while(0)

/*
 * Proposed new SAM header parsing

1 @SQ ID:foo LN:100
2 @SQ ID:bar LN:200
3 @SQ ID:ram LN:300 UR:xyz
4 @RG ID:r ...
5 @RG ID:s ...

Hash table for 2-char @keys without dup entries.
If dup lines, we form a circular linked list. Ie hash keys = {RG, SQ}.

HASH("SQ")--\
            |
    (3) <-> 1 <-> 2 <-> 3 <-> (1)

HASH("RG")--\
            |
    (5) <-> 4 <-> 5 <-> (4)

Items stored in the hash values also form their own linked lists:
Ie SQ->ID(foo)->LN(100)
   SQ->ID(bar)->LN(200)
   SQ->ID(ram)->LN(300)->UR(xyz)
   RG->ID(r)
 */

/*! A single key:value pair on a header line
 *
 * These form a linked list and hold strings. The strings are
 * allocated from a string_alloc_t pool referenced in the master
 * SAM_hdr structure. Do not attempt to free, malloc or manipulate
 * these strings directly.
 */
typedef struct SAM_hdr_tag_s {
    struct SAM_hdr_tag_s *next;
    char *str;
    int   len;
} SAM_hdr_tag;

/*! The parsed version of the SAM header string.
 * 
 * Each header type (SQ, RG, HD, etc) points to its own SAM_hdr_type
 * struct via the main hash table h in the SAM_hdr struct.
 *
 * These in turn consist of circular bi-directional linked lists (ie
 * rings) to hold the multiple instances of the same header type
 * code. For example if we have 5 \@SQ lines the primary hash table
 * will key on \@SQ pointing to the first SAM_hdr_type and that in turn
 * will be part of a ring of 5 elements.
 *
 * For each SAM_hdr_type structure we also point to a SAM_hdr_tag
 * structure which holds the tokenised attributes; the tab separated
 * key:value pairs per line.
 */
typedef struct SAM_hdr_item_s {
    struct SAM_hdr_item_s *next; // cirular
    struct SAM_hdr_item_s *prev;
    SAM_hdr_tag *tag;            // first tag
    int order;                   // 0 upwards
} SAM_hdr_type;

/*! Parsed \@SQ lines */
typedef struct {
    char *name;
    uint32_t len;
    SAM_hdr_type *ty;
    SAM_hdr_tag  *tag;
} SAM_SQ;

/*! Parsed \@RG lines */
typedef struct {
    char *name;
    SAM_hdr_type *ty;
    SAM_hdr_tag  *tag;
    int name_len;
    int id;           // numerical ID
} SAM_RG;

/*! Parsed \@PG lines */
typedef struct {
    char *name;
    SAM_hdr_type *ty;
    SAM_hdr_tag  *tag;
    int name_len;
    int id;           // numerical ID
    int prev_id;      // -1 if none
} SAM_PG;

KHASH_MAP_INIT_INT(sam_hdr, SAM_hdr_type*)
KHASH_MAP_INIT_STR(m_s2i, int)

/*! Primary structure for header manipulation
 *
 * The initial header text is held in the text kstring_t, but is also
 * parsed out into SQ, RG and PG arrays. These have a hash table
 * associated with each to allow lookup by ID or SN fields instead of
 * their numeric array indices. Additionally PG has an array to hold
 * the linked list start points (the last in a PP chain).
 *
 * Use the appropriate sam_hdr_* functions to edit the header, and 
 * call sam_hdr_rebuild() any time the textual form needs to be
 * updated again.
 */
typedef struct {
    kstring_t text;           //!< concatenated text, indexed by SAM_hdr_tag
    khash_t(sam_hdr) *h;
    string_alloc_t *str_pool; //!< Pool of SAM_hdr_tag->str strings
    pool_alloc_t   *type_pool;//!< Pool of SAM_hdr_type structs
    pool_alloc_t   *tag_pool; //!< Pool of SAM_hdr_tag structs

    // @SQ lines / references
    int nref;                 //!< Number of \@SQ lines
    SAM_SQ *ref;              //!< Array of parsed \@SQ lines
    khash_t(m_s2i) *ref_hash; //!< Maps SQ SN field to sq[] index

    // @RG lines / read-groups
    int nrg;                  //!< Number of \@RG lines
    SAM_RG *rg;               //!< Array of parsed \@RG lines
    khash_t(m_s2i) *rg_hash;  //!< Maps RG ID field to rg[] index

    // @PG lines / programs
    int npg;                  //!< Number of \@PG lines
    int npg_end;              //!< Number of terminating \@PG lines
    int npg_end_alloc;        //!< Size of pg_end field
    SAM_PG *pg;		      //!< Array of parsed \@PG lines
    khash_t(m_s2i) *pg_hash;  //!< Maps PG ID field to pg[] index
    int *pg_end;              //!< \@PG chain termination IDs

    // @cond internal
    char ID_buf[1024];  // temporary buffer
    int ID_cnt;
    int ref_count;      // number of uses of this SAM_hdr
    // @endcond
} SAM_hdr;

/*! Creates an empty SAM header, ready to be populated.
 * 
 * @return
 * Returns a SAM_hdr struct on success (free with sam_hdr_free())
 *         NULL on failure
 */
SAM_hdr *sam_hdr_new(void);

/*! Tokenises a SAM header into a hash table.
 *
 * Also extracts a few bits on specific data types, such as @RG lines.
 *
 * @return
 * Returns a SAM_hdr struct on success (free with sam_hdr_free());
 *         NULL on failure
 */
#ifdef SAMTOOLS
SAM_hdr *sam_hdr_parse_(const char *hdr, int len);
#else
SAM_hdr *sam_hdr_parse(const char *hdr, int len);
#endif


/*! Produces a duplicate copy of hdr and returns it.
 * @return
 * Returns NULL on failure
 */
SAM_hdr *sam_hdr_dup(SAM_hdr *hdr);


/*! Increments a reference count on hdr.
 *
 * This permits multiple files to share the same header, all calling
 * sam_hdr_free when done, without causing errors for other open  files.
 */
void sam_hdr_incr_ref(SAM_hdr *hdr);


/*! Increments a reference count on hdr.
 *
 * This permits multiple files to share the same header, all calling
 * sam_hdr_free when done, without causing errors for other open  files.
 *
 * If the reference count hits zero then the header is automatically
 * freed. This makes it a synonym for sam_hdr_free().
 */
void sam_hdr_decr_ref(SAM_hdr *hdr);


/*! Deallocates all storage used by a SAM_hdr struct.
 *
 * This also decrements the header reference count. If after decrementing 
 * it is still non-zero then the header is assumed to be in use by another
 * caller and the free is not done.
 *
 * This is a synonym for sam_hdr_dec_ref().
 */
void sam_hdr_free(SAM_hdr *hdr);

/*! Returns the current length of the SAM_hdr in text form.
 *
 * Call sam_hdr_rebuild() first if editing has taken place.
 */
int sam_hdr_length(SAM_hdr *hdr);

/*! Returns the string form of the SAM_hdr.
 *
 * Call sam_hdr_rebuild() first if editing has taken place.
 */
char *sam_hdr_str(SAM_hdr *hdr);

/*! Appends a formatted line to an existing SAM header.
 *
 * Line is a full SAM header record, eg "@SQ\tSN:foo\tLN:100", with
 * optional new-line. If it contains more than 1 line then multiple lines
 * will be added in order.
 *
 * Len is the length of the text data, or 0 if unknown (in which case
 * it should be null terminated).
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int sam_hdr_add_lines(SAM_hdr *sh, const char *lines, int len);

/*! Adds a single line to a SAM header.
 *
 * Specify type and one or more key,value pairs, ending with the NULL key.
 * Eg. sam_hdr_add(h, "SQ", "ID", "foo", "LN", "100", NULL).
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int sam_hdr_add(SAM_hdr *sh, const char *type, ...);

/*! Adds a single line to a SAM header.
 *
 * This is much like sam_hdr_add() but with the additional va_list
 * argument. This is followed by specifying type and one or more
 * key,value pairs, ending with the NULL key.
 *
 * Eg. sam_hdr_vadd(h, "SQ", args, "ID", "foo", "LN", "100", NULL).
 *
 * The purpose of the additional va_list parameter is to permit other
 * varargs functions to call this while including their own additional
 * parameters; an example is in sam_hdr_add_PG().
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int sam_hdr_vadd(SAM_hdr *sh, const char *type, va_list ap, ...);

/*!
 * @return
 * Returns the first header item matching 'type'. If ID is non-NULL it checks
 * for the tag ID: and compares against the specified ID.
 *
 * Returns NULL if no type/ID is found
 */
SAM_hdr_type *sam_hdr_find(SAM_hdr *hdr, char *type,
			   char *ID_key, char *ID_value);

/*!
 *
 * As per SAM_hdr_type, but returns a complete line of formatted text
 * for a specific head type/ID combination. If ID is NULL then it returns
 * the first line of the specified type.
 *
 * The returned string is malloced and should be freed by the calling
 * function with free().
 *
 * @return
 * Returns NULL if no type/ID is found.
 */
char *sam_hdr_find_line(SAM_hdr *hdr, char *type,
			char *ID_key, char *ID_value);

/*! Looks for a specific key in a single sam header line.
 *
 * If prev is non-NULL it also fills this out with the previous tag, to
 * permit use in key removal. *prev is set to NULL when the tag is the first
 * key in the list. When a tag isn't found, prev (if non NULL) will be the last
 * tag in the existing list.
 *
 * @return
 * Returns the tag pointer on success;
 *         NULL on failure
 */
SAM_hdr_tag *sam_hdr_find_key(SAM_hdr *sh,
			      SAM_hdr_type *type,
			      char *key,
			      SAM_hdr_tag **prev);

/*! Adds or updates tag key,value pairs in a header line.
 *
 * Eg for adding M5 tags to @SQ lines or updating sort order for the
 * @HD line (although use the sam_hdr_sort_order() function for
 * HD manipulation, which is a wrapper around this funuction).
 *
 * Specify multiple key,value pairs ending in NULL.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int sam_hdr_update(SAM_hdr *hdr, SAM_hdr_type *type, ...);

/*! Reconstructs the kstring from the header hash table.
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int sam_hdr_rebuild(SAM_hdr *hdr);

/*! Looks up a reference sequence by name and returns the numerical ID.
 * @return
 * Returns -1 if unknown reference.
 */
int sam_hdr_name2ref(SAM_hdr *hdr, const char *ref);

/*! Looks up a read-group by name and returns a pointer to the start of the
 * associated tag list.
 *
 * @return
 * Returns NULL on failure
 */
SAM_RG *sam_hdr_find_rg(SAM_hdr *hdr, const char *rg);

/*! Fixes any PP links in @PG headers.
 *
 * If the entries are in order then this doesn't need doing, but incase
 * our header is out of order this goes through the sh->pg[] array
 * setting the prev_id field.
 *
 * @return
 * Returns 0 on sucess;
 *        -1 on failure (indicating broken PG/PP records)
 */
int sam_hdr_link_pg(SAM_hdr *hdr);


/*! Add an @PG line.
 *
 * If we wish complete control over this use sam_hdr_add() directly. This
 * function uses that, but attempts to do a lot of tedious house work for
 * you too.
 *
 * - It will generate a suitable ID if the supplied one clashes.
 * - It will generate multiple @PG records if we have multiple PG chains.
 *
 * Call it as per sam_hdr_add() with a series of key,value pairs ending
 * in NULL.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int sam_hdr_add_PG(SAM_hdr *sh, const char *name, ...);

/*!
 * A function to help with construction of CL tags in @PG records.
 * Takes an argc, argv pair and returns a single space-separated string.
 * This string should be deallocated by the calling function.
 * 
 * @return
 * Returns malloced char * on success;
 *         NULL on failure
 */
char *stringify_argv(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif /* _SAM_HDR_H_ */
