/*
Copyright (c) 2013 Genome Research Ltd.
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

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <string.h>
#include <assert.h>

#include "cram/sam_header.h"
#include "cram/string_alloc.h"

#ifdef SAMTOOLS
#define sam_hdr_parse sam_hdr_parse_
#endif

static void sam_hdr_error(char *msg, char *line, int len, int lno) {
    int j;
    
    for (j = 0; j < len && line[j] != '\n'; j++)
	;
    fprintf(stderr, "%s at line %d: \"%.*s\"\n", msg, lno, j, line);
}

void sam_hdr_dump(SAM_hdr *hdr) {
    khint_t k;
    int i;

    printf("===DUMP===\n");
    for (k = kh_begin(hdr->h); k != kh_end(hdr->h); k++) {
	SAM_hdr_type *t1, *t2;
	char c[2];

	if (!kh_exist(hdr->h, k))
	    continue;

	t1 = t2 = kh_val(hdr->h, k);
	c[0] = kh_key(hdr->h, k)>>8;
	c[1] = kh_key(hdr->h, k)&0xff;
	printf("Type %.2s, count %d\n", c, t1->prev->order+1);

	do {
	    SAM_hdr_tag *tag;
	    printf(">>>%d ", t1->order);
	    for (tag = t1->tag; tag; tag=tag->next) {
		printf("\"%.2s\":\"%.*s\"\t",
		       tag->str, tag->len-3, tag->str+3);
	    }
	    putchar('\n');
	    t1 = t1->next;
	} while (t1 != t2);
    }

    /* Dump out PG chains */
    printf("\n@PG chains:\n");
    for (i = 0; i < hdr->npg_end; i++) {
	int j;
	printf("  %d:", i);
	for (j = hdr->pg_end[i]; j != -1; j = hdr->pg[j].prev_id) {
	    printf("%s%d(%.*s)", 
		   j == hdr->pg_end[i] ? " " : "->",
		   j, hdr->pg[j].name_len, hdr->pg[j].name);
	}
	printf("\n");
    }

    puts("===END DUMP===");
}

/* Updates the hash tables in the SAM_hdr structure.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
static int sam_hdr_update_hashes(SAM_hdr *sh,
				 int type,
				 SAM_hdr_type *h_type) {
    /* Add to reference hash? */
    if ((type>>8) == 'S' && (type&0xff) == 'Q') {
	SAM_hdr_tag *tag;
	int nref = sh->nref;

	sh->ref = realloc(sh->ref, (sh->nref+1)*sizeof(*sh->ref));
	if (!sh->ref)
	    return -1;

	tag = h_type->tag;
	sh->ref[nref].name = NULL;
	sh->ref[nref].len  = 0;
	sh->ref[nref].ty = h_type;
	sh->ref[nref].tag  = tag;

	while (tag) {
	    if (tag->str[0] == 'S' && tag->str[1] == 'N') {
		if (!(sh->ref[nref].name = malloc(tag->len)))
		    return -1;
		strncpy(sh->ref[nref].name, tag->str+3, tag->len-3);
		sh->ref[nref].name[tag->len-3] = 0;
	    } else if (tag->str[0] == 'L' && tag->str[1] == 'N') {
		sh->ref[nref].len = atoi(tag->str+3);
	    }
	    tag = tag->next;
	}

	if (sh->ref[nref].name) {
	    khint_t k;
	    int r;
	    k = kh_put(m_s2i, sh->ref_hash, sh->ref[nref].name, &r);
	    if (-1 == r) return -1;
	    kh_val(sh->ref_hash, k) = nref;
	}

	sh->nref++;
    }

    /* Add to read-group hash? */
    if ((type>>8) == 'R' && (type&0xff) == 'G') {
	SAM_hdr_tag *tag;
	int nrg = sh->nrg;

	sh->rg = realloc(sh->rg, (sh->nrg+1)*sizeof(*sh->rg));
	if (!sh->rg)
	    return -1;

	tag = h_type->tag;
	sh->rg[nrg].name = NULL;
	sh->rg[nrg].name_len = 0;
	sh->rg[nrg].ty   = h_type;
	sh->rg[nrg].tag  = tag;
	sh->rg[nrg].id   = nrg;

	while (tag) {
	    if (tag->str[0] == 'I' && tag->str[1] == 'D') {
		if (!(sh->rg[nrg].name = malloc(tag->len)))
		    return -1;
		strncpy(sh->rg[nrg].name, tag->str+3, tag->len-3);
		sh->rg[nrg].name[tag->len-3] = 0;
		sh->rg[nrg].name_len = strlen(sh->rg[nrg].name);
	    }
	    tag = tag->next;
	}

	if (sh->rg[nrg].name) {
	    khint_t k;
	    int r;
	    k = kh_put(m_s2i, sh->rg_hash, sh->rg[nrg].name, &r);
	    if (-1 == r) return -1;
	    kh_val(sh->rg_hash, k) = nrg;
	}

	sh->nrg++;
    }

    /* Add to program hash? */
    if ((type>>8) == 'P' && (type&0xff) == 'G') {
	SAM_hdr_tag *tag;
	int npg = sh->npg;

	sh->pg = realloc(sh->pg, (sh->npg+1)*sizeof(*sh->pg));
	if (!sh->pg)
	    return -1;

	tag = h_type->tag;
	sh->pg[npg].name = NULL;
	sh->pg[npg].name_len = 0;
	sh->pg[npg].ty  = h_type;
	sh->pg[npg].tag  = tag;
	sh->pg[npg].id   = npg;
	sh->pg[npg].prev_id = -1;

	while (tag) {
	    if (tag->str[0] == 'I' && tag->str[1] == 'D') {
		if (!(sh->pg[npg].name = malloc(tag->len)))
		    return -1;
		strncpy(sh->pg[npg].name, tag->str+3, tag->len-3);
		sh->pg[npg].name[tag->len-3] = 0;
		sh->pg[npg].name_len = strlen(sh->pg[npg].name);
	    } else if (tag->str[0] == 'P' && tag->str[1] == 'P') {
		// Resolve later if needed
		khint_t k;
		char tmp = tag->str[tag->len]; tag->str[tag->len] = 0;
		k = kh_get(m_s2i, sh->pg_hash, tag->str+3);
		tag->str[tag->len] = tmp;

		if (k != kh_end(sh->pg_hash)) {
		    int p_id = kh_val(sh->pg_hash, k);
		    sh->pg[npg].prev_id = sh->pg[p_id].id;

		    /* Unmark previous entry as a PG termination */
		    if (sh->npg_end > 0 &&
			sh->pg_end[sh->npg_end-1] == p_id) {
			sh->npg_end--;
		    } else {
			int i;
			for (i = 0; i < sh->npg_end; i++) {
			    if (sh->pg_end[i] == p_id) {
				memmove(&sh->pg_end[i], &sh->pg_end[i+1],
					(sh->npg_end-i-1)*sizeof(*sh->pg_end));
				sh->npg_end--;
			    }
			}
		    }
		} else {
		    sh->pg[npg].prev_id = -1;
		}
	    }
	    tag = tag->next;
	}

	if (sh->pg[npg].name) {
	    khint_t k;
	    int r;
	    k = kh_put(m_s2i, sh->pg_hash, sh->pg[npg].name, &r);
	    if (-1 == r) return -1;
	    kh_val(sh->pg_hash, k) = npg;
	}

	/* Add to npg_end[] array. Remove later if we find a PP line */
	if (sh->npg_end >= sh->npg_end_alloc) {
	    sh->npg_end_alloc = sh->npg_end_alloc
		? sh->npg_end_alloc*2
		: 4;
	    sh->pg_end = realloc(sh->pg_end,
				 sh->npg_end_alloc * sizeof(int));
	    if (!sh->pg_end)
		return -1;
	}
	sh->pg_end[sh->npg_end++] = npg;

	sh->npg++;
    }

    return 0;
}

/*
 * Appends a formatted line to an existing SAM header.
 * Line is a full SAM header record, eg "@SQ\tSN:foo\tLN:100", with
 * optional new-line. If it contains more than 1 line then multiple lines
 * will be added in order.
 *
 * Len is the length of the text data, or 0 if unknown (in which case
 * it should be null terminated).
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sam_hdr_add_lines(SAM_hdr *sh, const char *lines, int len) {
    int i, lno = 1, text_offset;
    char *hdr;

    if (!len)
	len = strlen(lines);

    text_offset = ks_len(&sh->text);
    if (EOF == kputsn(lines, len, &sh->text))
	return -1;
    hdr = ks_str(&sh->text) + text_offset;

    for (i = 0; i < len; i++) {
	khint32_t type;
	khint_t k;

	int l_start = i, new;
	SAM_hdr_type *h_type;
	SAM_hdr_tag *h_tag, *last;

	if (hdr[i] != '@') {
	    int j;
	    for (j = i; j < len && hdr[j] != '\n'; j++)
		;
	    sam_hdr_error("Header line does not start with '@'",
			  &hdr[l_start], len - l_start, lno);
	    return -1;
	}

	type = (hdr[i+1]<<8) | hdr[i+2];
	if (hdr[i+1] < 'A' || hdr[i+1] > 'z' ||
	    hdr[i+2] < 'A' || hdr[i+2] > 'z') {
	    sam_hdr_error("Header line does not have a two character key",
			  &hdr[l_start], len - l_start, lno);
	    return -1;
	}

	i += 3;
	if (hdr[i] == '\n')
	    continue;

	// Add the header line type
	if (!(h_type = pool_alloc(sh->type_pool)))
	    return -1;
	if (-1 == (k = kh_put(sam_hdr, sh->h, type, &new)))
	    return -1;

	// Form the ring, either with self or other lines of this type
	if (!new) {
	    SAM_hdr_type *t = kh_val(sh->h, k), *p;
	    p = t->prev;
	    
	    assert(p->next = t);
	    p->next = h_type;
	    h_type->prev = p;

	    t->prev = h_type;
	    h_type->next = t;
	    h_type->order = p->order+1;
	} else {
	    kh_val(sh->h, k) = h_type;
	    h_type->prev = h_type->next = h_type;
	    h_type->order = 0;
	}

	// Parse the tags on this line
	last = NULL;
	if ((type>>8) == 'C' && (type&0xff) == 'O') {
	    int j;
	    if (hdr[i] != '\t') {
		sam_hdr_error("Missing tab",
			      &hdr[l_start], len - l_start, lno);
		return -1;
	    }

	    for (j = ++i; j < len && hdr[j] != '\n'; j++)
		;

	    if (!(h_type->tag = h_tag = pool_alloc(sh->tag_pool)))
		return -1;
	    h_tag->str = string_ndup(sh->str_pool, &hdr[i], j-i);
	    h_tag->len = j-i;
	    h_tag->next = NULL;
	    if (!h_tag->str)
		return -1;

	    i = j;

	} else {
	    do {
		int j;
		if (hdr[i] != '\t') {
		    sam_hdr_error("Missing tab",
				  &hdr[l_start], len - l_start, lno);
		    return -1;
		}

		for (j = ++i; j < len && hdr[j] != '\n' && hdr[j] != '\t'; j++)
		    ;
	    
		if (!(h_tag = pool_alloc(sh->tag_pool)))
		    return -1;
		h_tag->str = string_ndup(sh->str_pool, &hdr[i], j-i);
		h_tag->len = j-i;
		h_tag->next = NULL;
		if (!h_tag->str)
		    return -1;

		if (h_tag->len < 3 || h_tag->str[2] != ':') {
		    sam_hdr_error("Malformed key:value pair",
				  &hdr[l_start], len - l_start, lno);
		    return -1;
		}
	    
		if (last)
		    last->next = h_tag;
		else
		    h_type->tag = h_tag;

		last = h_tag;
		i = j;
	    } while (i < len && hdr[i] != '\n');
	}

	/* Update RG/SQ hashes */
	if (-1 == sam_hdr_update_hashes(sh, type, h_type))
	    return -1;
    }

    return 0;
}

/*
 * Adds a single line to a SAM header.
 * Specify type and one or more key,value pairs, ending with the NULL key.
 * Eg. sam_hdr_add(h, "SQ", "ID", "foo", "LN", "100", NULL).
 *
 * Returns index for specific entry on success (eg 2nd SQ, 4th RG)
 *        -1 on failure
 */
int sam_hdr_add(SAM_hdr *sh, const char *type, ...) {
    va_list args;
    va_start(args, type);
    return sam_hdr_vadd(sh, type, args, NULL);
}

int sam_hdr_vadd(SAM_hdr *sh, const char *type, va_list ap, ...) {
    va_list args;
    SAM_hdr_type *h_type;
    SAM_hdr_tag *h_tag, *last;
    int new;
    khint32_t type_i = (type[0]<<8) | type[1], k;

#if defined(HAVE_VA_COPY)
    va_list ap_local;
#endif

    if (EOF == kputc_('@', &sh->text))
	return -1;
    if (EOF == kputsn(type, 2, &sh->text))
	return -1;

    if (!(h_type = pool_alloc(sh->type_pool)))
	return -1;
    if (-1 == (k = kh_put(sam_hdr, sh->h, type_i, &new)))
	return -1;
    kh_val(sh->h, k) = h_type;

    // Form the ring, either with self or other lines of this type
    if (!new) {
	SAM_hdr_type *t = kh_val(sh->h, k), *p;
	p = t->prev;
	    
	assert(p->next = t);
	p->next = h_type;
	h_type->prev = p;

	t->prev = h_type;
	h_type->next = t;
	h_type->order = p->order + 1;
    } else {
	h_type->prev = h_type->next = h_type;
	h_type->order = 0;
    }

    last = NULL;

    // Any ... varargs
    va_start(args, ap);
    for (;;) {
	char *k, *v;
	int idx;
	
	if (!(k = (char *)va_arg(args, char *)))
	    break;
	v = va_arg(args, char *);

	if (EOF == kputc_('\t', &sh->text))
	    return -1;

	if (!(h_tag = pool_alloc(sh->tag_pool)))
	    return -1;
	idx = ks_len(&sh->text);
	
	if (EOF == kputs(k, &sh->text))
	    return -1;
	if (EOF == kputc_(':', &sh->text))
	    return -1;
	if (EOF == kputs(v, &sh->text))
	    return -1;

	h_tag->len = ks_len(&sh->text) - idx;
	h_tag->str = string_ndup(sh->str_pool,
				 ks_str(&sh->text) + idx,
				 h_tag->len);
	h_tag->next = NULL;
	if (!h_tag->str)
	    return -1;

	if (last)
	    last->next = h_tag;
	else
	    h_type->tag = h_tag;
	
	last = h_tag;
    }
    va_end(args);

#if defined(HAVE_VA_COPY)
    va_copy(ap_local, ap);
#   define ap ap_local
#endif

    // Plus the specified va_list params
    for (;;) {
	char *k, *v;
	int idx;
	
	if (!(k = (char *)va_arg(ap, char *)))
	    break;
	v = va_arg(ap, char *);

	if (EOF == kputc_('\t', &sh->text))
	    return -1;

	if (!(h_tag = pool_alloc(sh->tag_pool)))
	    return -1;
	idx = ks_len(&sh->text);
	
	if (EOF == kputs(k, &sh->text))
	    return -1;
	if (EOF == kputc_(':', &sh->text))
	    return -1;
	if (EOF == kputs(v, &sh->text))
	    return -1;

	h_tag->len = ks_len(&sh->text) - idx;
	h_tag->str = string_ndup(sh->str_pool,
				 ks_str(&sh->text) + idx,
				 h_tag->len);
	h_tag->next = NULL;
	if (!h_tag->str)
	    return -1;

	if (last)
	    last->next = h_tag;
	else
	    h_type->tag = h_tag;
	
	last = h_tag;
    }
    va_end(ap);

    if (EOF == kputc('\n', &sh->text))
	return -1;

    int itype = (type[0]<<8) | type[1];
    if (-1 == sam_hdr_update_hashes(sh, itype, h_type))
	return -1;

    return h_type->order;
}

/*
 * Returns the first header item matching 'type'. If ID is non-NULL it checks
 * for the tag ID: and compares against the specified ID.
 *
 * Returns NULL if no type/ID is found
 */
SAM_hdr_type *sam_hdr_find(SAM_hdr *hdr, char *type,
			   char *ID_key, char *ID_value) {
    SAM_hdr_type *t1, *t2;
    int itype = (type[0]<<8)|(type[1]);
    khint_t k;

    /* Special case for types we have prebuilt hashes on */
    if (ID_key) {
	if (type[0]   == 'S' && type[1]   == 'Q' &&
	    ID_key[0] == 'S' && ID_key[1] == 'N') {
	    k = kh_get(m_s2i, hdr->ref_hash, ID_value);
	    return k != kh_end(hdr->ref_hash)
		? hdr->ref[kh_val(hdr->ref_hash, k)].ty
		: NULL;
	}

	if (type[0]   == 'R' && type[1]   == 'G' &&
	    ID_key[0] == 'I' && ID_key[1] == 'D') {
	    k = kh_get(m_s2i, hdr->rg_hash, ID_value);
	    return k != kh_end(hdr->rg_hash)
		? hdr->rg[kh_val(hdr->rg_hash, k)].ty
		: NULL;
	}

	if (type[0]   == 'P' && type[1]   == 'G' &&
	    ID_key[0] == 'I' && ID_key[1] == 'D') {
	    k = kh_get(m_s2i, hdr->pg_hash, ID_value);
	    return k != kh_end(hdr->pg_hash)
		? hdr->pg[kh_val(hdr->pg_hash, k)].ty
		: NULL;
	}
    }

    k = kh_get(sam_hdr, hdr->h, itype);
    if (k == kh_end(hdr->h))
	return NULL;
    
    if (!ID_key)
	return kh_val(hdr->h, k);

    t1 = t2 = kh_val(hdr->h, k);
    do {
	SAM_hdr_tag *tag;
	for (tag = t1->tag; tag; tag = tag->next) {
	    if (tag->str[0] == ID_key[0] && tag->str[1] == ID_key[1]) {
		char *cp1 = tag->str+3;
		char *cp2 = ID_value;
		while (*cp1 && *cp1 == *cp2)
		    cp1++, cp2++;
		if (*cp2 || *cp1)
		    continue;
		return t1;
	    }
	}
	t1 = t1->next;
    } while (t1 != t2);

    return NULL;
}

/*
 * As per SAM_hdr_type, but returns a complete line of formatted text
 * for a specific head type/ID combination. If ID is NULL then it returns
 * the first line of the specified type.
 *
 * The returned string is malloced and should be freed by the calling
 * function with free().
 *
 * Returns NULL if no type/ID is found.
 */
char *sam_hdr_find_line(SAM_hdr *hdr, char *type,
			char *ID_key, char *ID_value) {
    SAM_hdr_type *ty = sam_hdr_find(hdr, type, ID_key, ID_value);
    kstring_t ks = KS_INITIALIZER;
    SAM_hdr_tag *tag;
    int r = 0;

    if (!ty)
	return NULL;

    // Paste together the line from the hashed copy
    r |= (kputc_('@', &ks) == EOF);
    r |= (kputs(type, &ks) == EOF);
    for (tag = ty->tag; tag; tag = tag->next) {
	r |= (kputc_('\t', &ks) == EOF);
	r |= (kputsn(tag->str, tag->len, &ks) == EOF);
    }

    if (r) {
	KS_FREE(&ks);
	return NULL;
    }

    return ks_str(&ks);
}


/*
 * Looks for a specific key in a single sam header line.
 * If prev is non-NULL it also fills this out with the previous tag, to
 * permit use in key removal. *prev is set to NULL when the tag is the first
 * key in the list. When a tag isn't found, prev (if non NULL) will be the last
 * tag in the existing list.
 *
 * Returns the tag pointer on success
 *         NULL on failure
 */
SAM_hdr_tag *sam_hdr_find_key(SAM_hdr *sh,
			      SAM_hdr_type *type,
			      char *key,
			      SAM_hdr_tag **prev) {
    SAM_hdr_tag *tag, *p = NULL;

    for (tag = type->tag; tag; p = tag, tag = tag->next) {
	if (tag->str[0] == key[0] && tag->str[1] == key[1]) {
	    if (prev)
		*prev = p;
	    return tag;
	}
    }

    if (prev)
	*prev = p;

    return NULL;
}


/*
 * Adds or updates tag key,value pairs in a header line.
 * Eg for adding M5 tags to @SQ lines or updating sort order for the
 * @HD line (although use the sam_hdr_sort_order() function for
 * HD manipulation, which is a wrapper around this funuction).
 *
 * Specify multiple key,value pairs ending in NULL.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int sam_hdr_update(SAM_hdr *hdr, SAM_hdr_type *type, ...) {
    va_list ap;

    va_start(ap, type);
    
    for (;;) {
	char *k, *v;
	int idx;
	SAM_hdr_tag *tag, *prev;

	if (!(k = (char *)va_arg(ap, char *)))
	    break;
	v = va_arg(ap, char *);

	tag = sam_hdr_find_key(hdr, type, k, &prev);
	if (!tag) {
	    if (!(tag = pool_alloc(hdr->tag_pool)))
		return -1;
	    if (prev)
		prev->next = tag;
	    else
		type->tag = tag;

	    tag->next = NULL;
	}

	idx = ks_len(&hdr->text);
	if (ksprintf(&hdr->text, "%2.2s:%s", k, v) < 0)
	    return -1;
	tag->len = ks_len(&hdr->text) - idx;
	tag->str = string_ndup(hdr->str_pool,
			       ks_str(&hdr->text) + idx,
			       tag->len);
	if (!tag->str)
	    return -1;
    }

    va_end(ap);

    return 0;
}

#define K(a) (((a)[0]<<8)|((a)[1]))

/*
 * Reconstructs the kstring from the header hash table.
 * Returns 0 on success
 *        -1 on failure
 */
int sam_hdr_rebuild(SAM_hdr *hdr) {
    /* Order: HD then others */
    kstring_t ks = KS_INITIALIZER;
    khint_t k;


    k = kh_get(sam_hdr, hdr->h, K("HD"));
    if (k != kh_end(hdr->h)) {
	SAM_hdr_type *ty = kh_val(hdr->h, k);
	SAM_hdr_tag *tag;
	if (EOF == kputs("@HD", &ks))
	    return -1;
	for (tag = ty->tag; tag; tag = tag->next) {
	    if (EOF == kputc_('\t', &ks))
		return -1;
	    if (EOF == kputsn_(tag->str, tag->len, &ks))
		return -1;
	}
	if (EOF == kputc('\n', &ks))
	    return -1;
    }

    for (k = kh_begin(hdr->h); k != kh_end(hdr->h); k++) {
	SAM_hdr_type *t1, *t2;

	if (!kh_exist(hdr->h, k))
	    continue;

	if (kh_key(hdr->h, k) == K("HD"))
	    continue;

	t1 = t2 = kh_val(hdr->h, k);
	do {
	    SAM_hdr_tag *tag;
	    char c[2];

	    if (EOF == kputc_('@', &ks))
		return -1;
	    c[0] = kh_key(hdr->h, k)>>8;
	    c[1] = kh_key(hdr->h, k)&0xff;
	    if (EOF == kputsn_(c, 2, &ks))
		return -1;
	    for (tag = t1->tag; tag; tag=tag->next) {
		if (EOF == kputc_('\t', &ks))
		    return -1;
		if (EOF == kputsn_(tag->str, tag->len, &ks))
		    return -1;
	    }
	    if (EOF == kputc('\n', &ks))
		return -1;
	    t1 = t1->next;
	} while (t1 != t2);
    }

    if (ks_str(&hdr->text))
	KS_FREE(&hdr->text);

    hdr->text = ks;

    return 0;
}


/*
 * Creates an empty SAM header, ready to be populated.
 * 
 * Returns a SAM_hdr struct on success (free with sam_hdr_free())
 *         NULL on failure
 */
SAM_hdr *sam_hdr_new() {
    SAM_hdr *sh = calloc(1, sizeof(*sh));

    if (!sh)
	return NULL;
    
    sh->h = kh_init(sam_hdr);
    if (!sh->h)
	goto err;

    sh->ID_cnt = 1;
    sh->ref_count = 1;

    sh->nref = 0;
    sh->ref  = NULL;
    if (!(sh->ref_hash = kh_init(m_s2i)))
	goto err;

    sh->nrg = 0;
    sh->rg  = NULL;
    if (!(sh->rg_hash = kh_init(m_s2i)))
	goto err;

    sh->npg = 0;
    sh->pg  = NULL;
    sh->npg_end = sh->npg_end_alloc = 0;
    sh->pg_end = NULL;
    if (!(sh->pg_hash = kh_init(m_s2i)))
	goto err;

    KS_INIT(&sh->text);

    if (!(sh->tag_pool = pool_create(sizeof(SAM_hdr_tag))))
	goto err;

    if (!(sh->type_pool = pool_create(sizeof(SAM_hdr_type))))
	goto err;

    if (!(sh->str_pool = string_pool_create(8192)))
	goto err;

    return sh;

 err:
    if (sh->h)
	kh_destroy(sam_hdr, sh->h);

    if (sh->tag_pool)
	pool_destroy(sh->tag_pool);

    if (sh->type_pool)
	pool_destroy(sh->type_pool);

    if (sh->str_pool)
	string_pool_destroy(sh->str_pool);

    free(sh);

    return NULL;
}


/*
 * Tokenises a SAM header into a hash table.
 * Also extracts a few bits on specific data types, such as @RG lines.
 *
 * Returns a SAM_hdr struct on success (free with sam_hdr_free())
 *         NULL on failure
 */
SAM_hdr *sam_hdr_parse(const char *hdr, int len) {
    /* Make an empty SAM_hdr */
    SAM_hdr *sh;
    
    sh = sam_hdr_new();
    if (NULL == sh) return NULL;

    if (NULL == hdr) return sh; // empty header is permitted

    /* Parse the header, line by line */
    if (-1 == sam_hdr_add_lines(sh, hdr, len)) {
	sam_hdr_free(sh);
	return NULL;
    }

    //sam_hdr_dump(sh);
    //sam_hdr_add(sh, "RG", "ID", "foo", "SM", "bar", NULL);
    //sam_hdr_rebuild(sh);
    //printf(">>%s<<", ks_str(sh->text));

    //parse_references(sh);
    //parse_read_groups(sh);

    sam_hdr_link_pg(sh);
    //sam_hdr_dump(sh);

    return sh;
}

/*
 * Produces a duplicate copy of hdr and returns it.
 * Returns NULL on failure
 */
SAM_hdr *sam_hdr_dup(SAM_hdr *hdr) {
    if (-1 == sam_hdr_rebuild(hdr))
	return NULL;

    return sam_hdr_parse(sam_hdr_str(hdr), sam_hdr_length(hdr));
}

/*! Increments a reference count on hdr.
 *
 * This permits multiple files to share the same header, all calling
 * sam_hdr_free when done, without causing errors for other open  files.
 */
void sam_hdr_incr_ref(SAM_hdr *hdr) {
    hdr->ref_count++;
}

/*! Increments a reference count on hdr.
 *
 * This permits multiple files to share the same header, all calling
 * sam_hdr_free when done, without causing errors for other open  files.
 *
 * If the reference count hits zero then the header is automatically
 * freed. This makes it a synonym for sam_hdr_free().
 */
void sam_hdr_decr_ref(SAM_hdr *hdr) {
    sam_hdr_free(hdr);
}

/*! Deallocates all storage used by a SAM_hdr struct.
 *
 * This also decrements the header reference count. If after decrementing 
 * it is still non-zero then the header is assumed to be in use by another
 * caller and the free is not done.
 *
 * This is a synonym for sam_hdr_dec_ref().
 */
void sam_hdr_free(SAM_hdr *hdr) {
    if (!hdr)
	return;

    if (--hdr->ref_count > 0)
	return;

    if (ks_str(&hdr->text))
	KS_FREE(&hdr->text);

    if (hdr->h)
	kh_destroy(sam_hdr, hdr->h);

    if (hdr->ref_hash)
	kh_destroy(m_s2i, hdr->ref_hash);

    if (hdr->ref) {
	int i;
	for (i = 0; i < hdr->nref; i++)
	    if (hdr->ref[i].name)
		free(hdr->ref[i].name);
	free(hdr->ref);
    }

    if (hdr->rg_hash)
	kh_destroy(m_s2i, hdr->rg_hash);

    if (hdr->rg) {
	int i;
	for (i = 0; i < hdr->nrg; i++)
	    if (hdr->rg[i].name)
		free(hdr->rg[i].name);
	free(hdr->rg);
    }

    if (hdr->pg_hash)
	kh_destroy(m_s2i, hdr->pg_hash);

    if (hdr->pg) {
	int i;
	for (i = 0; i < hdr->npg; i++)
	    if (hdr->pg[i].name)
		free(hdr->pg[i].name);
	free(hdr->pg);
    }

    if (hdr->pg_end)
	free(hdr->pg_end);

    if (hdr->type_pool)
	pool_destroy(hdr->type_pool);

    if (hdr->tag_pool)
	pool_destroy(hdr->tag_pool);

    if (hdr->str_pool)
	string_pool_destroy(hdr->str_pool);

    free(hdr);
}

int sam_hdr_length(SAM_hdr *hdr) {
    return ks_len(&hdr->text);
}

char *sam_hdr_str(SAM_hdr *hdr) {
    return ks_str(&hdr->text);
}

/*
 * Looks up a reference sequence by name and returns the numerical ID.
 * Returns -1 if unknown reference.
 */
int sam_hdr_name2ref(SAM_hdr *hdr, const char *ref) {
    khint_t k = kh_get(m_s2i, hdr->ref_hash, ref);
    return k == kh_end(hdr->ref_hash) ? -1 : kh_val(hdr->ref_hash, k);
}

/*
 * Looks up a read-group by name and returns a pointer to the start of the
 * associated tag list.
 *
 * Returns NULL on failure
 */
SAM_RG *sam_hdr_find_rg(SAM_hdr *hdr, const char *rg) {
    khint_t k = kh_get(m_s2i, hdr->rg_hash, rg);
    return k == kh_end(hdr->rg_hash)
	? NULL
	: &hdr->rg[kh_val(hdr->rg_hash, k)];
}


/*
 * Fixes any PP links in @PG headers.
 * If the entries are in order then this doesn't need doing, but incase
 * our header is out of order this goes through the sh->pg[] array
 * setting the prev_id field.
 *
 * Note we can have multiple complete chains. This code should identify the
 * tails of these chains as these are the entries we have to link to in
 * subsequent PP records.
 *
 * Returns 0 on sucess
 *        -1 on failure (indicating broken PG/PP records)
 */
int sam_hdr_link_pg(SAM_hdr *hdr) {
    int i, j, ret = 0;

    hdr->npg_end_alloc = hdr->npg;
    hdr->pg_end = realloc(hdr->pg_end, hdr->npg * sizeof(*hdr->pg_end));
    if (!hdr->pg_end)
	return -1;

    for (i = 0; i < hdr->npg; i++)
	hdr->pg_end[i] = i;

    for (i = 0; i < hdr->npg; i++) {
	khint_t k;
	SAM_hdr_tag *tag;
	char tmp;

	for (tag = hdr->pg[i].tag; tag; tag = tag->next) {
	    if (tag->str[0] == 'P' && tag->str[1] == 'P')
		break;
	}
	if (!tag) {
	    /* Chain start points */
	    continue;
	}

	tmp = tag->str[tag->len]; tag->str[tag->len] = 0;
	k = kh_get(m_s2i, hdr->pg_hash, tag->str+3);
	tag->str[tag->len] = tmp;

	if (k == kh_end(hdr->pg_hash)) {
	    ret = -1;
	    continue;
	}

	hdr->pg[i].prev_id = hdr->pg[kh_val(hdr->pg_hash, k)].id;
	hdr->pg_end[kh_val(hdr->pg_hash, k)] = -1;
    }

    for (i = j = 0; i < hdr->npg; i++) {
	if (hdr->pg_end[i] != -1)
	    hdr->pg_end[j++] = hdr->pg_end[i];
    }
    hdr->npg_end = j;

    return ret;
}

/*
 * Returns a unique ID from a base name.
 *
 * The value returned is valid until the next call to
 * this function.
 */
const char *sam_hdr_PG_ID(SAM_hdr *sh, const char *name) {
    khint_t k = kh_get(m_s2i, sh->pg_hash, name);
    if (k == kh_end(sh->pg_hash))
	return name;

    do {
	sprintf(sh->ID_buf, "%.1000s.%d", name, sh->ID_cnt++);
	k = kh_get(m_s2i, sh->pg_hash, sh->ID_buf);
    } while (k == kh_end(sh->pg_hash));

    return sh->ID_buf;
}

/*
 * Add an @PG line.
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
 * Returns 0 on success
 *        -1 on failure
 */
int sam_hdr_add_PG(SAM_hdr *sh, const char *name, ...) {
    va_list args;
    va_start(args, name);

    if (sh->npg_end) {
	/* Copy ends array to avoid us looping while modifying it */
	int *end = malloc(sh->npg_end * sizeof(int));
	int i, nends = sh->npg_end;

	if (!end)
	    return -1;

	memcpy(end, sh->pg_end, nends * sizeof(*end));

	for (i = 0; i < nends; i++) {
	    if (-1 == sam_hdr_vadd(sh, "PG", args,
				   "ID", sam_hdr_PG_ID(sh, name),
				   "PN", name,
				   "PP", sh->pg[end[i]].name,
				   NULL)) {
		free(end);
		return  -1;
	    }
	}

	free(end);
    } else {
	if (-1 == sam_hdr_vadd(sh, "PG", args,
			       "ID", sam_hdr_PG_ID(sh, name),
			       "PN", name,
			       NULL))
	    return -1;
    }

    //sam_hdr_dump(sh);

    return 0;
}

/*
 * A function to help with construction of CL tags in @PG records.
 * Takes an argc, argv pair and returns a single space-separated string.
 * This string should be deallocated by the calling function.
 * 
 * Returns malloced char * on success
 *         NULL on failure
 */
char *stringify_argv(int argc, char *argv[]) {
    char *str, *cp;
    size_t nbytes = 1;
    int i, j;

    /* Allocate */
    for (i = 0; i < argc; i++) {
	nbytes += strlen(argv[i]) + 1;
    }
    if (!(str = malloc(nbytes)))
	return NULL;

    /* Copy */
    cp = str;
    for (i = 0; i < argc; i++) {
	j = 0;
	while (argv[i][j]) {
	    if (argv[i][j] == '\t')
		*cp++ = ' ';
	    else
		*cp++ = argv[i][j];
	    j++;
	}
	*cp++ = ' ';
    }
    *cp++ = 0;

    return str;
}
