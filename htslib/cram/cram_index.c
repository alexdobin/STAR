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

/*
 * The index is a gzipped tab-delimited text file with one line per slice.
 * The columns are:
 * 1: reference number (0 to N-1, as per BAM ref_id)
 * 2: reference position of 1st read in slice (1..?)
 * 3: number of reads in slice
 * 4: offset of container start (relative to end of SAM header, so 1st
 *    container is offset 0).
 * 5: slice number within container (ie which landmark).
 *
 * In memory, we hold this in a nested containment list. Each list element is
 * a cram_index struct. Each element in turn can contain its own list of
 * cram_index structs.
 *
 * Any start..end range which is entirely contained within another (and
 * earlier as it is sorted) range will be held within it. This ensures that
 * the outer list will never have containments and we can safely do a
 * binary search to find the first range which overlaps any given coordinate.
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
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <ctype.h>

#include "htslib/hfile.h"
#include "cram/cram.h"
#include "cram/os.h"
#include "cram/zfio.h"

#if 0
static void dump_index_(cram_index *e, int level) {
    int i, n;
    n = printf("%*s%d / %d .. %d, ", level*4, "", e->refid, e->start, e->end);
    printf("%*soffset %"PRId64"\n", MAX(0,50-n), "", e->offset);
    for (i = 0; i < e->nslice; i++) {
	dump_index_(&e->e[i], level+1);
    }
}

static void dump_index(cram_fd *fd) {
    int i;
    for (i = 0; i < fd->index_sz; i++) {
	dump_index_(&fd->index[i], 0);
    }
}
#endif

/*
 * Loads a CRAM .crai index into memory.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int cram_index_load(cram_fd *fd, const char *fn) {
    char fn2[PATH_MAX];
    char buf[65536];
    ssize_t len;
    kstring_t kstr = {0};
    hFILE *fp;
    cram_index *idx;
    cram_index **idx_stack = NULL, *ep, e;
    int idx_stack_alloc = 0, idx_stack_ptr = 0;
    size_t pos = 0;

    /* Check if already loaded */
    if (fd->index)
	return 0;

    fd->index = calloc((fd->index_sz = 1), sizeof(*fd->index));
    if (!fd->index)
	return -1;

    idx = &fd->index[0];
    idx->refid = -1;
    idx->start = INT_MIN;
    idx->end   = INT_MAX;

    idx_stack = calloc(++idx_stack_alloc, sizeof(*idx_stack));
    idx_stack[idx_stack_ptr] = idx;

    sprintf(fn2, "%s.crai", fn);
    if (!(fp = hopen(fn2, "r"))) {
	perror(fn2);
	free(idx_stack);
	return -1; 
    }

    // Load the file into memory
    while ((len = hread(fp, buf, 65536)) > 0)
	kputsn(buf, len, &kstr);
    if (len < 0 || kstr.l < 2) {
	if (kstr.s)
	    free(kstr.s);
	free(idx_stack);
	return -1;
    }

    if (hclose(fp)) {
	if (kstr.s)
	    free(kstr.s);
	free(idx_stack);
	return -1;
    }
	

    // Uncompress if required
    if (kstr.s[0] == 31 && (uc)kstr.s[1] == 139) {
	size_t l;
	char *s = zlib_mem_inflate(kstr.s, kstr.l, &l);
	free(kstr.s);
	if (!s) {
	    free(idx_stack);
	    return -1;
	}
	kstr.s = s;
	kstr.l = l;
    }


    // Parse it line at a time
    do {
	int nchars;
	char *line = &kstr.s[pos];

	/* 1.1 layout */
	if (sscanf(line, "%d\t%d\t%d\t%"PRId64"\t%d\t%d%n",
		   &e.refid,
		   &e.start,
		   &e.end,
		   &e.offset,
		   &e.slice,
		   &e.len,
		   &nchars) != 6) {
	    free(kstr.s);
	    free(idx_stack);
	    return -1;
	}

	e.end += e.start-1;
	//printf("%d/%d..%d\n", e.refid, e.start, e.end);

	if (e.refid < -1) {
	    free(kstr.s);
	    free(idx_stack);
	    fprintf(stderr, "Malformed index file, refid %d\n", e.refid);
	    return -1;
	}

	if (e.refid != idx->refid) {
	    if (fd->index_sz < e.refid+2) {
		size_t index_end = fd->index_sz * sizeof(*fd->index);
		fd->index_sz = e.refid+2;
		fd->index = realloc(fd->index,
				    fd->index_sz * sizeof(*fd->index));
		memset(((char *)fd->index) + index_end, 0,
		       fd->index_sz * sizeof(*fd->index) - index_end);
	    }
	    idx = &fd->index[e.refid+1];
	    idx->refid = e.refid;
	    idx->start = INT_MIN;
	    idx->end   = INT_MAX;
	    idx->nslice = idx->nalloc = 0;
	    idx->e = NULL;
	    idx_stack[(idx_stack_ptr = 0)] = idx;
	}

	while (!(e.start >= idx->start && e.end <= idx->end)) {
	    idx = idx_stack[--idx_stack_ptr];
	}

	// Now contains, so append
	if (idx->nslice+1 >= idx->nalloc) {
	    idx->nalloc = idx->nalloc ? idx->nalloc*2 : 16;
	    idx->e = realloc(idx->e, idx->nalloc * sizeof(*idx->e));
	}

	e.nalloc = e.nslice = 0; e.e = NULL;
	*(ep = &idx->e[idx->nslice++]) = e;
	idx = ep;

	if (++idx_stack_ptr >= idx_stack_alloc) {
	    idx_stack_alloc *= 2;
	    idx_stack = realloc(idx_stack, idx_stack_alloc*sizeof(*idx_stack));
	}
	idx_stack[idx_stack_ptr] = idx;

	pos += nchars;
	while (pos < kstr.l && kstr.s[pos] != '\n')
	    pos++;
	pos++;
    } while (pos < kstr.l);

    free(idx_stack);
    free(kstr.s);

    // dump_index(fd);

    return 0;
}

static void cram_index_free_recurse(cram_index *e) {
    if (e->e) {
	int i;
	for (i = 0; i < e->nslice; i++) {
	    cram_index_free_recurse(&e->e[i]);
	}
	free(e->e);
    }
}

void cram_index_free(cram_fd *fd) {
    int i;

    if (!fd->index)
	return;
    
    for (i = 0; i < fd->index_sz; i++) {
	cram_index_free_recurse(&fd->index[i]);
    }
    free(fd->index);

    fd->index = NULL;
}

/*
 * Searches the index for the first slice overlapping a reference ID
 * and position, or one immediately preceeding it if none is found in
 * the index to overlap this position. (Our index may have missing
 * entries, but we require at least one per reference.)
 *
 * If the index finds multiple slices overlapping this position we
 * return the first one only. Subsequent calls should specifying
 * "from" as the last slice we checked to find the next one. Otherwise
 * set "from" to be NULL to find the first one.
 *
 * Returns the cram_index pointer on sucess
 *         NULL on failure
 */
cram_index *cram_index_query(cram_fd *fd, int refid, int pos, 
			     cram_index *from) {
    int i, j, k;
    cram_index *e;

    if (refid+1 < 0 || refid+1 >= fd->index_sz)
	return NULL;

    i = 0, j = fd->index[refid+1].nslice-1;

    if (!from)
	from = &fd->index[refid+1];

    for (k = j/2; k != i; k = (j-i)/2 + i) {
	if (from->e[k].refid > refid) {
	    j = k;
	    continue;
	}

	if (from->e[k].refid < refid) {
	    i = k;
	    continue;
	}

	if (from->e[k].start >= pos) {
	    j = k;
	    continue;
	}

	if (from->e[k].start < pos) {
	    i = k;
	    continue;
	}
    }

    /* The above found *a* bin overlapping, but not necessarily the first */
    while (i > 0 && from->e[i-1].end >= pos)
	i--;

    /* Special case for matching a start pos */
    if (i+1 < from->nslice &&
	from->e[i+1].start == pos &&
	from->e[i+1].refid == refid)
	i++;

    e = &from->e[i];

    return e;
}


/*
 * Skips to a container overlapping the start coordinate listed in
 * cram_range.
 *
 * In theory we call cram_index_query multiple times, once per slice
 * overlapping the range. However slices may be absent from the index
 * which makes this problematic. Instead we find the left-most slice
 * and then read from then on, skipping decoding of slices and/or
 * whole containers when they don't overlap the specified cram_range.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek_to_refpos(cram_fd *fd, cram_range *r) {
    cram_index *e;

    // Ideally use an index, so see if we have one.
    if ((e = cram_index_query(fd, r->refid, r->start, NULL))) {
	if (0 != cram_seek(fd, e->offset, SEEK_SET))
	    if (0 != cram_seek(fd, e->offset - fd->first_container, SEEK_CUR))
		return -1;
    } else {
	fprintf(stderr, "Unknown reference ID. Missing from index?\n");
	return -1;
    }

    if (fd->ctr) {
	cram_free_container(fd->ctr);
	fd->ctr = NULL;
    }

    return 0;
}


/*
 * A specialised form of cram_index_build (below) that deals with slices
 * having multiple references in this (ref_id -2). In this scenario we
 * decode the slice to look at the RI data series instead.
 *
 * Returns 0 on success
 *        -1 on failure
 */
static int cram_index_build_multiref(cram_fd *fd,
				     cram_container *c,
				     cram_slice *s,
				     zfp *fp,
				     off_t cpos,
				     int32_t landmark,
				     int sz) {
    int i, ref = -2, ref_start = 0, ref_end;
    char buf[1024];

    if (0 != cram_decode_slice(fd, c, s, fd->header))
	return -1;

    ref_end = INT_MIN;
    for (i = 0; i < s->hdr->num_records; i++) {
	if (s->crecs[i].ref_id == ref) {
	    if (ref_end < s->crecs[i].aend)
		ref_end = s->crecs[i].aend;
	    continue;
	}

	if (ref != -2) {
	    sprintf(buf, "%d\t%d\t%d\t%"PRId64"\t%d\t%d\n",
		    ref, ref_start, ref_end - ref_start + 1,
		    (int64_t)cpos, landmark, sz);
	    zfputs(buf, fp);
	}

	ref = s->crecs[i].ref_id;
	ref_start = s->crecs[i].apos;
	ref_end = INT_MIN;
    }

    if (ref != -2) {
	sprintf(buf, "%d\t%d\t%d\t%"PRId64"\t%d\t%d\n",
		ref, ref_start, ref_end - ref_start + 1,
		(int64_t)cpos, landmark, sz);
	zfputs(buf, fp);
    }

    return 0;
}

/*
 * Builds an index file.
 *
 * fd is a newly opened cram file that we wish to index.
 * fn_base is the filename of the associated CRAM file. Internally we
 * add ".crai" to this to get the index filename.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_index_build(cram_fd *fd, const char *fn_base) {
    cram_container *c;
    off_t cpos, spos, hpos;
    zfp *fp;
    char fn_idx[PATH_MAX];

    if (strlen(fn_base) > PATH_MAX-6)
	return -1;

    sprintf(fn_idx, "%s.crai", fn_base);
    if (!(fp = zfopen(fn_idx, "wz"))) {
        perror(fn_idx);
        return -1;
    }

    cpos = htell(fd->fp);
    while ((c = cram_read_container(fd))) {
        int j;

        if (fd->err) {
            perror("Cram container read");
            return 1;
        }

        hpos = htell(fd->fp);

        if (!(c->comp_hdr_block = cram_read_block(fd)))
            return 1;
        assert(c->comp_hdr_block->content_type == COMPRESSION_HEADER);

        c->comp_hdr = cram_decode_compression_header(fd, c->comp_hdr_block);
        if (!c->comp_hdr)
            return -1;

        // 2.0 format
        for (j = 0; j < c->num_landmarks; j++) {
            char buf[1024];
            cram_slice *s;
            int sz;

            spos = htell(fd->fp);
            assert(spos - cpos - c->offset == c->landmark[j]);

            if (!(s = cram_read_slice(fd))) {
		zfclose(fp);
		return -1;
	    }

            sz = (int)(htell(fd->fp) - spos);

	    if (s->hdr->ref_seq_id == -2) {
		cram_index_build_multiref(fd, c, s, fp,
					  cpos, c->landmark[j], sz);
	    } else {
		sprintf(buf, "%d\t%d\t%d\t%"PRId64"\t%d\t%d\n",
			s->hdr->ref_seq_id, s->hdr->ref_seq_start,
			s->hdr->ref_seq_span, (int64_t)cpos,
			c->landmark[j], sz);
		zfputs(buf, fp);
	    }

            cram_free_slice(s);
        }

        cpos = htell(fd->fp);
        assert(cpos == hpos + c->length);

        cram_free_container(c);
    }
    if (fd->err) {
	zfclose(fp);
	return -1;
    }
	

    return zfclose(fp);
}
