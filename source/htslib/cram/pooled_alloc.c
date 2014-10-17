/*
Copyright (c) 2009 Genome Research Ltd.
Author: Rob Davies <rmd@sanger.ac.uk>

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

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "cram/pooled_alloc.h"

//#define TEST_MAIN

#define PSIZE 1024*1024

pool_alloc_t *pool_create(size_t dsize) {
    pool_alloc_t *p;

    if (NULL == (p = (pool_alloc_t *)malloc(sizeof(*p))))
	return NULL;

    /* Minimum size is a pointer, for free list */
    dsize = (dsize + sizeof(void *) - 1) & ~(sizeof(void *)-1);
    if (dsize < sizeof(void *))
	dsize = sizeof(void *);
    p->dsize = dsize;

    p->npools = 0;
    p->pools = NULL;
    p->free  = NULL;

    return p;
}

static pool_t *new_pool(pool_alloc_t *p) {
    size_t n = PSIZE / p->dsize;
    pool_t *pool;
    
    pool = realloc(p->pools, (p->npools + 1) * sizeof(*p->pools));
    if (NULL == pool) return NULL;
    p->pools = pool;
    pool = &p->pools[p->npools];

    pool->pool = malloc(n * p->dsize);
    if (NULL == pool->pool) return NULL;

    pool->used = 0;

    p->npools++;

    return pool;
}

void pool_destroy(pool_alloc_t *p) {
    size_t i;

    for (i = 0; i < p->npools; i++) {
        free(p->pools[i].pool);
    }
    free(p->pools);
    free(p);
}

void *pool_alloc(pool_alloc_t *p) {
    pool_t *pool;
    void *ret;

    /* Look on free list */
    if (NULL != p->free) {
        ret = p->free;
	p->free = *((void **)p->free);
	return ret;
    }

    /* Look for space in the last pool */
    if (p->npools) {
        pool = &p->pools[p->npools - 1];
        if (pool->used + p->dsize < PSIZE) {
	    ret = ((char *) pool->pool) + pool->used;
	    pool->used += p->dsize;
	    return ret;
	}
    }

    /* Need a new pool */
    pool = new_pool(p);
    if (NULL == pool) return NULL;

    pool->used = p->dsize;
    return pool->pool;
}

void pool_free(pool_alloc_t *p, void *ptr) {
    *(void **)ptr = p->free;
    p->free = ptr;
}

#ifdef TEST_MAIN
typedef struct {
    int x, y, z;
} xyz;

#define NP 10000
int main(void) {
    int i;
    xyz *item;
    xyz **items;
    pool_alloc_t *p = pool_create(sizeof(xyz));

    items = (xyz **)malloc(NP * sizeof(*items));

    for (i = 0; i < NP; i++) {
	item = pool_alloc(p);
	item->x = i;
	item->y = i+1;
	item->z = i+2;
	items[i] = item;
    }

    for (i = 0; i < NP; i++) {
	item = items[i];
	if (i % 3)
	    pool_free(p, item);
    }

    for (i = 0; i < NP; i++) {
	item = pool_alloc(p);
	item->x = 1000000+i;
	item->y = 1000000+i+1;
	item->z = 1000000+i+2;
    }

    for (i = 0; i < NP; i++) {
	item = items[i];
	printf("%d\t%d\t%d\t%d\n", i, item->x, item->y, item->z);
	pool_free(p, item);
    }

    return 0;
}
#endif
