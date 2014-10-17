/*
Copyright (c) 2010 Genome Research Ltd.
Author: Andrew Whitwham <aw7@sanger.ac.uk>

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
   A pooled string allocator intended to cut down on the
   memory overhead of many small string allocations.
   
   Andrew Whitwham, September 2010.
*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "string_alloc.h"

#define MIN_STR_SIZE 1024


/* creates the string pool. max_length is the initial size
   a single string can be.  Tha max_length can grow as
   needed */

string_alloc_t *string_pool_create(size_t max_length) {
    string_alloc_t *a_str;
    
    if (NULL == (a_str = (string_alloc_t *)malloc(sizeof(*a_str)))) {
    	return NULL;
    }
    
    if (max_length < MIN_STR_SIZE) max_length = MIN_STR_SIZE;
    
    a_str->nstrings    = 0;
    a_str->max_length  = max_length;
    a_str->strings     = NULL;
    
    return a_str;
}


/* internal function to do the actual memory allocation */

static string_t *new_string_pool(string_alloc_t *a_str) {
    string_t *str;
    
    str = realloc(a_str->strings, (a_str->nstrings + 1) * sizeof(*a_str->strings));
    
    if (NULL == str) return NULL;
    
    a_str->strings = str;
    str = &a_str->strings[a_str->nstrings];
    
    str->str = malloc(a_str->max_length);;
    
    if (NULL == str->str) return NULL;
    
    str->used = 0;
    a_str->nstrings++;
    
    return str;
}


/* free allocated memory */

void string_pool_destroy(string_alloc_t *a_str) {
    size_t i;
    
    for (i = 0; i < a_str->nstrings; i++) {
    	free(a_str->strings[i].str);
    }
    
    free(a_str->strings);
    free(a_str);
}


/* allocate space for a string */

char *string_alloc(string_alloc_t *a_str, size_t length) {
    string_t *str;
    char *ret;
    
    if (length <= 0) return NULL;
    
    // add to last string pool if we have space
    if (a_str->nstrings) {
    	str = &a_str->strings[a_str->nstrings - 1];
	
	if (str->used + length < a_str->max_length) {
	    ret = str->str + str->used;
	    str->used += length;
	    return ret;
	}
    }
    
    // increase the max length if needs be
    if (length > a_str->max_length) a_str->max_length = length;
	
    // need a new string pool 
    str = new_string_pool(a_str);
    
    if (NULL == str) return NULL;
    
    str->used = length;
    return str->str;
}


/* equivalent to strdup */

char *string_dup(string_alloc_t *a_str, char *instr) {
    return string_ndup(a_str, instr, strlen(instr));
}

char *string_ndup(string_alloc_t *a_str, char *instr, size_t len) {
    char *str = string_alloc(a_str, len + 1);
    
    if (NULL == str) return NULL;
    
    strncpy(str, instr, len);
    str[len] = 0;
    
    return str;
}
