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

#ifndef _STRING_ALLOC_H_
#define _STRING_ALLOC_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

/* 
 * A pooled string allocator intended to cut down on the
 * memory overhead of many small string allocations.
 *
 * Andrew Whitwham, September 2010.
 */

typedef struct {
    char *str;
    size_t used;
} string_t;

typedef struct {
    size_t max_length;
    size_t nstrings;
    string_t *strings;
} string_alloc_t;

string_alloc_t *string_pool_create(size_t max_length);
void string_pool_destroy(string_alloc_t *a_str);
char *string_alloc(string_alloc_t *a_str, size_t length);
char *string_dup(string_alloc_t *a_str, char *instr);
char *string_ndup(string_alloc_t *a_str, char *instr, size_t len);

#endif

#ifdef __cplusplus
}
#endif

