/*  khash_str2int.h -- C-string to integer hash table.

    Copyright (C) 2013 Genome Research Ltd.

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

#ifndef HTSLIB_KHASH_STR2INT_H
#define HTSLIB_KHASH_STR2INT_H

#include <htslib/khash.h>

KHASH_MAP_INIT_STR(str2int, int)

/*
 *  Wrappers for khash dictionaries used by mpileup.
 */

static inline void *khash_str2int_init(void)
{
    return kh_init(str2int);
}

/*
 *  Destroy the hash structure, but not the keys
 */
static inline void khash_str2int_destroy(void *_hash)
{
    khash_t(str2int) *hash = (khash_t(str2int)*)_hash;
    if (hash) kh_destroy(str2int, hash); // Note that strings are not freed.
}

/*
 *  Destroys both the hash structure and the keys
 */
static inline void khash_str2int_destroy_free(void *_hash)
{
    khash_t(str2int) *hash = (khash_t(str2int)*)_hash;
    khint_t k;
    if (hash == 0) return;
    for (k = 0; k < kh_end(hash); ++k)
        if (kh_exist(hash, k)) free((char*)kh_key(hash, k));
    kh_destroy(str2int, hash);
}

/*
 *  Returns 1 if key exists or 0 if not
 */
static inline int khash_str2int_has_key(void *_hash, const char *str)
{
    khash_t(str2int) *hash = (khash_t(str2int)*)_hash;
    khint_t k = kh_get(str2int, hash, str);
    if ( k == kh_end(hash) ) return 0;
    return 1;
}

/*
 *  Returns 0 on success and -1 when the key is not present. On success,
 *  *value is set, unless NULL is passed.
 */
static inline int khash_str2int_get(void *_hash, const char *str, int *value)
{
    khash_t(str2int) *hash = (khash_t(str2int)*)_hash;
    khint_t k;
    if ( !hash ) return -1;
    k = kh_get(str2int, hash, str);
    if ( k == kh_end(hash) ) return -1;
    if ( !value ) return 0;
    *value = kh_val(hash, k);
    return 0;
}

/*
 *  Add a new string to the dictionary, auto-incrementing the value.
 *  On success returns the newly inserted integer id, on error -1
 *  is returned. Note that the key must continue to exist throughout
 *  the whole life of _hash.
 */
static inline int khash_str2int_inc(void *_hash, const char *str)
{
    khint_t k;
    int ret;
    khash_t(str2int) *hash = (khash_t(str2int)*)_hash;
    if ( !hash ) return -1;
    k = kh_put(str2int, hash, str, &ret);
    if (ret == 0) return kh_val(hash, k);
    kh_val(hash, k) = kh_size(hash) - 1;
    return kh_val(hash, k);
}

/*
 *  Set a new key,value pair. On success returns the bin index, on
 *  error -1 is returned. Note that the key must contnue to exist
 *  throughout the whole life of _hash.
 */
static inline int khash_str2int_set(void *_hash, const char *str, int value)
{
    khint_t k;
    int ret;
    khash_t(str2int) *hash = (khash_t(str2int)*)_hash;
    if ( !hash ) return -1;
    k = kh_put(str2int, hash, str, &ret);
    kh_val(hash,k) = value;
    return k;
}

/*
 *  Return the number of keys in the hash table.
 */
static inline int khash_str2int_size(void *_hash)
{
    khash_t(str2int) *hash = (khash_t(str2int)*)_hash;
    return kh_size(hash);
}

#endif
