/* The MIT License

   Copyright (c) 2003-2006, 2008, 2009 by Heng Li <lh3@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef LH3_KPROBALN_H_
#define LH3_KPROBALN_H_

#include <stdint.h>

typedef struct {
	float d, e;
	int bw;
} kpa_par_t;

#ifdef __cplusplus
extern "C" {
#endif

	int kpa_glocal(const uint8_t *_ref, int l_ref, const uint8_t *_query, int l_query, const uint8_t *iqual,
				   const kpa_par_t *c, int *state, uint8_t *q);

#ifdef __cplusplus
}
#endif

extern kpa_par_t kpa_par_def, kpa_par_alt;

#endif
