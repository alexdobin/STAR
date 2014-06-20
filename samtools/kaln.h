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

#ifndef LH3_KALN_H_
#define LH3_KALN_H_

#include <stdint.h>

#define MINOR_INF -1073741823

typedef struct {
	int gap_open;
	int gap_ext;
	int gap_end_open;
	int gap_end_ext;

	int *matrix;
	int row;
	int band_width;
} ka_param_t;

typedef struct {
	int iio, iie, ido, ide;
	int eio, eie, edo, ede;
	int *matrix;
	int row;
	int band_width;
} ka_param2_t;

#ifdef __cplusplus
extern "C" {
#endif

	uint32_t *ka_global_core(uint8_t *seq1, int len1, uint8_t *seq2, int len2, const ka_param_t *ap,
							 int *_score, int *n_cigar);
	int ka_global_score(const uint8_t *_seq1, int len1, const uint8_t *_seq2, int len2, const ka_param2_t *ap);
#ifdef __cplusplus
}
#endif

extern ka_param_t ka_param_blast; /* = { 5, 2, 5, 2, aln_sm_blast, 5, 50 }; */
extern ka_param_t ka_param_qual; // only use this for global alignment!!!
extern ka_param2_t ka_param2_qual; // only use this for global alignment!!!

#endif
