/*  hts_defs.h -- Miscellaneous definitions.

    Copyright (C) 2013-2014 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

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

#ifndef HTSLIB_HTS_DEFS_H
#define HTSLIB_HTS_DEFS_H

#if __clang__major__ >= 2 || __GNUC__ >= 3
#define HTS_NORETURN __attribute__ ((__noreturn__))
#else
#define HTS_NORETURN
#endif

#if (defined __clang__ && __clang_major__ >= 3) || \
    (defined __GNUC__ && (__GNUC__ > 4 || (__GNUC__==4 && __GNUC_MINOR__ >= 5)))
#define HTS_RESULT_USED __attribute__ ((__warn_unused_result__))
#else
#define HTS_RESULT_USED
#endif

#if defined __clang__ || defined __GNUC__
#define HTS_UNUSED __attribute__ ((__unused__))
#else
#define HTS_UNUSED
#endif

#endif
