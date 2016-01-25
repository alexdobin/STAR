/*
Copyright (c) 1994-1997, 2001-2002 MEDICAL RESEARCH COUNCIL
All rights reserved

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   1 Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

   2 Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

   3 Neither the name of the MEDICAL RESEARCH COUNCIL, THE LABORATORY OF
MOLECULAR BIOLOGY nor the names of its contributors may be used to endorse or
promote products derived from this software without specific prior written
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
Copyright (c) 2003-2013 Genome Research Ltd.

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

#ifndef _misc_h
#define _misc_h

#include "cram/os.h"

#include <stdio.h>
#include <stdarg.h>  /* varargs needed for v*printf() prototypes */
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * This informs gcc that crash() doesn't return, so it doesn't need to
 * concern itself that code paths going via crash could mean some variables
 * being undefined and then issuing uninitialised variable warnings.
 * This particularly affected convert.
 */
#ifdef __GNUC__
#    define __NORETURN__ __attribute__ ((__noreturn__))
#else
#    define __NORETURN__
#endif

/*
 * Used for printf style argument checking. We can request a function such
 * as vTcl_SetResult does argument checking, avoiding bugs with using
 * %d and passing in a 64-bit record.
 */
#ifdef __GNUC__
#    define __PRINTF_FORMAT__(a,b) __attribute__ ((format (printf, a, b)))
#else
#    define __PRINTF_FORMAT__(a,b)
#endif

extern int is_directory(char * fn);
extern int is_file(char * fn);
extern int file_size(char * fn);

#define MIN(A,B) ( ( (A) < (B) ) ? (A) : (B) )
#define MAX(A,B) ( ( (A) > (B) ) ? (A) : (B) )

#ifdef __cplusplus
}
#endif

#endif /*_misc_h*/
