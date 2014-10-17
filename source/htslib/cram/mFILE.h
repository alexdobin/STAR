/*
Copyright (c) 2005-2006, 2008-2009 Genome Research Ltd.
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

#ifndef _MFILE_H_
#define _MFILE_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    FILE *fp;
    char *data;
    size_t alloced;
    int eof;
    int mode; /* open mode in MF_?? define bit pattern */
    size_t size;
    size_t offset;
    size_t flush_pos;
} mFILE;

#define MF_READ    1
#define MF_WRITE   2
#define MF_APPEND  4
#define MF_BINARY  8
#define MF_TRUNC  16
#define MF_MODEX  32

mFILE *mfreopen(const char *path, const char *mode, FILE *fp);
mFILE *mfopen(const char *path, const char *mode);
int mfdetach(mFILE *mf);
int mfclose(mFILE *mf);
int mfdestroy(mFILE *mf);
int mfseek(mFILE *mf, long offset, int whence);
long mftell(mFILE *mf);
void mrewind(mFILE *mf);
void mftruncate(mFILE *mf, long offset);
int mfeof(mFILE *mf);
size_t mfread(void *ptr, size_t size, size_t nmemb, mFILE *mf);
size_t mfwrite(void *ptr, size_t size, size_t nmemb, mFILE *mf);
int mfgetc(mFILE *mf);
int mungetc(int c, mFILE *mf);
mFILE *mfcreate(char *data, int size);
mFILE *mfcreate_from(const char *path, const char *mode_str, FILE *fp);
void mfrecreate(mFILE *mf, char *data, int size);
void *mfsteal(mFILE *mf, size_t *size_out);
char *mfgets(char *s, int size, mFILE *mf);
int mfflush(mFILE *mf);
int mfprintf(mFILE *mf, char *fmt, ...);
mFILE *mstdin(void);
mFILE *mstdout(void);
mFILE *mstderr(void);
void mfascii(mFILE *mf);

#ifdef __cplusplus
}
#endif

#endif /* _MFILE_H_ */
