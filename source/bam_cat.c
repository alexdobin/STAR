/*

bam_cat -- efficiently concatenates bam files

bam_cat can be used to concatenate BAM files. Under special
circumstances, it can be used as an alternative to 'samtools merge' to
concatenate multiple sorted files into a single sorted file. For this
to work each file must be sorted, and the sorted files must be given
as command line arguments in order such that the final read in file i
is less than or equal to the first read in file i+1.

This code is derived from the bam_reheader function in samtools 0.1.8
and modified to perform concatenation by Chris Saunders on behalf of
Illumina.

2014-06-27: Alex Dobin modified the code slighlty:
* to compile with only htslib, no need for samtools package
* removed the samtools interface function (main)
* added header file "bam_cat.h"


########## License:

The MIT License

Original SAMtools work copyright (c) 2008-2009 Genome Research Ltd.
Modified SAMtools work copyright (c) 2010 Illumina, Inc.

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/

#include "bam_cat.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "htslib/htslib/bgzf.h"
#include "htslib/htslib/sam.h"
#include <cstring>

#define BUF_SIZE 0x10000

#define GZIPID1 31
#define GZIPID2 139

#define BGZF_EMPTY_BLOCK_SIZE 28


int bam_cat(int nfn, char * const *fn, const bam_hdr_t *h, const char* outbam)
{
    BGZF *fp;
    uint8_t *buf;
    uint8_t ebuf[BGZF_EMPTY_BLOCK_SIZE];
    const int es=BGZF_EMPTY_BLOCK_SIZE;
    int i;

    fp = strcmp(outbam, "-")? bgzf_open(outbam, "w") : bgzf_fdopen(fileno(stdout), "w");
    if (fp == 0) {
        fprintf(stderr, "[%s] ERROR: fail to open output file '%s'.\n", __func__, outbam);
        return 1;
    }
    if (h) bam_hdr_write(fp, h);

    buf = (uint8_t*) malloc(BUF_SIZE);
    for(i = 0; i < nfn; ++i){
        BGZF *in;
        bam_hdr_t *old;
        int len,j;

        in = strcmp(fn[i], "-")? bgzf_open(fn[i], "r") : bgzf_fdopen(fileno(stdin), "r");
        if (in == 0) {
            fprintf(stderr, "[%s] ERROR: fail to open file '%s'.\n", __func__, fn[i]);
            return -1;
        }
        if (in->is_write) return -1;

        old = bam_hdr_read(in);
        if (h == 0 && i == 0) bam_hdr_write(fp, old);

        if (in->block_offset < in->block_length) {
            bgzf_write(fp, (void*)((char*)in->uncompressed_block + in->block_offset), in->block_length - in->block_offset);
            bgzf_flush(fp);
        }

        j=0;
        while ((len = bgzf_raw_read(in, buf, BUF_SIZE)) > 0) {
            if(len<es){
                int diff=es-len;
                if(j==0) {
                    fprintf(stderr, "[%s] ERROR: truncated file?: '%s'.\n", __func__, fn[i]);
                    return -1;
                }
                bgzf_raw_write(fp, ebuf, len);
                memcpy(ebuf,ebuf+len,diff);
                memcpy(ebuf+diff,buf,len);
            } else {
                if(j!=0) bgzf_raw_write(fp, ebuf, es);
                len-= es;
                memcpy(ebuf,buf+len,es);
                bgzf_raw_write(fp, buf, len);
            }
            j=1;
        }

        /* check final gzip block */
        {
            const uint8_t gzip1=ebuf[0];
            const uint8_t gzip2=ebuf[1];
            const uint32_t isize=*((uint32_t*)(ebuf+es-4));
            if(((gzip1!=GZIPID1) || (gzip2!=GZIPID2)) || (isize!=0)) {
                fprintf(stderr, "[%s] WARNING: Unexpected block structure in file '%s'.", __func__, fn[i]);
                fprintf(stderr, " Possible output corruption.\n");
                bgzf_raw_write(fp, ebuf, es);
            }
        }
        bam_hdr_destroy(old);
        bgzf_close(in);
    }
    free(buf);
    bgzf_close(fp);
    return 0;
}



