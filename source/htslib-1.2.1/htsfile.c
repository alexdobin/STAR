/*  htsfile.c -- file identifier and minimal viewer.

    Copyright (C) 2014-2015 Genome Research Ltd.

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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
#include "common\getopt\getopt.h"
#else
#include <getopt.h>
#endif
#include <unistd.h>
#include "common\cross_platform.h"
#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

enum { identify, view_headers, view_all } mode = identify;
int show_headers = 1;

static htsFile *dup_stdout(const char *mode)
{
    int fd = dup(STDOUT_FILENO);
    hFILE *hfp = (fd >= 0)? hdopen(fd, mode) : NULL;
    return hfp? hts_hopen(hfp, "-", mode) : NULL;
}

static int view_sam(hFILE *hfp, const char *filename)
{
    samFile *in = hts_hopen(hfp, filename, "r");
    if (in == NULL) return 0;
    samFile *out = dup_stdout("w");
    bam_hdr_t *hdr = sam_hdr_read(in);

    if (show_headers) sam_hdr_write(out, hdr);
    if (mode == view_all) {
        bam1_t *b = bam_init1();
        while (sam_read1(in, hdr, b) >= 0)
            sam_write1(out, hdr, b);
        bam_destroy1(b);
    }

    bam_hdr_destroy(hdr);
    hts_close(out);
    hts_close(in);
    return 1;
}

static int view_vcf(hFILE *hfp, const char *filename)
{
    vcfFile *in = hts_hopen(hfp, filename, "r");
    if (in == NULL) return 0;
    vcfFile *out = dup_stdout("w");
    bcf_hdr_t *hdr = bcf_hdr_read(in);

    if (show_headers) bcf_hdr_write(out, hdr);
    if (mode == view_all) {
        bcf1_t *rec = bcf_init();
        while (bcf_read(in, hdr, rec) >= 0)
            bcf_write(out, hdr, rec);
        bcf_destroy(rec);
    }

    bcf_hdr_destroy(hdr);
    hts_close(out);
    hts_close(in);
    return 1;
}

static void usage(FILE *fp, int status)
{
    fprintf(fp,
"Usage: htsfile [-chH] FILE...\n"
"Options:\n"
"  -c, --view         Write textual form of FILEs to standard output\n"
"  -h, --header-only  Display only headers in view mode, not records\n"
"  -H, --no-header    Suppress header display in view mode\n");
    exit(status);
}

int main(int argc, char **argv)
{
    static const struct option options[] = {
        { "header-only", no_argument, NULL, 'h' },
        { "no-header", no_argument, NULL, 'H' },
        { "view", no_argument, NULL, 'c' },
        { "help", no_argument, NULL, '?' },
        { "version", no_argument, NULL, 1 },
        { NULL, 0, NULL, 0 }
    };

    int status = EXIT_SUCCESS;
    int c, i;
    while ((c = getopt_long(argc, argv, "chH?", options, NULL)) >= 0)
        switch (c) {
        case 'c': mode = view_all; break;
        case 'h': mode = view_headers; show_headers = 1; break;
        case 'H': show_headers = 0; break;
        case 1:
            printf(
"htsfile (htslib) %s\n"
"Copyright (C) 2015 Genome Research Ltd.\n",
                   hts_version());
            exit(EXIT_SUCCESS);
            break;
        case '?': usage(stdout, EXIT_SUCCESS); break;
        default:  usage(stderr, EXIT_FAILURE); break;
        }

    if (optind == argc) usage(stderr, EXIT_FAILURE);

    for (i = optind; i < argc; i++) {
        htsFormat fmt;
        hFILE *fp = hopen(argv[i], "r");
        if (fp == NULL) {
            fprintf(stderr, "htsfile: can't open \"%s\": %s\n", argv[i], strerror(errno));
            status = EXIT_FAILURE;
            continue;
        }

        if (hts_detect_format(fp, &fmt) < 0) {
            fprintf(stderr, "htsfile: detecting \"%s\" format failed: %s\n", argv[i], strerror(errno));
            hclose_abruptly(fp);
            status = EXIT_FAILURE;
            continue;
        }

        if (mode == identify) {
            char *description = hts_format_description(&fmt);
            printf("%s:\t%s\n", argv[i], description);
            free(description);
        }
        else
            switch (fmt.category) {
            case sequence_data: if (view_sam(fp, argv[i])) fp = NULL; break;
            case variant_data:  if (view_vcf(fp, argv[i])) fp = NULL; break;
            default:
                fprintf(stderr, "htsfile: can't view %s: unknown format\n", argv[i]);
                status = EXIT_FAILURE;
                break;
            }

        if (fp && hclose(fp) < 0) {
            fprintf(stderr, "htsfile: closing %s failed\n", argv[i]);
            status = EXIT_FAILURE;
        }
    }

    return status;
}
