#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include "htslib/tbx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"

typedef struct
{
    int min_shift;
}
args_t;

static void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(EXIT_FAILURE);
}


#define IS_GFF (1<<0)
#define IS_BED (1<<1)
#define IS_SAM (1<<2)
#define IS_VCF (1<<3)
#define IS_BCF (1<<4)
#define IS_BAM (1<<5)
#define IS_TXT (IS_GFF|IS_BED|IS_SAM|IS_VCF)

int file_type(const char *fname)
{
    int l = strlen(fname);
    int strcasecmp(const char *s1, const char *s2);
    if (l>=7 && strcasecmp(fname+l-7, ".gff.gz") == 0) return IS_GFF;
    else if (l>=7 && strcasecmp(fname+l-7, ".bed.gz") == 0) return IS_BED;
    else if (l>=7 && strcasecmp(fname+l-7, ".sam.gz") == 0) return IS_SAM;
    else if (l>=7 && strcasecmp(fname+l-7, ".vcf.gz") == 0) return IS_VCF;
    else if (l>=4 && strcasecmp(fname+l-4, ".bcf") == 0) return IS_BCF;
    else if (l>=4 && strcasecmp(fname+l-4, ".bam") == 0) return IS_BAM;
    return 0;
}

#define PRINT_HEADER 1
#define HEADER_ONLY  2
static int query_regions(char **argv, int argc, int mode)
{
    char *fname = argv[0];
    int i, ftype = file_type(fname);

    if ( ftype & IS_TXT || !ftype )
    {
        htsFile *fp = hts_open(fname,"r");
        if ( !fp ) error("Could not read %s\n", fname);
        tbx_t *tbx = tbx_index_load(fname);
        if ( !tbx ) error("Could not load .tbi index of %s\n", fname);
        kstring_t str = {0,0,0};
        if ( mode )
        {
            while ( hts_getline(fp, KS_SEP_LINE, &str) >= 0 )
            {
                if ( !str.l || str.s[0]!=tbx->conf.meta_char ) break;
                puts(str.s);
            }
        }
        if ( mode!=HEADER_ONLY )
        {
            for (i=1; i<argc; i++)
            {
                hts_itr_t *itr = tbx_itr_querys(tbx, argv[i]);
                if ( !itr ) continue;
                while (tbx_itr_next(fp, tbx, itr, &str) >= 0) puts(str.s);
                tbx_itr_destroy(itr);
            }
        }
        free(str.s);
        if ( hts_close(fp) ) error("hts_close returned non-zero status: %s\n", fname);
        tbx_destroy(tbx);
    }
    else if ( ftype==IS_BCF )   // output uncompressed VCF
    {
        htsFile *fp = hts_open(fname,"r");
        if ( !fp ) error("Could not read %s\n", fname);
        htsFile *out = hts_open("-","w"); 
        if ( !out ) error("Could not open stdout\n", fname);
        hts_idx_t *idx = bcf_index_load(fname);
        if ( !idx ) error("Could not load .csi index of %s\n", fname);
        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        if ( !hdr ) error("Could not read the header: %s\n", fname);
        if ( mode )
        {
            bcf_hdr_write(out,hdr);
        }
        if ( mode!=HEADER_ONLY )
        {
            bcf1_t *rec = bcf_init();
            for (i=1; i<argc; i++)
            {
                hts_itr_t *itr = bcf_itr_querys(idx,hdr,argv[i]);
                if ( !itr ) continue;
                while ( bcf_itr_next(fp, itr, rec) >=0 ) bcf_write(out,hdr,rec);
                tbx_itr_destroy(itr);
            }
            bcf_destroy(rec);
        }
        if ( hts_close(fp) ) error("hts_close returned non-zero status: %s\n", fname);
        if ( hts_close(out) ) error("hts_close returned non-zero status for stdout\n");
        bcf_hdr_destroy(hdr);
        hts_idx_destroy(idx);
    }
    else if ( ftype==IS_BAM )   // todo: BAM
        error("Please use \"samtools view\" for querying BAM files.\n");
    return 0;
}
static int query_chroms(char *fname)
{
    const char **seq;
    int i, nseq, ftype = file_type(fname);
    if ( ftype & IS_TXT || !ftype )
    {
        tbx_t *tbx = tbx_index_load(fname);
        if ( !tbx ) error("Could not load .tbi index of %s\n", fname);
        seq = tbx_seqnames(tbx, &nseq);
        for (i=0; i<nseq; i++)
            printf("%s\n", seq[i]);
        free(seq);
        tbx_destroy(tbx);
    }
    else if ( ftype==IS_BCF )
    {
        htsFile *fp = hts_open(fname,"r");
        if ( !fp ) error("Could not read %s\n", fname);
        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        if ( !hdr ) error("Could not read the header: %s\n", fname);
        hts_close(fp);
        hts_idx_t *idx = bcf_index_load(fname);
        if ( !idx ) error("Could not load .csi index of %s\n", fname);
        seq = bcf_index_seqnames(idx, hdr, &nseq);
        for (i=0; i<nseq; i++) 
            printf("%s\n", seq[i]);
        free(seq);
        bcf_hdr_destroy(hdr);
        hts_idx_destroy(idx);
    }
    else if ( ftype==IS_BAM )   // todo: BAM
        error("BAM: todo\n");
    return 0;
}

int reheader_file(const char *fname, const char *header, int ftype, tbx_conf_t *conf)
{
    if ( ftype & IS_TXT || !ftype )
    {
        BGZF *fp = bgzf_open(fname,"r");
        if ( !fp || bgzf_read_block(fp) != 0 || !fp->block_length ) return -1;

        char *buffer = fp->uncompressed_block;
        int skip_until = 0;

        // Skip the header: find out the position of the data block
        if ( buffer[0]==conf->meta_char )
        {
            skip_until = 1;
            while (1)
            {
                if ( buffer[skip_until]=='\n' )
                {
                    skip_until++;
                    if ( skip_until>=fp->block_length )
                    {
                        if ( bgzf_read_block(fp) != 0 || !fp->block_length ) error("FIXME: No body in the file: %s\n", fname);
                        skip_until = 0;
                    }
                    // The header has finished
                    if ( buffer[skip_until]!=conf->meta_char ) break;
                }
                skip_until++;
                if ( skip_until>=fp->block_length )
                {
                    if (bgzf_read_block(fp) != 0 || !fp->block_length) error("FIXME: No body in the file: %s\n", fname);
                    skip_until = 0;
                }
            }
        }

        // Output the new header
        FILE *hdr  = fopen(header,"r");
        if ( !hdr ) error("%s: %s", header,strerror(errno));
        int page_size = getpagesize();
        char *buf = valloc(page_size);
        BGZF *bgzf_out = bgzf_dopen(fileno(stdout), "w");
        ssize_t nread;
        while ( (nread=fread(buf,1,page_size-1,hdr))>0 )
        {
            if ( nread<page_size-1 && buf[nread-1]!='\n' ) buf[nread++] = '\n';
            if (bgzf_write(bgzf_out, buf, nread) < 0) error("Error: %d\n",bgzf_out->errcode);
        }
        if ( fclose(hdr) ) error("close failed: %s\n", header);

        // Output all remainig data read with the header block
        if ( fp->block_length - skip_until > 0 )
        {
            if (bgzf_write(bgzf_out, buffer+skip_until, fp->block_length-skip_until) < 0) error("Error: %d\n",fp->errcode);
        }
        if (bgzf_flush(bgzf_out) < 0) error("Error: %d\n",bgzf_out->errcode);

        while (1)
        {
            nread = bgzf_raw_read(fp, buf, page_size);
            if ( nread<=0 ) break;

            int count = bgzf_raw_write(bgzf_out, buf, nread);
            if (count != nread) error("Write failed, wrote %d instead of %d bytes.\n", count,(int)nread);
        }
        if (bgzf_close(bgzf_out) < 0) error("Error: %d\n",bgzf_out->errcode);
        if (bgzf_close(fp) < 0) error("Error: %d\n",fp->errcode);
    }
    else
        error("todo: reheader BCF, BAM\n");  // BCF is difficult, records contain pointers to the header.
    return 0;
}

static int usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Version: %s\n", hts_version());
	fprintf(stderr, "Usage:   tabix [OPTIONS] [FILE] [REGION [...]]\n");
	fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -0, --zero-based        coordinates are zero-based\n");
    fprintf(stderr, "   -b, --begin INT         column number for region start [4]\n");
    fprintf(stderr, "   -c, --comment CHAR      skip comment lines starting with CHAR [null]\n");
    fprintf(stderr, "   -e, --end INT           column number for region end (if no end, set INT to -b) [5]\n");
    fprintf(stderr, "   -f, --force             overwrite existing index without asking\n");
    fprintf(stderr, "   -h, --print-header      print also the header lines\n");
    fprintf(stderr, "   -H, --only-header       print only the header lines\n");
    fprintf(stderr, "   -l, --list-chroms       list chromosome names\n");
    fprintf(stderr, "   -m, --min-shift INT     set the minimal interval size to 1<<INT; 0 for the old tabix index [0]\n");
    fprintf(stderr, "   -p, --preset STR        gff, bed, sam, vcf, bcf, bam\n");
    fprintf(stderr, "   -r, --reheader FILE     replace the header with the content of FILE\n");
    fprintf(stderr, "   -s, --sequence INT      column number for sequence names (suppressed by -p) [1]\n");
    fprintf(stderr, "   -S, --skip-lines INT    skip first INT lines [0]\n");
    fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	int c, min_shift = -1, is_force = 0, list_chroms = 0, mode = 0;
	tbx_conf_t conf = tbx_conf_gff, *conf_ptr = NULL;
    char *reheader = NULL;

    static struct option loptions[] = 
    {
        {"help",0,0,'h'},
        {"zero-based",0,0,'0'},
        {"print-header",0,0,'h'},
        {"only-header",0,0,'H'},
        {"begin",1,0,'b'},
        {"comment",1,0,'c'},
        {"end",1,0,'e'},
        {"force",0,0,'f'},
        {"preset",1,0,'p'},
        {"sequence",1,0,'s'},
        {"skip-lines",1,0,'S'},
        {"list-chroms",0,0,'l'},
        {"reheader",1,0,'r'},
        {0,0,0,0}
    };

	while ((c = getopt_long(argc, argv, "hH?0b:c:e:fm:p:s:S:lr:", loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'r': reheader = optarg; break;
            case 'h': mode = PRINT_HEADER; break;
            case 'H': mode = HEADER_ONLY; break;
            case 'l': list_chroms = 1; break;
            case '0': conf.preset |= TBX_UCSC; break;
            case 'b': conf.bc = atoi(optarg); break;
            case 'e': conf.ec = atoi(optarg); break;
            case 'c': conf.meta_char = *optarg; break;
            case 'f': is_force = 1; break;
            case 'm': min_shift = atoi(optarg); break;
            case 'p':
                      if (strcmp(optarg, "gff") == 0) conf_ptr = &tbx_conf_gff;
                      else if (strcmp(optarg, "bed") == 0) conf_ptr = &tbx_conf_bed;
                      else if (strcmp(optarg, "sam") == 0) conf_ptr = &tbx_conf_sam;
                      else if (strcmp(optarg, "vcf") == 0) conf_ptr = &tbx_conf_vcf;
                      else error("The preset string not recognised: '%s'\n", optarg);
                      break;
            case 's': conf.sc = atoi(optarg); break;
            case 'S': conf.line_skip = atoi(optarg); break;
            default: return usage();
        }
    }

	if ( optind==argc ) return usage();

    if ( list_chroms )
        return query_chroms(argv[optind]);

    if ( argc > optind+1 || mode==HEADER_ONLY )
        return query_regions(&argv[optind], argc-optind, mode);

    char *fname = argv[optind];
    int ftype = file_type(fname);
    if ( !conf_ptr )    // no preset given
    {
        if ( ftype==IS_GFF ) conf_ptr = &tbx_conf_gff; 
        else if ( ftype==IS_BED ) conf_ptr = &tbx_conf_bed;
        else if ( ftype==IS_SAM ) conf_ptr = &tbx_conf_sam;
        else if ( ftype==IS_VCF ) conf_ptr = &tbx_conf_vcf;
        else if ( ftype==IS_BCF )
        { 
            if ( min_shift <= 0 ) min_shift = 14;
        }
        else if ( ftype==IS_BAM )
        { 
            if ( min_shift <= 0 ) min_shift = 14;
        }
    }
    if ( reheader )
        return reheader_file(fname, reheader, ftype, conf_ptr);

    if ( conf_ptr )
        conf = *conf_ptr;

    char *suffix = min_shift <= 0 ? ".tbi" : (ftype==IS_BAM ? ".bai" : ".csi");
    char *idx_fname = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(idx_fname, fname), suffix);

    struct stat stat_tbi, stat_file;
    if ( !is_force && stat(idx_fname, &stat_tbi)==0 )
    {
        // Before complaining about existing index, check if the VCF file isn't
        // newer. This is a common source of errors, people tend not to notice
        // that tabix failed
        stat(fname, &stat_file);
        if ( stat_file.st_mtime <= stat_tbi.st_mtime )
            error("[tabix] the index file exists. Please use '-f' to overwrite.\n");
    }
    free(idx_fname);

    if ( min_shift > 0 ) // CSI index
    {
        if ( ftype==IS_BCF ) 
        {
            if ( bcf_index_build(fname, min_shift)!=0 ) error("bcf_index_build failed: %s\n", fname);
            return 0;
        }
        if ( ftype==IS_BAM )
        {
            if ( bam_index_build(fname, min_shift)!=0 ) error("bam_index_build failed: %s\n", fname);
            return 0;
        }
        if ( tbx_index_build(fname, min_shift, &conf)!=0 ) error("tbx_index_build failed: %s\n", fname);
        return 0;
    }
    else
    {
        if ( tbx_index_build(fname, min_shift, &conf) ) error("tbx_index_build failed: %s\n", fname);
        return 0;
    }
	return 0;
}
