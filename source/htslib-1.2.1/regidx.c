/* 
    Copyright (C) 2014 Genome Research Ltd.

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
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

#ifdef _WIN32
#include <string.h>
#endif
#include "common\cross_platform.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/khash_str2int.h"
#include "htslib/regidx.h"


#define LIDX_SHIFT 13   // number of insignificant index bits

// List of regions for one chromosome
typedef struct
{
    int *idx, nidx;
    int nregs, mregs;   // n:used, m:alloced
    reg_t *regs;
    void *payload;
}
reglist_t;

// Container of all sequences
struct _regidx_t
{
    int nseq, mseq;     // n:used, m:alloced
    reglist_t *seq;     // regions for each sequence
    void *seq2regs;     // hash for fast lookup from chr name to regions
    char **seq_names;
    regidx_free_f free;     // function to free any data allocated by regidx_parse_f
    regidx_parse_f parse;   // parse one input line
    void *usr;              // user data to pass to regidx_parse_f

    // temporary data for index initialization
    kstring_t str;
    int rid_prev, start_prev, end_prev;
    int payload_size;
    void *payload;
};

int regidx_seq_nregs(regidx_t *idx, const char *seq)
{
    int iseq;
    if ( khash_str2int_get(idx->seq2regs, seq, &iseq)!=0 ) return 0; // no such sequence
    return idx->seq[iseq].nregs;
}

int regidx_nregs(regidx_t *idx)
{
    int i, nregs = 0;
    for (i=0; i<idx->nseq; i++) nregs += idx->seq[i].nregs;
    return nregs;
}

char **regidx_seq_names(regidx_t *idx, int *n)
{
    *n = idx->nseq;
    return idx->seq_names;
}

int _regidx_build_index(regidx_t *idx)
{
    int iseq;
    for (iseq=0; iseq<idx->nseq; iseq++)
    {
        reglist_t *list = &idx->seq[iseq];
        int j,k, imax = 0;   // max index bin
        for (j=0; j<list->nregs; j++)
        {
            int ibeg = list->regs[j].start >> LIDX_SHIFT;
            int iend = list->regs[j].end >> LIDX_SHIFT;
            if ( imax < iend + 1 )
            {
                int old_imax = imax; 
                imax = iend + 1;
                kroundup32(imax);
                list->idx = (int*) realloc(list->idx, imax*sizeof(int));
                for (k=old_imax; k<imax; k++) list->idx[k] = -1;
            }
            if ( ibeg==iend )
            {
                if ( list->idx[ibeg]<0 ) list->idx[ibeg] = j;
            }
            else
            {
                for (k=ibeg; k<=iend; k++)
                    if ( list->idx[k]<0 ) list->idx[k] = j;
            }
            list->nidx = iend + 1;
        }
    }
    return 0;
}

int regidx_insert(regidx_t *idx, char *line)
{
    if ( !line )
        return _regidx_build_index(idx);

    char *chr_from, *chr_to;
    reg_t reg;
    int ret = idx->parse(line,&chr_from,&chr_to,&reg,idx->payload,idx->usr);
    if ( ret==-2 ) return -1;   // error
    if ( ret==-1 ) return 0;    // skip the line

    int rid;
    idx->str.l = 0;
    kputsn(chr_from, chr_to-chr_from+1, &idx->str);
    if ( khash_str2int_get(idx->seq2regs, idx->str.s, &rid)!=0 )
    {
        idx->nseq++;
        int m_prev = idx->mseq;
        hts_expand0(reglist_t,idx->nseq,idx->mseq,idx->seq);
        hts_expand0(char*,idx->nseq,m_prev,idx->seq_names);
        idx->seq_names[idx->nseq-1] = strdup(idx->str.s);
        rid = khash_str2int_inc(idx->seq2regs, idx->seq_names[idx->nseq-1]);
    }

    reglist_t *list = &idx->seq[rid];
    list->nregs++;
    int m_prev = list->mregs;
    hts_expand(reg_t,list->nregs,list->mregs,list->regs);
    list->regs[list->nregs-1] = reg;
    if ( idx->payload_size )
    {
        if ( m_prev < list->mregs ) list->payload = realloc(list->payload,idx->payload_size*list->mregs);
        memcpy((char*)list->payload + idx->payload_size*(list->nregs-1), idx->payload, idx->payload_size);
    }

    if ( idx->rid_prev==rid )
    {
        if ( idx->start_prev > reg.start || (idx->start_prev==reg.start && idx->end_prev>reg.end) ) 
        { 
            fprintf(stderr,"The regions are not sorted: %s:%d-%d is before %s:%d-%d\n", 
                idx->str.s,idx->start_prev+1,idx->end_prev+1,idx->str.s,reg.start+1,reg.end+1); 
            return -1;
        }
    }
    idx->rid_prev = rid;
    idx->start_prev = reg.start;
    idx->end_prev = reg.end;
    return 0;
}

regidx_t *regidx_init(const char *fname, regidx_parse_f parser, regidx_free_f free_f, size_t payload_size, void *usr_dat)
{
    if ( !parser )
    {
        if ( !fname ) parser = regidx_parse_tab;
        else
        {
            int len = strlen(fname);
            if ( len>=7 && !strcasecmp(".bed.gz",fname+len-7) )
                parser = regidx_parse_bed;
            else if ( len>=8 && !strcasecmp(".bed.bgz",fname+len-8) )
                parser = regidx_parse_bed;
            else if ( len>=4 && !strcasecmp(".bed",fname+len-4) )
                parser = regidx_parse_bed;
            else
                parser = regidx_parse_tab;
        }
    }

    regidx_t *idx = (regidx_t*) calloc(1,sizeof(regidx_t));
    idx->free  = free_f;
    idx->parse = parser;
    idx->usr   = usr_dat;
    idx->seq2regs = khash_str2int_init();
    idx->rid_prev   = -1;
    idx->start_prev = -1;
    idx->end_prev   = -1;
    idx->payload_size = payload_size;
    if ( payload_size ) idx->payload = malloc(payload_size);

    if ( !fname ) return idx;
    
    kstring_t str = {0,0,0};

    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) goto error;

    while ( hts_getline(fp, KS_SEP_LINE, &str) > 0 )
    {
        if ( regidx_insert(idx, str.s) ) goto error;
    }
    regidx_insert(idx, NULL);
    
    free(str.s);
    hts_close(fp);
    return idx;

error:
    free(str.s);
    if ( fp ) hts_close(fp);
    regidx_destroy(idx);
    return NULL;
}

void regidx_destroy(regidx_t *idx)
{
    int i, j;
    for (i=0; i<idx->nseq; i++)
    {
        reglist_t *list = &idx->seq[i];
        if ( idx->free )
        {
            for (j=0; j<list->nregs; j++)
                idx->free((char*)list->payload + idx->payload_size*j);
        }
        free(list->payload);
        free(list->regs);
        free(list->idx);
    }
    free(idx->seq_names);
    free(idx->seq);
    free(idx->str.s);
    free(idx->payload);
    khash_str2int_destroy_free(idx->seq2regs);
    free(idx);
}

int regidx_overlap(regidx_t *idx, const char *chr, uint32_t from, uint32_t to, regitr_t *itr)
{
    if ( itr ) itr->i = itr->n = 0;

    int iseq;
    if ( khash_str2int_get(idx->seq2regs, chr, &iseq)!=0 ) return 0; // no such sequence

    reglist_t *list = &idx->seq[iseq];
    if ( !list->nregs ) return 0;

    int i, ibeg = from>>LIDX_SHIFT; 
    int ireg = ibeg < list->nidx ? list->idx[ibeg] : list->idx[ list->nidx - 1 ];
    if ( ireg < 0 )
    {
        // linear search; if slow, replace with binary search
        if ( ibeg > list->nidx ) ibeg = list->nidx;
        for (i=ibeg - 1; i>=0; i--)
            if ( list->idx[i] >=0 ) break;
        ireg = i>=0 ? list->idx[i] : 0;
    }
    for (i=ireg; i<list->nregs; i++)
    {
        if ( list->regs[i].start > to ) return 0;   // no match
        if ( list->regs[i].end >= from && list->regs[i].start <= to ) break; // found
    }

    if ( i>=list->nregs ) return 0;   // no match

    if ( !itr ) return 1;

    itr->i = 0;
    itr->n = list->nregs - i;
    itr->reg = &idx->seq[iseq].regs[i];
    if ( idx->payload_size )
        itr->payload = (char*)idx->seq[iseq].payload + i*idx->payload_size;
    else
        itr->payload = NULL;

    return 1;
}

int regidx_parse_bed(const char *line, char **chr_beg, char **chr_end, reg_t *reg, void *payload, void *usr)
{
    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments
    
    char *se = ss;
    while ( *se && !isspace(*se) ) se++;
    if ( !*se ) { fprintf(stderr,"Could not parse bed line: %s\n", line); return -2; }

    *chr_beg = ss;
    *chr_end = se-1;

    ss = se+1;
    reg->start = strtol(ss, &se, 10);
    if ( ss==se ) { fprintf(stderr,"Could not parse bed line: %s\n", line); return -2; }

    ss = se+1;
    reg->end = strtol(ss, &se, 10) - 1;
    if ( ss==se ) { fprintf(stderr,"Could not parse bed line: %s\n", line); return -2; }
    
    return 0;
}

int regidx_parse_tab(const char *line, char **chr_beg, char **chr_end, reg_t *reg, void *payload, void *usr)
{
    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments
    
    char *se = ss;
    while ( *se && !isspace(*se) ) se++;
    if ( !*se ) { fprintf(stderr,"Could not parse bed line: %s\n", line); return -2; }

    *chr_beg = ss;
    *chr_end = se-1;

    ss = se+1;
    reg->start = strtol(ss, &se, 10) - 1;
    if ( ss==se ) { fprintf(stderr,"Could not parse bed line: %s\n", line); return -2; }

    if ( !se[0] || !se[1] )
        reg->end = reg->start;
    else
    {
        ss = se+1;
        reg->end = strtol(ss, &se, 10);
        if ( ss==se ) reg->end = reg->start;
        else reg->end--;
    }
    
    return 0;
}

