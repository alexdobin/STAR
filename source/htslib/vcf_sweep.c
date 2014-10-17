#include "htslib/vcf_sweep.h"
#include "htslib/bgzf.h"

#define SW_FWD 0
#define SW_BWD 1

struct _bcf_sweep_t
{
    htsFile *file;
    bcf_hdr_t *hdr;
    BGZF *fp;

    int direction;          // to tell if the direction has changed
    int block_size;         // the size of uncompressed data to hold in memory
    bcf1_t *rec;            // bcf buffer
    int nrec, mrec;         // number of used records; total size of the buffer
    int lrid, lpos, lnals, lals_len, mlals;   // to check uniqueness of a record
    char *lals;

    uint64_t *idx;          // uncompressed offsets of VCF/BCF records
    int iidx, nidx, midx;   // i: current offset; n: used; m: allocated
    int idx_done;           // the index is built during the first pass
};

BGZF *hts_get_bgzfp(htsFile *fp);
int hts_useek(htsFile *file, long uoffset, int where);
long hts_utell(htsFile *file);

static inline int sw_rec_equal(bcf_sweep_t *sw, bcf1_t *rec)
{
    if ( sw->lrid!=rec->rid ) return 0;
    if ( sw->lpos!=rec->pos ) return 0;
    if ( sw->lnals!=rec->n_allele ) return 0;

    char *t = rec->d.allele[sw->lnals-1];
    int len = t - rec->d.allele[0] + 1;
    while ( *t ) { t++; len++; }
    if ( sw->lals_len!=len ) return 0;
    if ( memcmp(sw->lals,rec->d.allele[0],len) ) return 0;
    return 1;
}

static void sw_rec_save(bcf_sweep_t *sw, bcf1_t *rec)
{
    sw->lrid  = rec->rid;
    sw->lpos  = rec->pos;
    sw->lnals = rec->n_allele;

    char *t = rec->d.allele[sw->lnals-1];
    int len = t - rec->d.allele[0] + 1;
    while ( *t ) { t++; len++; }
    sw->lals_len = len;
    hts_expand(char, len, sw->mlals, sw->lals);
    memcpy(sw->lals, rec->d.allele[0], len);
}

static void sw_fill_buffer(bcf_sweep_t *sw)
{
    if ( !sw->iidx ) return;
    sw->iidx--;

    int ret = hts_useek(sw->file, sw->idx[sw->iidx], 0);
    assert( ret==0 );

    sw->nrec = 0;
    bcf1_t *rec = &sw->rec[sw->nrec];
    while ( (ret=bcf_read1(sw->file, sw->hdr, rec))==0 )
    {
        bcf_unpack(rec, BCF_UN_STR);

        // if not in the last block, stop at the saved record
        if ( sw->iidx+1 < sw->nidx && sw_rec_equal(sw,rec) ) break;

        sw->nrec++;
        hts_expand0(bcf1_t, sw->nrec+1, sw->mrec, sw->rec);
        rec = &sw->rec[sw->nrec];
    }
    sw_rec_save(sw, &sw->rec[0]);
}

bcf_sweep_t *bcf_sweep_init(const char *fname)
{
    bcf_sweep_t *sw = (bcf_sweep_t*) calloc(1,sizeof(bcf_sweep_t));
    sw->file = hts_open(fname, "r");
    sw->fp   = hts_get_bgzfp(sw->file);
    bgzf_index_build_init(sw->fp);
    sw->hdr  = bcf_hdr_read(sw->file);
    sw->mrec = 1;
    sw->rec  = (bcf1_t*) calloc(sw->mrec,(sizeof(bcf1_t)));
    sw->block_size = 1024*1024*3;
    sw->direction = SW_FWD;
    return sw;
}

void bcf_empty1(bcf1_t *v);
void bcf_sweep_destroy(bcf_sweep_t *sw)
{
    int i;
    for (i=0; i<sw->mrec; i++) bcf_empty1(&sw->rec[i]);
    free(sw->idx);
    free(sw->rec);
    free(sw->lals);
    bcf_hdr_destroy(sw->hdr);
    hts_close(sw->file);
    free(sw);
}

static void sw_seek(bcf_sweep_t *sw, int direction)
{
    sw->direction = direction;
    if ( direction==SW_FWD )
        hts_useek(sw->file, sw->idx[0], 0);
    else
    {
        sw->iidx = sw->nidx;
        sw->nrec = 0;
    }
}

bcf1_t *bcf_sweep_fwd(bcf_sweep_t *sw)
{
    if ( sw->direction==SW_BWD ) sw_seek(sw, SW_FWD);

    long pos = hts_utell(sw->file);

    bcf1_t *rec = &sw->rec[0];
    int ret = bcf_read1(sw->file, sw->hdr, rec);

    if ( ret!=0 )   // last record, get ready for sweeping backwards
    {
        sw->idx_done = 1;
        sw->fp->idx_build_otf = 0;
        sw_seek(sw, SW_BWD);
        return NULL;
    }

    if ( !sw->idx_done )
    {
        if ( !sw->nidx || pos - sw->idx[sw->nidx-1] > sw->block_size )
        {
            sw->nidx++;
            hts_expand(uint64_t, sw->nidx, sw->midx, sw->idx);
            sw->idx[sw->nidx-1] = pos;
        }
    }
    return rec;
}

bcf1_t *bcf_sweep_bwd(bcf_sweep_t *sw)
{
    if ( sw->direction==SW_FWD ) sw_seek(sw, SW_BWD);
    if ( !sw->nrec ) sw_fill_buffer(sw);
    if ( !sw->nrec ) return NULL;
    return &sw->rec[ --sw->nrec ];
}

bcf_hdr_t *bcf_sweep_hdr(bcf_sweep_t *sw) { return sw->hdr; }

