#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <ctype.h>
#include <sys/stat.h>
#include "htslib/synced_bcf_reader.h"
#include "htslib/kseq.h"
#include "htslib/khash_str2int.h"

#define MAX_CSI_COOR 0x7fffffff     // maximum indexable coordinate of .csi

typedef struct
{
    uint32_t start, end;
}
region1_t;

typedef struct _region_t 
{
    region1_t *regs;
    int nregs, mregs, creg;
}
region_t;

static void _regions_add(bcf_sr_regions_t *reg, const char *chr, int start, int end);
static bcf_sr_regions_t *_regions_init_string(const char *str);
static int _regions_match_alleles(bcf_sr_regions_t *reg, int als_idx, bcf1_t *rec);

static int *init_filters(bcf_hdr_t *hdr, const char *filters, int *nfilters)
{
    kstring_t str = {0,0,0};
    const char *tmp = filters, *prev = filters;
    int nout = 0, *out = NULL;
    while ( 1 )
    {
        if ( *tmp==',' || !*tmp )
        {
            out = (int*) realloc(out, sizeof(int));
            if ( tmp-prev==1 && *prev=='.' )
                out[nout] = -1;
            else
            {
                str.l = 0;
                kputsn(prev, tmp-prev, &str);
                out[nout] = bcf_hdr_id2int(hdr, BCF_DT_ID, str.s);
            }
            nout++;
            if ( !*tmp ) break;
            prev = tmp+1;
        }
        tmp++;
    }
    if ( str.m ) free(str.s);
    *nfilters = nout;
    return out;
}

int bcf_sr_set_regions(bcf_srs_t *readers, const char *regions, int is_file)
{
    assert( !readers->regions );
    if ( readers->nreaders ) 
    {
        fprintf(stderr,"[%s:%d %s] Error: bcf_sr_set_regions() must be called before bcf_sr_add_reader()\n", __FILE__,__LINE__,__FUNCTION__);
        return -1;
    }
    readers->regions = bcf_sr_regions_init(regions,is_file,0,1,-2);
    if ( !readers->regions ) return -1;
    readers->explicit_regs = 1;
    readers->require_index = 1;
    return 0;
}
int bcf_sr_set_targets(bcf_srs_t *readers, const char *targets, int is_file, int alleles)
{
    assert( !readers->targets );
    readers->targets = bcf_sr_regions_init(targets,is_file,0,1,-2);
    if ( !readers->targets ) return -1;
    readers->targets_als = alleles;
    return 0;
}

int bcf_sr_add_reader(bcf_srs_t *files, const char *fname)
{
    files->has_line = (int*) realloc(files->has_line, sizeof(int)*(files->nreaders+1));
    files->has_line[files->nreaders] = 0;
    files->readers  = (bcf_sr_t*) realloc(files->readers, sizeof(bcf_sr_t)*(files->nreaders+1));
    bcf_sr_t *reader = &files->readers[files->nreaders++];
    memset(reader,0,sizeof(bcf_sr_t));

    reader->file = hts_open(fname, "r");
    if ( !reader->file ) return 0;

    reader->type = reader->file->is_bin? FT_BCF : FT_VCF;
    if (reader->file->is_compressed) reader->type |= FT_GZ;

    if ( files->require_index )
    {
        if ( reader->type==FT_VCF_GZ ) 
        {
            reader->tbx_idx = tbx_index_load(fname);
            if ( !reader->tbx_idx )
            {
                fprintf(stderr,"[add_reader] Could not load the index of %s\n", fname);
                return 0;
            }

            reader->header = bcf_hdr_read(reader->file);
        }
        else if ( reader->type==FT_BCF_GZ ) 
        {
            reader->header = bcf_hdr_read(reader->file);

            reader->bcf_idx = bcf_index_load(fname);
            if ( !reader->bcf_idx ) 
            {
                fprintf(stderr,"[add_reader] Could not load the index of %s\n", fname);
                return 0;   // not indexed..?
            }
        }
        else
        {
            fprintf(stderr,"Index required, expected .vcf.gz or .bcf file: %s\n", fname);
            return 0;
        }
    }
    else 
    {
        if ( reader->type & FT_BCF )
        {
            reader->header = bcf_hdr_read(reader->file);
        }
        else if ( reader->type & FT_VCF )
        {
            reader->header = bcf_hdr_read(reader->file);
        }
        else
        {
            fprintf(stderr,"File type not recognised: %s\n", fname);
            return 0;
        }
        files->streaming = 1;
    }
    if ( files->streaming && files->nreaders>1 )
    {
        fprintf(stderr,"[%s:%d %s] Error: %d readers, yet require_index not set\n", __FILE__,__LINE__,__FUNCTION__,files->nreaders);
        return 0;
    }
    if ( files->streaming && files->regions )
    {
        fprintf(stderr,"[%s:%d %s] Error: cannot tabix-jump in streaming mode\n", __FILE__,__LINE__,__FUNCTION__);
        return 0;
    }
    if ( !reader->header ) return 0;

    reader->fname = fname;
    if ( files->apply_filters )
        reader->filter_ids = init_filters(reader->header, files->apply_filters, &reader->nfilter_ids);

    // Update list of chromosomes
    if ( !files->explicit_regs && !files->streaming )
    {
        int n,i;
        const char **names = reader->tbx_idx ? tbx_seqnames(reader->tbx_idx, &n) : bcf_hdr_seqnames(reader->header, &n);
        for (i=0; i<n; i++)
        {
            if ( !files->regions )
                files->regions = _regions_init_string(names[i]);
            else
                _regions_add(files->regions, names[i], -1, -1);
        }
        free(names);
    }

    return 1;
}

bcf_srs_t *bcf_sr_init(void)
{
    bcf_srs_t *files = (bcf_srs_t*) calloc(1,sizeof(bcf_srs_t));
    return files;
}

static void bcf_sr_destroy1(bcf_sr_t *reader)
{
    if ( reader->tbx_idx ) tbx_destroy(reader->tbx_idx);
    if ( reader->bcf_idx ) hts_idx_destroy(reader->bcf_idx);
    bcf_hdr_destroy(reader->header);
    hts_close(reader->file);
    if ( reader->itr ) tbx_itr_destroy(reader->itr);
    int j;
    for (j=0; j<reader->mbuffer; j++)
        bcf_destroy1(reader->buffer[j]);
    free(reader->buffer);
    free(reader->samples);
    free(reader->filter_ids);
}
void bcf_sr_destroy(bcf_srs_t *files)
{
    int i;
    for (i=0; i<files->nreaders; i++)
        bcf_sr_destroy1(&files->readers[i]);
    free(files->has_line);
    free(files->readers);
    for (i=0; i<files->n_smpl; i++) free(files->samples[i]);
    free(files->samples);
    if (files->targets) bcf_sr_regions_destroy(files->targets);
    if (files->regions) bcf_sr_regions_destroy(files->regions);
    if ( files->tmps.m ) free(files->tmps.s);
    free(files);
}

void bcf_sr_remove_reader(bcf_srs_t *files, int i)
{
    assert( !files->samples );  // not ready for this yet
    bcf_sr_destroy1(&files->readers[i]);
    if ( i+1 < files->nreaders ) 
    {
        memmove(&files->readers[i], &files->readers[i+1], (files->nreaders-i-1)*sizeof(bcf_sr_t));
        memmove(&files->has_line[i], &files->has_line[i+1], (files->nreaders-i-1)*sizeof(int));
    }
    files->nreaders--;
}


/*
   Removes duplicate records from the buffer. The meaning of "duplicate" is
   controlled by the $collapse variable, which can cause that from multiple
   <indel|snp|any> lines only the first is considered and the rest is ignored.
   The removal is done by setting the redundant lines' positions to -1 and
   moving these lines at the end of the buffer.
 */
static void collapse_buffer(bcf_srs_t *files, bcf_sr_t *reader)
{
    int irec,jrec, has_snp=0, has_indel=0, has_any=0;
    for (irec=1; irec<=reader->nbuffer; irec++)
    {
        bcf1_t *line = reader->buffer[irec];
        if ( line->pos != reader->buffer[1]->pos ) break;
        if ( files->collapse&COLLAPSE_ANY )
        {
            if ( !has_any ) has_any = 1;
            else line->pos = -1;
        }
        int line_type = bcf_get_variant_types(line);
        if ( files->collapse&COLLAPSE_SNPS && line_type&(VCF_SNP|VCF_MNP) )
        {
            if ( !has_snp ) has_snp = 1;
            else line->pos = -1;
        }
        if ( files->collapse&COLLAPSE_INDELS && line_type&VCF_INDEL )
        {
            if ( !has_indel ) has_indel = 1;
            else line->pos = -1;
        }
    }
    bcf1_t *tmp;
    irec = jrec = 1;
    while ( irec<=reader->nbuffer && jrec<=reader->nbuffer )
    {
        if ( reader->buffer[irec]->pos != -1 ) { irec++; continue; }
        if ( jrec<=irec ) jrec = irec+1;
        while ( jrec<=reader->nbuffer && reader->buffer[jrec]->pos==-1 ) jrec++;
        if ( jrec<=reader->nbuffer )
        {
            tmp = reader->buffer[irec]; reader->buffer[irec] = reader->buffer[jrec]; reader->buffer[jrec] = tmp;
        }
    }
    reader->nbuffer = irec - 1;
}

void debug_buffer(FILE *fp, bcf_sr_t *reader)
{
    int j;
    for (j=0; j<=reader->nbuffer; j++)
    {
        bcf1_t *line = reader->buffer[j];
        fprintf(fp,"%s%s\t%s:%d\t%s ", reader->fname,j==0?"*":"",reader->header->id[BCF_DT_CTG][line->rid].key,line->pos+1,line->n_allele?line->d.allele[0]:"");
        int k;
        for (k=1; k<line->n_allele; k++) fprintf(fp," %s", line->d.allele[k]);
        fprintf(fp,"\n");
    }
}

void debug_buffers(FILE *fp, bcf_srs_t *files)
{
    int i;
    for (i=0; i<files->nreaders; i++)
    {
        fprintf(fp, "has_line: %d\t%s\n", bcf_sr_has_line(files,i),files->readers[i].fname);
        debug_buffer(fp, &files->readers[i]);
    }
    fprintf(fp,"\n");
}

static inline int has_filter(bcf_sr_t *reader, bcf1_t *line)
{
    int i, j;
    if ( !line->d.n_flt )
    {
        for (j=0; j<reader->nfilter_ids; j++)
            if ( reader->filter_ids[j]<0 ) return 1;
        return 0;
    }
    for (i=0; i<line->d.n_flt; i++)
    {
        for (j=0; j<reader->nfilter_ids; j++)
            if ( line->d.flt[i]==reader->filter_ids[j] ) return 1;
    }
    return 0;
}

static int _reader_seek(bcf_sr_t *reader, const char *seq, int start, int end)
{
    if ( end>=MAX_CSI_COOR )
    {
        fprintf(stderr,"The coordinate is out of csi index limit: %d\n", end+1);
        exit(1);
    }
    if ( reader->itr ) 
    {
        hts_itr_destroy(reader->itr); 
        reader->itr = NULL; 
    }
    reader->nbuffer = 0;
    if ( reader->tbx_idx )
    {
        int tid = tbx_name2id(reader->tbx_idx, seq);
        if ( tid==-1 ) return -1;    // the sequence not present in this file
        reader->itr = tbx_itr_queryi(reader->tbx_idx,tid,start,end+1);
    }
    else
    {
        int tid = bcf_hdr_name2id(reader->header, seq);
        if ( tid==-1 ) return -1;    // the sequence not present in this file
        reader->itr = bcf_itr_queryi(reader->bcf_idx,tid,start,end+1);
    }
    assert(reader->itr);
    return 0;
}

/*
 *  _readers_next_region() - jumps to next region if necessary
 *  Returns 0 on success or -1 when there are no more regions left
 */
static int _readers_next_region(bcf_srs_t *files)
{
    // Need to open new chromosome? Check number of lines in all readers' buffers
    int i, eos = 0;
    for (i=0; i<files->nreaders; i++)
        if ( !files->readers[i].itr && !files->readers[i].nbuffer ) eos++;

    if ( eos!=files->nreaders )
    {
        // Some of the readers still has buffered lines
        return 0;
    }

    // No lines in the buffer, need to open new region or quit
    if ( bcf_sr_regions_next(files->regions)<0 ) return -1;

    for (i=0; i<files->nreaders; i++)
        _reader_seek(&files->readers[i],files->regions->seq_names[files->regions->iseq],files->regions->start,files->regions->end);

    return 0;
}

/*
 *  _reader_fill_buffer() - buffers all records with the same coordinate
 */
static void _reader_fill_buffer(bcf_srs_t *files, bcf_sr_t *reader)
{
    // Return if the buffer is full: the coordinate of the last buffered record differs
    if ( reader->nbuffer && reader->buffer[reader->nbuffer]->pos != reader->buffer[1]->pos ) return;

    // No iterator (sequence not present in this file) and not streaming
    if ( !reader->itr && !files->streaming ) return;

    // Fill the buffer with records starting at the same position
    int i, ret = 0;
    while (1)
    {
        if ( reader->nbuffer+1 >= reader->mbuffer ) 
        {
            // Increase buffer size
            reader->mbuffer += 8;
            reader->buffer = (bcf1_t**) realloc(reader->buffer, sizeof(bcf1_t*)*reader->mbuffer);
            for (i=8; i>0; i--)     // initialize
            {
                reader->buffer[reader->mbuffer-i] = bcf_init1();
                reader->buffer[reader->mbuffer-i]->max_unpack = files->max_unpack;
                reader->buffer[reader->mbuffer-i]->pos = -1;    // for rare cases when VCF starts from 1
            }
        }
        if ( files->streaming )
        {
            if ( reader->type & FT_VCF )
            {
                if ( (ret=hts_getline(reader->file, KS_SEP_LINE, &files->tmps)) < 0 ) break;   // no more lines
                int ret = vcf_parse1(&files->tmps, reader->header, reader->buffer[reader->nbuffer+1]);
                if ( ret<0 ) break;
            }
            else if ( reader->type & FT_BCF )
            {
                if ( (ret=bcf_read1(reader->file, reader->header, reader->buffer[reader->nbuffer+1])) < 0 ) break; // no more lines
            }
            else
            {
                fprintf(stderr,"[%s:%d %s] fixme: not ready for this\n", __FILE__,__LINE__,__FUNCTION__);
                exit(1);
            }
        }
        else if ( reader->tbx_idx )
        {
            if ( (ret=tbx_itr_next(reader->file, reader->tbx_idx, reader->itr, &files->tmps)) < 0 ) break;  // no more lines
            vcf_parse1(&files->tmps, reader->header, reader->buffer[reader->nbuffer+1]);
        }
        else
        {
            if ( (ret=bcf_itr_next(reader->file, reader->itr, reader->buffer[reader->nbuffer+1])) < 0 ) break; // no more lines
            bcf_subset_format(reader->header,reader->buffer[reader->nbuffer+1]);
        }

        // apply filter
        if ( !reader->nfilter_ids )
            bcf_unpack(reader->buffer[reader->nbuffer+1], BCF_UN_STR);
        else
        {
            bcf_unpack(reader->buffer[reader->nbuffer+1], BCF_UN_STR|BCF_UN_FLT);
            if ( !has_filter(reader, reader->buffer[reader->nbuffer+1]) ) continue;
        }
        reader->nbuffer++;

        if ( reader->buffer[reader->nbuffer]->pos != reader->buffer[1]->pos ) break;    // the buffer is full
    }
    if ( ret<0 ) 
    { 
        // done for this region
        tbx_itr_destroy(reader->itr);
        reader->itr = NULL; 
    }
    if ( files->collapse && reader->nbuffer>=2 && reader->buffer[1]->pos==reader->buffer[2]->pos )
        collapse_buffer(files, reader);
}

/*
 *  _readers_shift_buffer() - removes the first line and all subsequent lines with the same position
 */
static void _reader_shift_buffer(bcf_sr_t *reader)
{
    int i;
    for (i=2; i<=reader->nbuffer; i++)
        if ( reader->buffer[i]->pos!=reader->buffer[1]->pos ) break;
    if ( i<=reader->nbuffer )
    {
        // A record with a different position follows, swap it. Because of the reader's logic,
        // only one such line can be present.
        bcf1_t *tmp = reader->buffer[1]; reader->buffer[1] = reader->buffer[i]; reader->buffer[i] = tmp;
        reader->nbuffer = 1;
    }
    else 
        reader->nbuffer = 0;    // no other line
}

/*
 *  _reader_match_alleles() - from multiple buffered lines selects the one which
 *  corresponds best to the template line. The logic is controlled by COLLAPSE_*
 *  Returns 0 on success or -1 when no good matching line is found.
 */
static int _reader_match_alleles(bcf_srs_t *files, bcf_sr_t *reader, bcf1_t *tmpl)
{
    int i, irec = -1;

    // if no template given, use the first available record
    if ( !tmpl )
        irec = 1;
    else
    {
        int tmpl_type = bcf_get_variant_types(tmpl);
        for (i=1; i<=reader->nbuffer; i++)
        {
            bcf1_t *line = reader->buffer[i];
            if ( line->pos != reader->buffer[1]->pos ) break;  // done with this reader

            // Easiest case: matching by position only
            if ( files->collapse&COLLAPSE_ANY ) { irec=i; break; }

            int line_type = bcf_get_variant_types(line);

            // No matter what the alleles are, as long as they are both SNPs
            if ( files->collapse&COLLAPSE_SNPS && tmpl_type&VCF_SNP && line_type&VCF_SNP ) { irec=i; break; }
            // ... or indels
            if ( files->collapse&COLLAPSE_INDELS && tmpl_type&VCF_INDEL && line_type&VCF_INDEL ) { irec=i; break; }

            // More thorough checking: REFs must match
            if ( tmpl->rlen != line->rlen ) continue;  // different length
            if ( strcmp(tmpl->d.allele[0], line->d.allele[0]) ) continue; // the strings do not match

            int ial,jal;
            if ( files->collapse==COLLAPSE_NONE )
            {
                // Exact match, all alleles must be identical
                if ( tmpl->n_allele!=line->n_allele ) continue;   // different number of alleles, skip

                int nmatch = 1; // REF has been already checked
                for (ial=1; ial<tmpl->n_allele; ial++)
                {
                    for (jal=1; jal<line->n_allele; jal++)
                        if ( !strcmp(tmpl->d.allele[ial], line->d.allele[jal]) ) { nmatch++; break; }
                }
                if ( nmatch==tmpl->n_allele ) { irec=i; break; }    // found: exact match
                continue;
            }
            
            // COLLAPSE_SOME: at least some ALTs must match
            for (ial=1; ial<tmpl->n_allele; ial++)
            {
                for (jal=1; jal<line->n_allele; jal++)
                    if ( !strcmp(tmpl->d.allele[ial], line->d.allele[jal]) ) { irec=i; break; }
                if ( irec>=1 ) break;
            }
            if ( irec>=1 ) break;
        }
        if ( irec==-1 ) return -1;  // no matching line was found
    }

    // Set the selected line (irec) as active: set it to buffer[0], move the remaining lines forward
    // and put the old bcf1_t record at the end.
    bcf1_t *tmp = reader->buffer[0];
    reader->buffer[0] = reader->buffer[irec];
    for (i=irec+1; i<=reader->nbuffer; i++) reader->buffer[i-1] = reader->buffer[i];
    reader->buffer[ reader->nbuffer ] = tmp;
    reader->nbuffer--;

    return 0;
}

int _reader_next_line(bcf_srs_t *files)
{
    int i, min_pos = INT_MAX;

    // Loop until next suitable line is found or all readers have finished
    while ( 1 )
    {
        // Get all readers ready for the next region.
        if ( files->regions && _readers_next_region(files)<0 ) break;

        // Fill buffers
        const char *chr = NULL;
        for (i=0; i<files->nreaders; i++)
        {
            _reader_fill_buffer(files, &files->readers[i]);

            // Update the minimum coordinate
            if ( !files->readers[i].nbuffer ) continue;
            if ( min_pos > files->readers[i].buffer[1]->pos ) 
            {
                min_pos = files->readers[i].buffer[1]->pos; 
                chr = bcf_seqname(files->readers[i].header, files->readers[i].buffer[1]);
            }
        }
        if ( min_pos==INT_MAX ) 
        {
            if ( !files->regions ) break;
            continue;
        }

        // Skip this position if not present in targets
        if ( files->targets )
        {
            if ( bcf_sr_regions_overlap(files->targets, chr, min_pos, min_pos)<0 ) 
            {
                // Remove all lines with this position from the buffer
                for (i=0; i<files->nreaders; i++)
                    if ( files->readers[i].nbuffer && files->readers[i].buffer[1]->pos==min_pos ) 
                        _reader_shift_buffer(&files->readers[i]);
                min_pos = INT_MAX;
                continue;
            }
        }
        
        break;  // done: min_pos is set 
    }

    // There can be records with duplicate positions. Set the active line intelligently so that
    // the alleles match.
    int nret = 0;   // number of readers sharing the position
    bcf1_t *first = NULL;   // record which will be used for allele matching
    for (i=0; i<files->nreaders; i++)
    {
        files->has_line[i] = 0;
        
        // Skip readers with no records at this position
        if ( !files->readers[i].nbuffer || files->readers[i].buffer[1]->pos!=min_pos ) continue;

        // Until now buffer[0] of all reader was empty and the lines started at buffer[1].
        // Now lines which are ready to be output will be moved to buffer[0].
        if ( _reader_match_alleles(files, &files->readers[i], first) < 0 ) continue;
        if ( !first ) first = files->readers[i].buffer[0];

        nret++;
        files->has_line[i] = 1;
    }
    return nret;
}

int bcf_sr_next_line(bcf_srs_t *files)
{
    if ( !files->targets_als ) 
        return _reader_next_line(files);

    while (1)
    {
        int i, ret = _reader_next_line(files);
        if ( !ret ) return ret;
        
        for (i=0; i<files->nreaders; i++)
            if ( files->has_line[i] ) break;

        if ( _regions_match_alleles(files->targets, files->targets_als-1, files->readers[i].buffer[0]) ) return ret;
        
        // Check if there are more duplicate lines in the buffers. If not, return this line as if it
        // matched the targets, even if there is a type mismatch
        for (i=0; i<files->nreaders; i++)
        {
            if ( !files->has_line[i] ) continue;
            if ( files->readers[i].nbuffer==0 || files->readers[i].buffer[1]->pos!=files->readers[i].buffer[0]->pos ) continue;
            break;
        }
        if ( i==files->nreaders ) return ret;   // no more lines left, output even if target alleles are not of the same type
    }
}

static void bcf_sr_seek_start(bcf_srs_t *readers)
{
    bcf_sr_regions_t *reg = readers->regions;
    int i;
    for (i=0; i<reg->nseqs; i++)
        reg->regs[i].creg = -1;
    reg->iseq = 0;
}


int bcf_sr_seek(bcf_srs_t *readers, const char *seq, int pos)
{
    if ( !seq && !pos ) 
    {
        // seek to start
        bcf_sr_seek_start(readers);
        return 0;
    }

    bcf_sr_regions_overlap(readers->regions, seq, pos, pos);
    int i, nret = 0;
    for (i=0; i<readers->nreaders; i++) 
    {
        nret += _reader_seek(&readers->readers[i],seq,pos,MAX_CSI_COOR-1);
    }
    return nret;
}

int bcf_sr_set_samples(bcf_srs_t *files, const char *fname, int is_file)
{
    int i, j, nsmpl, free_smpl = 0;
    char **smpl = NULL;
    
    void *exclude = (fname[0]=='^') ? khash_str2int_init() : NULL;
    if ( exclude || strcmp("-",fname) ) // "-" stands for all samples
    {
        smpl = hts_readlist(fname, is_file, &nsmpl);
        if ( !smpl ) 
        {
            fprintf(stderr,"Could not read the file: \"%s\"\n", fname);
            return 0;
        }
        if ( exclude )
        {
            for (i=0; i<nsmpl; i++)
                khash_str2int_inc(exclude, smpl[i]);
        }
        free_smpl = 1;
    }
    if ( !smpl )
    {
        smpl  = files->readers[0].header->samples;   // intersection of all samples
        nsmpl = bcf_hdr_nsamples(files->readers[0].header);
    }

    files->samples = NULL;
    files->n_smpl  = 0;
    for (i=0; i<nsmpl; i++)
    {
        if ( exclude && khash_str2int_has_key(exclude,smpl[i])  ) continue;

        int n_isec = 0;
        for (j=0; j<files->nreaders; j++)
        {
            if ( bcf_hdr_id2int(files->readers[j].header, BCF_DT_SAMPLE, smpl[i])<0 ) break;
            n_isec++;
        }
        if ( n_isec!=files->nreaders )
        {
            fprintf(stderr,"Warning: The sample \"%s\" was not found in %s, skipping\n", smpl[i], files->readers[n_isec].fname);
            continue;
        }

        files->samples = (char**) realloc(files->samples, (files->n_smpl+1)*sizeof(const char*));
        files->samples[files->n_smpl++] = strdup(smpl[i]);
    }

    if ( exclude ) khash_str2int_destroy(exclude);
    if ( free_smpl ) 
    {
        for (i=0; i<nsmpl; i++) free(smpl[i]);
        free(smpl);
    }

    if ( !files->n_smpl ) 
    {
        if ( files->nreaders>1 ) 
            fprintf(stderr,"No samples in common.\n");
        return 0;
    }
    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        reader->samples  = (int*) malloc(sizeof(int)*files->n_smpl);
        reader->n_smpl   = files->n_smpl;
        for (j=0; j<files->n_smpl; j++)
            reader->samples[j] = bcf_hdr_id2int(reader->header, BCF_DT_SAMPLE, files->samples[j]);
    }
    return 1;
}

// Add a new region into a list sorted by start,end. On input the coordinates
// are 1-based, stored 0-based, inclusive.
static void _regions_add(bcf_sr_regions_t *reg, const char *chr, int start, int end)
{
    if ( start==-1 && end==-1 )
    {
        start = 0; end = MAX_CSI_COOR-1;
    }
    else
    {
        start--; end--; // store 0-based coordinates
    }

    if ( !reg->seq_hash )
         reg->seq_hash = khash_str2int_init();

    int iseq;
    if ( khash_str2int_get(reg->seq_hash, chr, &iseq)<0 )
    {
        // the chromosome block does not exist
        iseq = reg->nseqs++;
        reg->seq_names = (char**) realloc(reg->seq_names,sizeof(char*)*reg->nseqs);
        reg->regs = (region_t*) realloc(reg->regs,sizeof(region_t)*reg->nseqs);
        memset(&reg->regs[reg->nseqs-1],0,sizeof(region_t));
        reg->seq_names[iseq] = strdup(chr);
        reg->regs[iseq].creg = -1;
        khash_str2int_set(reg->seq_hash,reg->seq_names[iseq],iseq);
    }

    region_t *creg = &reg->regs[iseq];

    // the regions may not be sorted on input: binary search
    int i, min = 0, max = creg->nregs - 1;
    while ( min<=max )
    {
        i = (max+min)/2;
        if ( start < creg->regs[i].start ) max = i - 1;
        else if ( start > creg->regs[i].start ) min = i + 1;
        else break;
    }
    if ( min>max || creg->regs[i].start!=start || creg->regs[i].end!=end )
    {
        // no such region, insert a new one just after max
        hts_expand(region1_t,creg->nregs+1,creg->mregs,creg->regs);
        if ( ++max < creg->nregs )
            memmove(&creg->regs[max+1],&creg->regs[max],(creg->nregs - max)*sizeof(region1_t));
        creg->regs[max].start = start;
        creg->regs[max].end   = end;
        creg->nregs++;
    }
}

// File name or a list of genomic locations. If file name, NULL is returned.
static bcf_sr_regions_t *_regions_init_string(const char *str)
{
    bcf_sr_regions_t *reg = (bcf_sr_regions_t *) calloc(1, sizeof(bcf_sr_regions_t));
    reg->start = reg->end = -1;
    reg->prev_start = reg->prev_seq = -1;

    kstring_t tmp = {0,0,0};
    const char *sp = str, *ep = str;
    int from, to;
    while ( 1 )
    {
        while ( *ep && *ep!=',' && *ep!=':' ) ep++;
        tmp.l = 0;
        kputsn(sp,ep-sp,&tmp);
        if ( *ep==':' )
        {
            sp = ep+1;
            from = strtol(sp,(char**)&ep,10);
            if ( sp==ep )
            {
                fprintf(stderr,"[%s:%d %s] Could not parse the region(s): %s\n", __FILE__,__LINE__,__FUNCTION__,str);
                free(reg); free(tmp.s); return NULL;
            }
            if ( !*ep || *ep==',' )
            {
                _regions_add(reg, tmp.s, from, from);
                sp = ep;
                continue;
            }
            if ( *ep!='-' )
            {
                fprintf(stderr,"[%s:%d %s] Could not parse the region(s): %s\n", __FILE__,__LINE__,__FUNCTION__,str);
                free(reg); free(tmp.s); return NULL;
            }
            ep++;
            sp = ep;
            to = strtol(sp,(char**)&ep,10);
            if ( *ep && *ep!=',' )
            {
                fprintf(stderr,"[%s:%d %s] Could not parse the region(s): %s\n", __FILE__,__LINE__,__FUNCTION__,str);
                free(reg); free(tmp.s); return NULL;
            }
            if ( sp==ep ) to = MAX_CSI_COOR-1;
            _regions_add(reg, tmp.s, from, to);
            if ( !*ep ) break;
            sp = ep;
        }
        else
        {
            if ( tmp.l ) _regions_add(reg, tmp.s, -1, -1);
            if ( !*ep ) break;
            sp = ++ep;
        }
    }
    free(tmp.s);
    return reg;
}

// ichr,ifrom,ito are 0-based;
// returns -1 on error, 0 if the line is a comment line, 1 on success
static int _regions_parse_line(char *line, int ichr,int ifrom,int ito, char **chr,char **chr_end,int *from,int *to)
{
    *chr_end = NULL;

    if ( line[0]=='#' ) return 0;

    int k,l;    // index of the start and end column of the tab-delimited file
    if ( ifrom <= ito ) 
        k = ifrom, l = ito;
    else 
        l = ifrom, k = ito;

    int i;
    char *se = line, *ss = NULL; // start and end 
    char *tmp;
    for (i=0; i<=k && *se; i++)
    {
        ss = i==0 ? se++ : ++se;
        while (*se && *se!='\t') se++;
    }
    if ( i<=k ) return -1;
    if ( k==l )
    {
        *from = *to = strtol(ss, &tmp, 10);
        if ( tmp==ss ) return -1;
    }
    else
    {
        if ( k==ifrom ) 
            *from = strtol(ss, &tmp, 10);
        else
            *to = strtol(ss, &tmp, 10);
        if ( ss==tmp ) return -1;

        for (i=k; i<l && *se; i++)
        {
            ss = ++se;
            while (*se && *se!='\t') se++;
        }
        if ( i<l ) return -1;
        if ( k==ifrom ) 
            *to = strtol(ss, &tmp, 10);
        else
            *from = strtol(ss, &tmp, 10);
        if ( ss==tmp ) return -1;
    }

    ss = se = line;
    for (i=0; i<=ichr && *se; i++)
    {
        if ( i>0 ) ss = ++se;
        while (*se && *se!='\t') se++;
    }
    if ( i<=ichr ) return -1;
    *chr_end = se;
    *chr = ss;
    return 1;
}

bcf_sr_regions_t *bcf_sr_regions_init(const char *regions, int is_file, int ichr, int ifrom, int ito)
{
    bcf_sr_regions_t *reg;
    if ( !is_file ) return _regions_init_string(regions);

    reg = (bcf_sr_regions_t *) calloc(1, sizeof(bcf_sr_regions_t));
    reg->start = reg->end = -1;
    reg->prev_start = reg->prev_seq = -1;

    reg->file = hts_open(regions, "rb");
    if ( !reg->file )
    {
        fprintf(stderr,"[%s:%d %s] Could not open file: %s\n", __FILE__,__LINE__,__FUNCTION__,regions);
        free(reg);
        return NULL;
    }

    reg->tbx = tbx_index_load(regions);
    if ( !reg->tbx ) 
    {
        int len = strlen(regions);
        int is_bed  = strcasecmp(".bed",regions+len-4) ? 0 : 1;
        if ( !is_bed && !strcasecmp(".bed.gz",regions+len-7) ) is_bed = 1;
        int ft_type = hts_file_type(regions);
        if ( ft_type & FT_VCF ) ito = 1;

        // read the whole file, tabix index is not present
        while ( hts_getline(reg->file, KS_SEP_LINE, &reg->line) > 0 )
        {
            char *chr, *chr_end;
            int from, to, ret;
            ret = _regions_parse_line(reg->line.s, ichr,ifrom,abs(ito), &chr,&chr_end,&from,&to);
            if ( ret < 0 ) 
            {
                if ( ito<0 )
                    ret = _regions_parse_line(reg->line.s, ichr,ifrom,ifrom, &chr,&chr_end,&from,&to);
                if ( ret<0 )
                {
                    fprintf(stderr,"[%s:%d] Could not parse the file %s, using the columns %d,%d[,%d]\n", __FILE__,__LINE__,regions,ichr+1,ifrom+1,ito+1);
                    hts_close(reg->file); reg->file = NULL; free(reg); 
                    return NULL;
                }
            }
            if ( !ret ) continue;
            if ( is_bed ) from++;
            *chr_end = 0;
            _regions_add(reg, chr, from, to);
            *chr_end = '\t';
        }
        hts_close(reg->file); reg->file = NULL;
        if ( !reg->nseqs ) { free(reg); return NULL; }
        return reg;
    }

    reg->seq_names = (char**) tbx_seqnames(reg->tbx, &reg->nseqs);
    if ( !reg->seq_hash )
        reg->seq_hash = khash_str2int_init();
    int i;
    for (i=0; i<reg->nseqs; i++)
    {
        khash_str2int_set(reg->seq_hash,reg->seq_names[i],i);
    }
    reg->fname  = strdup(regions);
    reg->is_bin = 1;
    return reg;
}

void bcf_sr_regions_destroy(bcf_sr_regions_t *reg)
{
    int i;
    free(reg->fname);
    if ( reg->itr ) tbx_itr_destroy(reg->itr);
    if ( reg->tbx ) tbx_destroy(reg->tbx);
    if ( reg->file ) hts_close(reg->file);
    if ( reg->als ) free(reg->als);
    if ( reg->als_str.s ) free(reg->als_str.s);
    free(reg->line.s);
    if ( reg->regs ) 
    {
         // free only in-memory names, tbx names are const
        for (i=0; i<reg->nseqs; i++) 
        {   
            free(reg->seq_names[i]);
            free(reg->regs[i].regs);
        }
    }
    free(reg->regs);
    free(reg->seq_names);
    khash_str2int_destroy(reg->seq_hash);
    free(reg);
}

int bcf_sr_regions_seek(bcf_sr_regions_t *reg, const char *seq)
{
    reg->iseq = reg->start = reg->end = -1;
    if ( khash_str2int_get(reg->seq_hash, seq, &reg->iseq) < 0 ) return -1;  // sequence seq not in regions

    // using in-memory regions
    if ( reg->regs ) return 0;

    // reading regions from tabix
    if ( reg->itr ) tbx_itr_destroy(reg->itr);
    reg->itr = tbx_itr_querys(reg->tbx, seq);
    if ( reg->itr ) return 0;

    return -1;
}

int bcf_sr_regions_next(bcf_sr_regions_t *reg)
{
    if ( reg->iseq<0 ) return -1;
    reg->start = reg->end = -1;
    reg->nals = 0;

    // using in-memory regions
    if ( reg->regs )
    {
        while ( reg->iseq < reg->nseqs )
        {
            reg->regs[reg->iseq].creg++;
            if ( reg->regs[reg->iseq].creg < reg->regs[reg->iseq].nregs ) break;
            reg->iseq++;
        }
        if ( reg->iseq >= reg->nseqs ) { reg->iseq = -1; return -1; } // no more regions left
        region1_t *creg = &reg->regs[reg->iseq].regs[reg->regs[reg->iseq].creg];
        reg->start = creg->start;
        reg->end   = creg->end;
        return 0;
    }

    // reading from tabix
    char *chr, *chr_end;
    int ichr = 0, ifrom = 1, ito = 2, is_bed = 0, from, to;
    if ( reg->tbx )
    {
        ichr   = reg->tbx->conf.sc-1;
        ifrom  = reg->tbx->conf.bc-1;
        ito    = reg->tbx->conf.ec-1;
        if ( ito<0 ) ito = ifrom;
        is_bed = reg->tbx->conf.preset==TBX_UCSC ? 1 : 0;
    }

    int ret = 0;
    while ( !ret )
    {
        if ( reg->itr )
        {
            // tabix index present, reading a chromosome block
            ret = tbx_itr_next(reg->file, reg->tbx, reg->itr, &reg->line);
            if ( ret<0 ) { reg->iseq = -1; return -1; }
        }
        else
        {
            if ( reg->is_bin )
            {
                // Waited for seek which never came. Reopen in text mode and stream
                // through the regions, otherwise hts_getline would fail
                hts_close(reg->file);
                reg->file = hts_open(reg->fname, "r");
                if ( !reg->file )
                {
                    fprintf(stderr,"[%s:%d %s] Could not open file: %s\n", __FILE__,__LINE__,__FUNCTION__,reg->fname);
                    reg->file = NULL;
                    bcf_sr_regions_destroy(reg);
                    return -1;
                }
                reg->is_bin = 0;
            }

            // tabix index absent, reading the whole file
            ret = hts_getline(reg->file, KS_SEP_LINE, &reg->line);
            if ( ret<0 ) { reg->iseq = -1; return -1; }
        }
        ret = _regions_parse_line(reg->line.s, ichr,ifrom,ito, &chr,&chr_end,&from,&to);
        if ( ret<0 ) 
        {
            fprintf(stderr,"[%s:%d] Could not parse the file %s, using the columns %d,%d,%d\n", __FILE__,__LINE__,reg->fname,ichr+1,ifrom+1,ito+1);
            return -1;
        }
    }
    if ( is_bed ) from++;

    *chr_end = 0;
    if ( khash_str2int_get(reg->seq_hash, chr, &reg->iseq)<0 )
    {
        fprintf(stderr,"Broken tabix index? The sequence \"%s\" not in dictionary [%s]\n", chr,reg->line.s);
        exit(1);
    }
    *chr_end = '\t';

    reg->start = from - 1;
    reg->end   = to - 1;
    return 0;
}

static int _regions_match_alleles(bcf_sr_regions_t *reg, int als_idx, bcf1_t *rec)
{
    int i = 0, max_len = 0;
    if ( !reg->nals )
    {
        char *ss = reg->line.s;
        while ( i<als_idx && *ss )
        {
            if ( *ss=='\t' ) i++;
            ss++;
        }
        char *se = ss;
        reg->nals = 1;
        while ( *se && *se!='\t' ) 
        { 
            if ( *se==',' ) reg->nals++;
            se++; 
        }
        ks_resize(&reg->als_str, se-ss+1+reg->nals);
        reg->als_str.l = 0;
        hts_expand(char*,reg->nals,reg->mals,reg->als);
        reg->nals = 0;

        se = ss;
        while ( *(++se) )
        {
            if ( *se=='\t' ) break;
            if ( *se!=',' ) continue;
            reg->als[reg->nals] = &reg->als_str.s[reg->als_str.l];
            kputsn(ss,se-ss,&reg->als_str);
            if ( &reg->als_str.s[reg->als_str.l] - reg->als[reg->nals] > max_len ) max_len = &reg->als_str.s[reg->als_str.l] - reg->als[reg->nals];
            reg->als_str.l++;
            reg->nals++;
            ss = ++se;
        }
        reg->als[reg->nals] = &reg->als_str.s[reg->als_str.l];
        kputsn(ss,se-ss,&reg->als_str);
        if ( &reg->als_str.s[reg->als_str.l] - reg->als[reg->nals] > max_len ) max_len = &reg->als_str.s[reg->als_str.l] - reg->als[reg->nals];
        reg->nals++;
        reg->als_type = max_len > 1 ? VCF_INDEL : VCF_SNP;  // this is a simplified check, see vcf.c:bcf_set_variant_types
    }
    int type = bcf_get_variant_types(rec);
    if ( reg->als_type & VCF_INDEL )
        return type & VCF_INDEL ? 1 : 0;
    return !(type & VCF_INDEL) ? 1 : 0;
}

int bcf_sr_regions_overlap(bcf_sr_regions_t *reg, const char *seq, int start, int end)
{
    int iseq;
    if ( khash_str2int_get(reg->seq_hash, seq, &iseq)<0 ) return -1;    // no such sequence

    if ( reg->prev_seq==-1 || iseq!=reg->prev_seq || reg->prev_start > start ) // new chromosome or after a seek
    {
        // flush regions left on previous chromosome
        if ( reg->missed_reg_handler && reg->prev_seq!=-1 && reg->iseq!=-1 )
            bcf_sr_regions_flush(reg);

        bcf_sr_regions_seek(reg, seq);
        reg->start = reg->end = -1;
    }
    if ( reg->prev_seq==iseq && reg->iseq!=iseq ) return -2;    // no more regions on this chromosome
    reg->prev_seq = reg->iseq;
    reg->prev_start = start;

    while ( iseq==reg->iseq && reg->end < start )
    {
        if ( bcf_sr_regions_next(reg) < 0 ) return -2;  // no more regions left
        if ( reg->iseq != iseq ) return -1; // does not overlap any regions
        if ( reg->missed_reg_handler && reg->end < start ) reg->missed_reg_handler(reg, reg->missed_reg_data);
    }
    if ( reg->start <= end ) return 0;    // region overlap
    return -1;  // no overlap
}

void bcf_sr_regions_flush(bcf_sr_regions_t *reg)
{
    if ( !reg->missed_reg_handler || reg->prev_seq==-1 ) return;
    while ( !bcf_sr_regions_next(reg) ) reg->missed_reg_handler(reg, reg->missed_reg_data);
    return;
}

