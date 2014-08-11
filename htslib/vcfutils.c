#include "htslib/vcfutils.h"

int bcf_calc_ac(const bcf_hdr_t *header, bcf1_t *line, int *ac, int which)
{
	int i;
	for (i=0; i<line->n_allele; i++) ac[i]=0;

	// Use INFO/AC,AN field only when asked
	if ( which&BCF_UN_INFO )
	{
		bcf_unpack(line, BCF_UN_INFO);
		int an_id = bcf_hdr_id2int(header, BCF_DT_ID, "AN");
		int ac_id = bcf_hdr_id2int(header, BCF_DT_ID, "AC");
        int i, an=-1, ac_len=0, ac_type=0;
        uint8_t *ac_ptr=NULL;
		if ( an_id>=0 && ac_id>=0 )
		{
			for (i=0; i<line->n_info; i++)
			{
				bcf_info_t *z = &line->d.info[i];
				if ( z->key == an_id ) an = z->v1.i;
				else if ( z->key == ac_id ) { ac_ptr = z->vptr; ac_len = z->len; ac_type = z->type; }
			}
        }
        if ( an>=0 && ac_ptr )
        {
			int nac = 0;
            #define BRANCH_INT(type_t) {        \
                type_t *p = (type_t *) ac_ptr;  \
                for (i=0; i<ac_len; i++)        \
                {                               \
                    ac[i+1] = p[i];             \
                    nac += p[i];                \
                }                               \
            }
            switch (ac_type) {
                case BCF_BT_INT8:  BRANCH_INT(int8_t); break;
                case BCF_BT_INT16: BRANCH_INT(int16_t); break;
                case BCF_BT_INT32: BRANCH_INT(int32_t); break;
                default: fprintf(stderr, "[E::%s] todo: %d at %s:%d\n", __func__, ac_type, header->id[BCF_DT_CTG][line->rid].key, line->pos+1); exit(1); break;
            }
            #undef BRANCH_INT
            assert( an>=nac );  // sanity check for missing values
			ac[0] = an - nac;
			return 1;
        }
	}

	// Split genotype fields only when asked
	if ( which&BCF_UN_FMT )
	{
		int i, gt_id = bcf_hdr_id2int(header,BCF_DT_ID,"GT");
		if ( gt_id<0 ) return 0;
		bcf_unpack(line, BCF_UN_FMT);
		bcf_fmt_t *fmt_gt = NULL;
		for (i=0; i<(int)line->n_fmt; i++) 
			if ( line->d.fmt[i].id==gt_id ) { fmt_gt = &line->d.fmt[i]; break; }
		if ( !fmt_gt ) return 0;
        #define BRANCH_INT(type_t,missing,vector_end) { \
		    for (i=0; i<line->n_sample; i++) \
		    { \
                type_t *p = (type_t*) (fmt_gt->p + i*fmt_gt->size); \
		    	int ial; \
		    	for (ial=0; ial<fmt_gt->n; ial++) \
		    	{ \
                    if ( p[ial]==vector_end ) break; /* smaller ploidy */ \
                    if ( !(p[ial]>>1) || p[ial]==missing ) continue; /* missing allele */ \
		    		ac[(p[ial]>>1)-1]++; \
		    	} \
		    } \
        }
        switch (fmt_gt->type) {
            case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
            case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
            case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
            default: fprintf(stderr, "[E::%s] todo: %d at %s:%d\n", __func__, fmt_gt->type, header->id[BCF_DT_CTG][line->rid].key, line->pos+1); exit(1); break;
        }
        #undef BRANCH_INT
		return 1;
	}
	return 0;
}

int bcf_gt_type(bcf_fmt_t *fmt_ptr, int isample, int *_ial, int *_jal)
{
    int i, nals = 0, has_ref = 0, has_alt = 0, ial = 0, jal = 0;
    #define BRANCH_INT(type_t,missing,vector_end) { \
        type_t *p = (type_t*) (fmt_ptr->p + isample*fmt_ptr->size); \
        for (i=0; i<fmt_ptr->n; i++) \
        { \
            if ( p[i] == vector_end ) break; /* smaller ploidy */ \
            if ( !p[i] || p[i] == missing ) continue; /* missing allele */ \
            int tmp = p[i]>>1; \
            if ( tmp>1 ) \
            { \
                if ( !ial ) { ial = tmp; has_alt = 1; } \
                else if ( tmp!=ial ) \
                { \
                    if ( tmp<ial ) \
                    { \
                        jal = ial; \
                        ial = tmp; \
                    } \
                    else \
                    { \
                        jal = tmp; \
                    } \
                    has_alt = 2; \
                } \
            } \
            else has_ref = 1; \
            nals++; \
        } \
    }
    switch (fmt_ptr->type) {
        case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
        default: fprintf(stderr, "[E::%s] todo: fmt_type %d\n", __func__, fmt_ptr->type); exit(1); break;
    }
    #undef BRANCH_INT

    if ( _ial ) *_ial = ial>0 ? ial-1 : ial;
    if ( _jal ) *_jal = jal>0 ? jal-1 : jal;
    if ( !nals ) return GT_UNKN;
    if ( nals==1 )
        return has_ref ? GT_HAPL_R : GT_HAPL_A;
    if ( !has_ref ) 
        return has_alt==1 ? GT_HOM_AA : GT_HET_AA;
    if ( !has_alt )
        return GT_HOM_RR;
    return GT_HET_RA;
}

int bcf_trim_alleles(const bcf_hdr_t *header, bcf1_t *line)
{
    int i;
    bcf_fmt_t *gt = bcf_get_fmt(header, line, "GT");
    if ( !gt ) return 0;

    int *ac = (int*) calloc(line->n_allele,sizeof(int));

    // check if all alleles are populated
    #define BRANCH(type_t,missing,vector_end) { \
        for (i=0; i<line->n_sample; i++) \
        { \
            type_t *p = (type_t*) (gt->p + i*gt->size); \
            int ial; \
            for (ial=0; ial<gt->size; ial++) \
            { \
                if ( p[ial]==vector_end ) break; /* smaller ploidy */ \
                if ( !(p[ial]>>1) || p[ial]==missing ) continue; /* missing allele */ \
                assert((p[ial]>>1)-1<line->n_allele); \
                ac[(p[ial]>>1)-1]++; \
            } \
        } \
    }
    switch (gt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
        default: fprintf(stderr, "[E::%s] todo: %d at %s:%d\n", __func__, gt->type, header->id[BCF_DT_CTG][line->rid].key, line->pos+1); exit(1); break;
    }
    #undef BRANCH

    int rm_als = 0, nrm = 0;
    for (i=1; i<line->n_allele; i++) 
    {
        if ( !ac[i] ) { rm_als |= 1<<i; nrm++; }
    }
    free(ac);

    if ( nrm ) bcf_remove_alleles(header, line, rm_als);
    return nrm;
}

void bcf_remove_alleles(const bcf_hdr_t *header, bcf1_t *line, int rm_mask)
{
    int *map = (int*) calloc(line->n_allele, sizeof(int));

    // create map of indexes from old to new ALT numbering and modify ALT
    kstring_t str = {0,0,0};
    kputs(line->d.allele[0], &str);

    int nrm = 0, i,j;  // i: ori alleles, j: new alleles
    for (i=1, j=1; i<line->n_allele; i++) 
    {
        if ( rm_mask & 1<<i )
        {
            // remove this allele
            line->d.allele[i] = NULL;
            nrm++;
            continue;
        }
        kputc(',', &str);
        kputs(line->d.allele[i], &str);
        map[i] = j;
        j++;
    }
    if ( !nrm ) { free(map); free(str.s); return; }

    int nR_ori = line->n_allele;
    int nR_new = line->n_allele-nrm;
    assert(nR_new > 0); // should not be able to remove reference allele
    int nA_ori = nR_ori-1;
    int nA_new = nR_new-1;

    int nG_ori = nR_ori*(nR_ori + 1)/2;
    int nG_new = nR_new*(nR_new + 1)/2;

    bcf_update_alleles_str(header, line, str.s);

    // remove from Number=G, Number=R and Number=A INFO fields.
    uint8_t *dat = NULL;
    int mdat = 0, ndat = 0, mdat_bytes = 0, nret;
    for (i=0; i<line->n_info; i++)
    {
        bcf_info_t *info = &line->d.info[i];
        int vlen = bcf_hdr_id2length(header,BCF_HL_INFO,info->key);
        
        if ( vlen!=BCF_VL_A && vlen!=BCF_VL_G && vlen!=BCF_VL_R ) continue; // no need to change

        int type = bcf_hdr_id2type(header,BCF_HL_INFO,info->key);
        if ( type==BCF_HT_FLAG ) continue;
        int size = 1;
        if ( type==BCF_HT_REAL || type==BCF_HT_INT ) size = 4;

        mdat = mdat_bytes / size;
        nret = bcf_get_info_values(header, line, bcf_hdr_int2id(header,BCF_DT_ID,info->key), (void**)&dat, &mdat, type);
        mdat_bytes = mdat * size;
        if ( nret<0 ) 
        { 
            fprintf(stderr,"[%s:%d %s] Could not access INFO/%s at %s:%d [%d]\n", __FILE__,__LINE__,__FUNCTION__, 
                bcf_hdr_int2id(header,BCF_DT_ID,info->key), bcf_seqname(header,line), line->pos+1, nret); 
            exit(1);
        }
        if ( type==BCF_HT_STR ) 
        { 
            str.l = 0;
            char *ss = (char*) dat, *se = (char*) dat;
            if ( vlen==BCF_VL_A || vlen==BCF_VL_R )
            {
                int nexp, inc = 0;
                if ( vlen==BCF_VL_A )
                {
                    nexp = nA_ori;
                    inc  = 1;
                }
                else
                    nexp = nR_ori;
                for (j=0; j<nexp; j++)
                {
                    if ( !*se ) break;
                    while ( *se && *se!=',' ) se++;
                    if ( rm_mask & 1<<(j+inc) ) 
                    { 
                        if ( *se ) se++;
                        ss = se; 
                        continue; 
                    }
                    if ( str.l ) kputc(',',&str);
                    kputsn(ss,se-ss,&str);
                    if ( *se ) se++;
                    ss = se;
                }
                assert( j==nexp );
            }
            else    // Number=G, assuming diploid genotype
            {
                int k = 0, n = 0;
                for (j=0; j<nR_ori; j++)
                {
                    for (k=0; k<=j; k++)
                    {
                        if ( !*se ) break;
                        while ( *se && *se!=',' ) se++;
                        n++;
                        if ( rm_mask & 1<<j || rm_mask & 1<<k ) 
                        { 
                            if ( *se ) se++;
                            ss = se; 
                            continue; 
                        }
                        if ( str.l ) kputc(',',&str);
                        kputsn(ss,se-ss,&str);
                        if ( *se ) se++;
                        ss = se;
                    }
                    if ( !*se ) break;
                }
                assert( n=nG_ori );
            }

            nret = bcf_update_info(header, line, bcf_hdr_int2id(header,BCF_DT_ID,info->key), (void*)str.s, str.l, type);
            if ( nret<0 )
            {
                fprintf(stderr,"[%s:%d %s] Could not update INFO/%s at %s:%d [%d]\n", __FILE__,__LINE__,__FUNCTION__,
                        bcf_hdr_int2id(header,BCF_DT_ID,info->key), bcf_seqname(header,line), line->pos+1, nret);
                exit(1);
            }
            continue; 
        }
        
        if ( vlen==BCF_VL_A || vlen==BCF_VL_R )
        {
            int inc = 0, ntop;
            if ( vlen==BCF_VL_A )
            {
                assert( nret==nA_ori );
                ntop = nA_ori;
                ndat = nA_new;
                inc  = 1;
            }
            else
            {
                assert( nret==nR_ori );
                ntop = nR_ori;
                ndat = nR_new;
            }
            int k = 0;

            #define BRANCH(type_t,is_vector_end) \
            { \
                type_t *ptr = (type_t*) dat; \
                int size = sizeof(type_t); \
                for (j=0; j<ntop; j++) /* j:ori, k:new */ \
                { \
                    if ( is_vector_end ) { memcpy(dat+k*size, dat+j*size, size); break; } \
                    if ( rm_mask & 1<<(j+inc) ) continue; \
                    if ( j!=k ) memcpy(dat+k*size, dat+j*size, size); \
                    k++; \
                } \
            }
            switch (type) 
            {
                case BCF_HT_INT:  BRANCH(int32_t,ptr[j]==bcf_int32_vector_end); break;
                case BCF_HT_REAL: BRANCH(float,bcf_float_is_vector_end(ptr[j])); break;
            }
            #undef BRANCH
        }
        else    // Number=G
        {
            assert( nret==nG_ori );
            int k, l_ori = -1, l_new = 0;
            ndat = nG_new;

            #define BRANCH(type_t,is_vector_end) \
            { \
                type_t *ptr = (type_t*) dat; \
                int size = sizeof(type_t); \
                for (j=0; j<nR_ori; j++) \
                { \
                    for (k=0; k<=j; k++) \
                    { \
                        l_ori++; \
                        if ( is_vector_end ) { memcpy(dat+l_new*size, dat+l_ori*size, size); break; } \
                        if ( rm_mask & 1<<j || rm_mask & 1<<k ) continue; \
                        if ( l_ori!=l_new ) memcpy(dat+l_new*size, dat+l_ori*size, size); \
                        l_new++; \
                    } \
                } \
            }
            switch (type) 
            {
                case BCF_HT_INT:  BRANCH(int32_t,ptr[l_ori]==bcf_int32_vector_end); break;
                case BCF_HT_REAL: BRANCH(float,bcf_float_is_vector_end(ptr[l_ori])); break;
            }
            #undef BRANCH
        }

        nret = bcf_update_info(header, line, bcf_hdr_int2id(header,BCF_DT_ID,info->key), (void*)dat, ndat, type);
        if ( nret<0 )
        {
            fprintf(stderr,"[%s:%d %s] Could not update INFO/%s at %s:%d [%d]\n", __FILE__,__LINE__,__FUNCTION__,
                    bcf_hdr_int2id(header,BCF_DT_ID,info->key), bcf_seqname(header,line), line->pos+1, nret);
            exit(1);
        }
    }

    // Update GT fields, the allele indexes might have changed
    for (i=1; i<line->n_allele; i++) if ( map[i]!=i ) break;
    if ( i<line->n_allele )
    {
        mdat = mdat_bytes / 4;  // sizeof(int32_t)
        nret = bcf_get_genotypes(header,line,(void**)&dat,&mdat);
        mdat_bytes = mdat * 4;
        if ( nret>0 )
        {
            nret /= line->n_sample;
            int32_t *ptr = (int32_t*) dat;
            for (i=0; i<line->n_sample; i++)
            {
                for (j=0; j<nret; j++)
                {
                    if ( ptr[j]==bcf_gt_missing ) continue;
                    if ( ptr[j]==bcf_int32_vector_end ) break;
                    int al = bcf_gt_allele(ptr[j]);
                    assert( al<nR_ori && map[al]>=0 );
                    ptr[j] = (map[al]+1)<<1 | (ptr[j]&1);
                }
                ptr += nret;
            }
            bcf_update_genotypes(header, line, (void*)dat, nret*line->n_sample);
        }
    }

    // Remove from Number=G, Number=R and Number=A FORMAT fields. 
    // Assuming haploid or diploid GTs
    for (i=0; i<line->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &line->d.fmt[i];
        int vlen = bcf_hdr_id2length(header,BCF_HL_FMT,fmt->id);

        if ( vlen!=BCF_VL_A && vlen!=BCF_VL_G && vlen!=BCF_VL_R ) continue; // no need to change

        int type = bcf_hdr_id2type(header,BCF_HL_FMT,fmt->id);
        if ( type==BCF_HT_FLAG ) continue;

        int size = 1;
        if ( type==BCF_HT_REAL || type==BCF_HT_INT ) size = 4;

        mdat = mdat_bytes / size;
        nret = bcf_get_format_values(header, line, bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), (void**)&dat, &mdat, type);
        mdat_bytes = mdat * size;
        if ( nret<0 ) 
        { 
            fprintf(stderr,"[%s:%d %s] Could not access FORMAT/%s at %s:%d [%d]\n", __FILE__,__LINE__,__FUNCTION__, 
                    bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), bcf_seqname(header,line), line->pos+1, nret); 
            exit(1);
        }

        if ( type==BCF_HT_STR ) 
        {
            int size = nret/line->n_sample;     // number of bytes per sample
            str.l = 0;
            if ( vlen==BCF_VL_A || vlen==BCF_VL_R )
            {
                int nexp, inc = 0;
                if ( vlen==BCF_VL_A )
                {
                    nexp = nA_ori;
                    inc  = 1;
                }
                else
                    nexp = nR_ori;
                for (j=0; j<line->n_sample; j++)
                {
                    char *ss = ((char*)dat) + j*size, *se = ss + size, *ptr = ss;
                    int k_src = 0, k_dst = 0, l = str.l;
                    for (k_src=0; k_src<nexp; k_src++)
                    {
                        if ( ptr>=se || !*ptr) break;
                        while ( ptr<se && *ptr && *ptr!=',' ) ptr++;
                        if ( rm_mask & 1<<(k_src+inc) )
                        {
                            ss = ++ptr;
                            continue;
                        }
                        if ( k_dst ) kputc(',',&str);
                        kputsn(ss,ptr-ss,&str);
                        ss = ++ptr;
                        k_dst++;
                    }
                    assert( k_src==nexp );
                    l = str.l - l;
                    for (; l<size; l++) kputc(0, &str);
                }
            }
            else    // Number=G, diploid or haploid
            {
                for (j=0; j<line->n_sample; j++)
                {
                    char *ss = ((char*)dat) + j*size, *se = ss + size, *ptr = ss;
                    int k_src = 0, k_dst = 0, l = str.l;
                    int nexp = 0; // diploid or haploid?
                    while ( ptr<se )
                    {
                        if ( !*ptr ) break;
                        if ( *ptr==',' ) nexp++;
                        ptr++;
                    }
                    if ( ptr!=ss ) nexp++;
                    assert( nexp==nG_ori || nexp==nR_ori );
                    ptr = ss;
                    if ( nexp==nG_ori ) // diploid
                    {
                        int ia, ib;
                        for (ia=0; ia<nR_ori; ia++)
                        {
                            for (ib=0; ib<=ia; ib++)
                            {
                                if ( ptr>=se || !*ptr ) break;
                                while ( ptr<se && *ptr && *ptr!=',' ) ptr++;
                                if ( rm_mask & 1<<ia || rm_mask & 1<<ib )
                                {
                                    ss = ++ptr;
                                    continue;
                                }
                                if ( k_dst ) kputc(',',&str);
                                kputsn(ss,ptr-ss,&str);
                                ss = ++ptr;
                                k_dst++;
                            }
                            if ( ptr>=se || !*ptr ) break;
                        }
                    }
                    else    // haploid
                    {
                        for (k_src=0; k_src<nR_ori; k_src++)
                        {
                            if ( ptr>=se || !*ptr ) break;
                            while ( ptr<se && *ptr && *ptr!=',' ) ptr++;
                            if ( rm_mask & 1<<k_src )
                            {
                                ss = ++ptr;
                                continue;
                            }
                            if ( k_dst ) kputc(',',&str);
                            kputsn(ss,ptr-ss,&str);
                            ss = ++ptr;
                            k_dst++;
                        }
                        assert( k_src==nR_ori );
                        l = str.l - l;
                        for (; l<size; l++) kputc(0, &str);
                    }
                }
            }
            nret = bcf_update_format(header, line, bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), (void*)str.s, str.l, type);
            if ( nret<0 )
            {
                fprintf(stderr,"[%s:%d %s] Could not update FORMAT/%s at %s:%d [%d]\n", __FILE__,__LINE__,__FUNCTION__,
                        bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), bcf_seqname(header,line), line->pos+1, nret);
                exit(1);
            }
            continue;
        }

        int nori = nret / line->n_sample;
        if ( vlen==BCF_VL_A || vlen==BCF_VL_R || (vlen==BCF_VL_G && nori==nR_ori) ) // Number=A, R or haploid Number=G
        {
            int ntop, inc = 0;
            if ( vlen==BCF_VL_A )
            {
                assert( nori==nA_ori );     // todo: will fail if all values are missing
                ntop = nA_ori;
                ndat = nA_new*line->n_sample;
                inc  = 1;
            }
            else
            {
                assert( nori==nR_ori );     // todo: will fail if all values are missing
                ntop = nR_ori;
                ndat = nR_new*line->n_sample;
            }

            #define BRANCH(type_t,is_vector_end) \
            { \
                for (j=0; j<line->n_sample; j++) \
                { \
                    type_t *ptr_src = ((type_t*)dat) + j*nori; \
                    type_t *ptr_dst = ((type_t*)dat) + j*nA_new; \
                    int size = sizeof(type_t); \
                    int k_src, k_dst = 0; \
                    for (k_src=0; k_src<ntop; k_src++) \
                    { \
                        if ( is_vector_end ) { memcpy(ptr_dst+k_dst, ptr_src+k_src, size); break; } \
                        if ( rm_mask & 1<<(k_src+inc) ) continue; \
                        if ( k_src!=k_dst ) memcpy(ptr_dst+k_dst, ptr_src+k_src, size); \
                        k_dst++; \
                    } \
                } \
            }
            switch (type) 
            {
                case BCF_HT_INT:  BRANCH(int32_t,ptr_src[k_src]==bcf_int32_vector_end); break;
                case BCF_HT_REAL: BRANCH(float,bcf_float_is_vector_end(ptr_src[k_src])); break;
            }
            #undef BRANCH
        }
        else    // Number=G, diploid or mixture of haploid+diploid
        {
            assert( nori==nG_ori );
            ndat = nG_new*line->n_sample;

            #define BRANCH(type_t,is_vector_end) \
            { \
                for (j=0; j<line->n_sample; j++) \
                { \
                    type_t *ptr_src = ((type_t*)dat) + j*nori; \
                    type_t *ptr_dst = ((type_t*)dat) + j*nG_new; \
                    int size = sizeof(type_t); \
                    int ia, ib, k_dst = 0, k_src; \
                    int nset = 0;   /* haploid or diploid? */ \
                    for (k_src=0; k_src<nG_ori; k_src++) { if ( is_vector_end ) break; nset++; } \
                    if ( nset==nR_ori ) /* haploid */ \
                    { \
                        for (k_src=0; k_src<nR_ori; k_src++) \
                        { \
                            if ( rm_mask & 1<<k_src ) continue; \
                            if ( k_src!=k_dst ) memcpy(ptr_dst+k_dst, ptr_src+k_src, size); \
                            k_dst++; \
                        } \
                        memcpy(ptr_dst+k_dst, ptr_src+k_src, size); \
                    } \
                    else /* diploid */ \
                    { \
                        k_src = -1; \
                        for (ia=0; ia<nR_ori; ia++) \
                        { \
                            for (ib=0; ib<=ia; ib++) \
                            { \
                                k_src++; \
                                if ( is_vector_end ) { memcpy(ptr_dst+k_dst, ptr_src+k_src, size); ia = nR_ori; break; } \
                                if ( rm_mask & 1<<ia || rm_mask & 1<<ib ) continue; \
                                if ( k_src!=k_dst ) memcpy(ptr_dst+k_dst, ptr_src+k_src, size); \
                                k_dst++; \
                            } \
                        } \
                    } \
                } \
            }
            switch (type) 
            {
                case BCF_HT_INT:  BRANCH(int32_t,ptr_src[k_src]==bcf_int32_vector_end); break;
                case BCF_HT_REAL: BRANCH(float,bcf_float_is_vector_end(ptr_src[k_src])); break;
            }
            #undef BRANCH
        }
        nret = bcf_update_format(header, line, bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), (void*)dat, ndat, type);
        if ( nret<0 )
        {
            fprintf(stderr,"[%s:%d %s] Could not update FORMAT/%s at %s:%d [%d]\n", __FILE__,__LINE__,__FUNCTION__,
                    bcf_hdr_int2id(header,BCF_DT_ID,fmt->id), bcf_seqname(header,line), line->pos+1, nret);
            exit(1);
        }
    }
    free(dat);
    free(str.s);
    free(map);
}

