#include "BAMfunctions.h"
#include "htslib/htslib/kstring.h"

string bam_cigarString (bam1_t *b) {//output CIGAR string
//    kstring_t strK;
//    kstring_t *str=&strK;
    const bam1_core_t *c = &b->core;

    string cigarString("");
    if ( c->n_cigar > 0 ) {
      uint32_t *cigar = bam_get_cigar(b);
      for (int i = 0; i < c->n_cigar; ++i) {
        cigarString+=to_string((uint)bam_cigar_oplen(cigar[i]))+bam_cigar_opchr(cigar[i]);
      };
    };


//	if (c->n_cigar) { // cigar
//		for (int i = 0; i < c->n_cigar; ++i) {
//			kputw(bam_cigar_oplen(cigar[i]), str);
//			kputc(bam_cigar_opchr(cigar[i]), str);
//		}
//	} else kputc('*', str);
//
//    string cigarString (str->s,str->l);
    return cigarString;
};

int bam_read1_fromArray(char *bamChar, bam1_t *b) //modified from samtools bam_read1 to assign BAM record in mmemry to bam structure
{
	bam1_core_t *c = &b->core;
	int32_t block_len; //, ret, i;
// // 	uint32_t x[8];
// // 	if ((ret = bgzf_read(fp, &block_len, 4)) != 4) {
// // 		if (ret == 0) return -1; // normal end-of-file
// // 		else return -2; // truncated
// // 	}
 	uint32_t *x;

    uint32_t *bamU32=(uint32_t*) bamChar;
    block_len=bamU32[0];

// // 	if (bgzf_read(fp, x, 32) != 32) return -3;
// // 	if (fp->is_be) {
// // 		ed_swap_4p(&block_len);
// // 		for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
// // 	}
    x=bamU32+1;

	c->tid = x[0]; c->pos = x[1];
	c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
	c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
	c->l_qseq = x[4];
	c->mtid = x[5]; c->mpos = x[6]; c->isize = x[7];
	b->l_data = block_len - 32;
	if (b->l_data < 0 || c->l_qseq < 0) return -4;
	if ((char *)bam_get_aux(b) - (char *)b->data > b->l_data)
		return -4;
	if (b->m_data < b->l_data) {
		b->m_data = b->l_data;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
		if (!b->data)
			return -4;
	}
// // 	if (bgzf_read(fp, b->data, b->l_data) != b->l_data) return -4;
// // 	//b->l_aux = b->l_data - c->n_cigar * 4 - c->l_qname - c->l_qseq - (c->l_qseq+1)/2;
// // 	if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    b->data=(uint8_t*) bamChar+4*9;

	return 4 + block_len;
}


void outBAMwriteHeader (BGZF* fp, const string &samh, const vector <string> &chrn, const vector <uint> &chrl) {
    bgzf_write(fp,"BAM\001",4);
    int32 hlen=samh.size();
    bgzf_write(fp,(char*) &hlen,sizeof(hlen));
    bgzf_write(fp,samh.c_str(),hlen);
    int32 nchr=(int32) chrn.size();
    bgzf_write(fp,(char*) &nchr,sizeof(nchr));
    for (int32 ii=0;ii<nchr;ii++) {
        int32 rlen = (int32) (chrn.at(ii).size()+1);
        int32 slen = (int32) chrl[ii];
        bgzf_write(fp,(char*) &rlen,sizeof(rlen));
        bgzf_write(fp,chrn.at(ii).data(),rlen); //this includes \0 at the end of the string
        bgzf_write(fp,(char*) &slen,sizeof(slen));
    };
    bgzf_flush(fp);
};

template <class TintType>
TintType bamAttributeInt(const char *bamAux, const char *attrName) {//not tested!!!
    const char *attrStart=strstr(bamAux,attrName);
    if (attrStart==NULL) return (TintType) -1;
    switch (attrStart[2]) {
        case ('c'):
            return (TintType) *(int8_t*)(attrStart+3);
        case ('s'):
            return (TintType) *(int16_t*)(attrStart+3);
        case ('i'):
            return (TintType) *(int32_t*)(attrStart+3);
        case ('C'):
            return (TintType) *(uint8_t*)(attrStart+3);
        case ('S'):
            return (TintType) *(uint16_t*)(attrStart+3);
        case ('I'):
            return (TintType) *(uint32_t*)(attrStart+3);
    };
};
