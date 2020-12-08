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
// no need to realloc b->data, because it is overwritten later with bamChar
//		b->data = (uint8_t*)realloc(b->data, b->m_data);
//		if (!b->data)
//			return -4;
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

// calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open)
int reg2bin(int beg, int end)
{
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
    return 0;
};

int bamAttrArrayWrite(int32 attr, const char* tagName, char* attrArray ) {
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    attrArray[2]='i';
    *( (int32*) (attrArray+3))=attr;
    return 3+sizeof(int32);
};
int bamAttrArrayWrite(float attr, const char* tagName, char* attrArray ) {
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    attrArray[2]='f';
    *( (float*) (attrArray+3))=attr;
    return 3+sizeof(int32);
};
int bamAttrArrayWrite(char attr, const char* tagName, char* attrArray ) {
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    attrArray[2]='A';
    attrArray[3]=attr;
    return 3+sizeof(char);
};
int bamAttrArrayWrite(string &attr, const char* tagName, char* attrArray ) {
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    attrArray[2]='Z';
    memcpy(attrArray+3,attr.c_str(),attr.size()+1);//copy string data including \0
    return 3+attr.size()+1;
};
int bamAttrArrayWrite(const vector<char> &attr, const char* tagName, char* attrArray ) {
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    attrArray[2]='B';
    attrArray[3]='c';
    *( (int32*) (attrArray+4))=attr.size();
    memcpy(attrArray+4+sizeof(int32),attr.data(),attr.size());//copy array data
    return 4+sizeof(int32)+attr.size();
};
int bamAttrArrayWrite(const vector<int32> &attr, const char* tagName, char* attrArray ) {
    attrArray[0]=tagName[0];attrArray[1]=tagName[1];
    attrArray[2]='B';
    attrArray[3]='i';
    *( (int32*) (attrArray+4))=attr.size();
    memcpy(attrArray+4+sizeof(int32),attr.data(),sizeof(int32)*attr.size());//copy array data
    return 4+sizeof(int32)+sizeof(int32)*attr.size();
};

int bamAttrArrayWriteSAMtags(string &attrStr, char *attrArray, Parameters &P) {//write bam record into attrArray for string attribute attString
    size_t pos1=0, pos2=0;
    int nattr=0;
    do {//cycle over multiple tags separated by tab
        pos2 = attrStr.find('\t',pos1);
        string attr1 = attrStr.substr(pos1, pos2-pos1); //substring containing one tag
        pos1=pos2+1; //for the next search

        if (attr1.empty())
            continue; //extra tab at the beginning, or consecutive tabs
            
        uint16_t tagn = * ( (uint16_t*) attr1.c_str() );
        if ( !P.readFiles.samAttrKeepAll && P.readFiles.samAttrKeep.count(tagn)==0 )
            continue; //skip tags not contained the list            

        switch (attr1.at(3)) {//TODO: support other tag types (vector tags)
            case 'i':
            {
                int32 a1=stol(attr1.substr(5));
                nattr += bamAttrArrayWrite(a1,attr1.c_str(),attrArray+nattr);
                break;
            };
            case 'A':
            {
                char a1=attr1.at(5);
                nattr += bamAttrArrayWrite(a1,attr1.c_str(),attrArray+nattr);
                break;
            };
                break;
            case 'Z':
            {
                string a1=attr1.substr(5);
                nattr += bamAttrArrayWrite(a1,attr1.c_str(),attrArray+nattr);
                break;
            };
            case 'f':
            {
                float a1=stof(attr1.substr(5));
                nattr += bamAttrArrayWrite(a1,attr1.c_str(),attrArray+nattr);
                break;
            };
        };
    } while (pos2 != string::npos);

    return nattr;
};


