#include <unordered_map>
#include "bamRemoveDuplicates.h"
#include <iostream>
#include "htslib/htslib/sam.h"
#include "IncludeDefine.h"
#include SAMTOOLS_BGZF_H
#include "ErrorWarning.h"

#define compareReturn(a,b) if(a>b) {return 1;} else if (a<b) {return -1;}

uint g_bamRemoveDuplicatesMate2basesN;

int funCompareNames(const void *a, const void *b) {//compare read names
    uint32* pa=(uint32*) *(uint32**) a;
    uint32* pb=(uint32*) *(uint32**) b;

    uint32 la=(pa[3]<<24)>>24;
    uint32 lb=(pb[3]<<24)>>24;

    compareReturn(la,lb) else {
        char* ca=(char*) (pa+9);
        char* cb=(char*) (pb+9);
        for (uint32 ii=0;ii<la;ii++) {
            compareReturn(ca[ii],cb[ii]);
        };
        uint32 fa=pa[4]>>16;
        uint32 fb=pb[4]>>16;

        compareReturn((fa&0x80), (fb&0x80));
        return 0;
    };
};

uint32 funStartExtendS(const uint32* const p) {//calculates align start extending right S operation
    uint32* cig=(uint32*) (((char*) p)+9*4+((p[3]<<24)>>24));
    if ( ((cig[0]<<28)>>28) == 4 ) {//first (right) operation is S
        return p[2]-(cig[0]>>4);
    } else {
        return p[2];
    };
};

uint32 funCigarExtendS(const uint32* const p, uint32* cout) {
    uint32* cig=(uint32*) (((char*) p)+9*4+((p[3]<<24)>>24));
    uint32 n=(p[4]<<16)>>16, n1=n;

    if (((cig[0]<<28)>>28) == 4) {
        --n1;
        memcpy((char*) cout, (char*) (cig+1), n1*sizeof(uint32));//copy CIGAR starting from the 2nd operation
        cout[0]+=(cig[0]>>4)<<4;
    } else {
        memcpy((char*) cout, (char*) cig, n*sizeof(uint32));//copy full CIGAR
    };
    if (((cig[n-1]<<28)>>28) == 4) {//remove last S opeartion add length to previous M
        --n1;
        cout[n1-1]+=(cig[n-1]>>4)<<4;
    };
    return n1;
};

int funCompareCigarsExtendS(const uint32* const pa, const uint32* const pb){
    uint32 ca[100], cb[100];
    uint32 na=funCigarExtendS(pa,ca);
    uint32 nb=funCigarExtendS(pb,cb);
    compareReturn(na,nb);
    for (uint32 ii=0; ii<na; ii++) {
        compareReturn(ca[ii],cb[ii]);
    };
    return 0;
};

int funCompareCoordFlagCigarSeq(const void *a, const void *b) {
    uint32* pa1=(uint32*) *(uint32**) a;
    uint32* pa2=(uint32*) *(uint32**) ((char*)a+8);

    uint32* pb1=(uint32*) *(uint32**) b;
    uint32* pb2=(uint32*) *(uint32**) ((char*)b+8);

    compareReturn(funStartExtendS(pa1),funStartExtendS(pb1));//position match
    compareReturn(funStartExtendS(pa2),funStartExtendS(pb2));//2nd mate position match
    compareReturn(pa1[4]>>16,pb1[4]>>16);//FLAG match
    compareReturn(pa2[4]>>16,pb2[4]>>16);//FLAG match - 2nd mate

    int ret1=funCompareCigarsExtendS(pa1,pb1);
    if (ret1!=0) return ret1;
    ret1=funCompareCigarsExtendS(pa2,pb2);
    if (ret1!=0) return ret1;

    //compare sequences
    uint8_t* sa=((uint8_t*) pa2)+9*4+((pa2[3]<<24)>>24)+((pa2[4]<<16)>>16)*4;
    uint8_t* sb=((uint8_t*) pb2)+9*4+((pb2[3]<<24)>>24)+((pb2[4]<<16)>>16)*4;
    if (((pa2[4]>>16) & 0x10) == 0) {//not reverse complemented
        uint ii=1;
        for (; ii<g_bamRemoveDuplicatesMate2basesN; ii+=2) {
            compareReturn(sa[ii/2],sb[ii/2]);
        };
        if (g_bamRemoveDuplicatesMate2basesN%2>0) {
            compareReturn((sa[ii/2]>>4),(sb[ii/2]>>4));
        };
    } else {
        uint32 ii=pa2[5]-g_bamRemoveDuplicatesMate2basesN;
        if (ii%2>0) {
            compareReturn((sa[ii/2]&15),(sb[ii/2]&15));
            ++ii;
        };
        for (; ii<pa2[5]; ii+=2) {
            compareReturn(sa[ii/2],sb[ii/2]);
        };
    };

    return 0;
};

void bamRemoveDuplicates(const string bamFileName, const string bamFileNameOut, Parameters &P) {
    g_bamRemoveDuplicatesMate2basesN=P.removeDuplicates.mate2basesN;

    bam1_t *bamA;
    bamA=bam_init1();

    BGZF *bamIn=bgzf_open(bamFileName.c_str(),"r");
    bam_hdr_t *bamHeader=bam_hdr_read(bamIn);

    BGZF *bgzfOut;
    bgzfOut=bgzf_open(bamFileNameOut.c_str(),("w"+to_string((long long) P.outBAMcompression)).c_str());
    bam_hdr_write(bgzfOut, bamHeader);

    uint bamLengthMax=P.limitBAMsortRAM; //max length to load
    if (bamLengthMax==0) bamLengthMax=16000000000;
    uint grNmax=1000000;//max number of alignments

    char *bamRaw=new char[bamLengthMax];
    uint *aD=new uint[grNmax*2];

    uint bamLength=0;
    uint bamS=0, bamE=0, bamE1=1; //start/end/next-end position for read group search
    uint32 rightMax=0;
    uint grN=0;//number of reads in group
    bool bamFileEnd=false;//when the last alignment of rhe file was reached
    while (true) {
        if (bamE1>bamLength) {//reached end of loaded BAM block, add BAM data
            if (bamLength<bamLengthMax && bamLength>0) {//reached end of BAM file, cannot load more
                bamFileEnd=true;
            } else {
                if (bamS==0 && bamLength>0) {//TODO
                    ostringstream errOut;
                    errOut <<"EXITING because of fatal ERROR: not enough memory for marking duplicates \n";
                    errOut <<"SOLUTION: re-run STAR with at least --limitBAMsortRAM " <<bamLengthMax*2;
                    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
                };

                //write out processed block
                bgzf_write(bgzfOut,bamRaw,bamS);

                bamLength-=bamS;
                memmove(bamRaw, bamRaw+bamS,bamLength); //move the non-processed part of the block to the beginning of bamRaw
                bamLength+=bgzf_read(bamIn, bamRaw+bamLength, bamLengthMax-bamLength);//marks the end of the BAM block that has been read
                //restart search for the group
                bamS=0;
                bamE=0;
                bamE1=bamE+*(uint32*)(bamRaw+bamE)+4;//next alignment
                rightMax=0;
                grN=0;
            };
        };

        int nMult=0;
        uint32 chrE=0;
        uint32 leftE=0;
        uint32 rightE=0;
        uint32 chrS=0;        
        
        if (!bamFileEnd)
        {
            uint32* bamP=(uint32*) (bamRaw+bamE);//pointer to the 1st mate of the pair
            
            bamA->data=((uint8_t*) bamP)+9*4+((bamP[3]<<24)>>24)+((bamP[4]<<16)>>16)*4+(bamP[5]+1)/2+bamP[5];//add length for: core, name, cigar, seq, qual
            bamA->l_data=((uint8_t*) bamP)+bamP[0]+1-bamA->data;
                
            nMult=bam_aux2i(bam_aux_get(bamA,"NH"));
        
            if (nMult==1 || (nMult>1 && P.removeDuplicates.markMulti)) 
            {
                bamP[4] |= (0x400<<16);//mark all aligns as duplicate, will unmark. If multimappers, onyl mark if markMult=true
            };
            
            chrE=bamP[1];
            leftE=bamP[2];
            rightE=bamP[7];

            chrS=*(uint32*)(bamRaw+bamS+4*1);            
        };

        if ( chrE !=chrS ||  (rightMax>0 && leftE>rightMax) || bamFileEnd ) {//found new group of reads to be processed, start collapsing procedure
            qsort((void*) aD, grN, sizeof(uint), funCompareNames);
            qsort((void*) aD, grN/2, 2*sizeof(uint), funCompareCoordFlagCigarSeq);
            //go through the list and select non-duplicates
            int bScore=-999, bP=0;
            for (uint pp=0; pp<grN/2; pp++) {
                uint32* bamP1=(uint32*) aD[pp*2];//pointer to the 1st mate of the pair
                bamA->data=((uint8_t*) bamP1)+9*4+((bamP1[3]<<24)>>24)+((bamP1[4]<<16)>>16)*4+(bamP1[5]+1)/2+bamP1[5];//add length for: core, name, cigar, seq, qual
                bamA->l_data=((uint8_t*) bamP1)+bamP1[0]+1-bamA->data;
                int score1=bam_aux2i(bam_aux_get(bamA,"AS"));
                if (score1>bScore) {
                    bScore=score1;
                    bP=pp;
                };

                if ( pp==(grN/2-1) || funCompareCoordFlagCigarSeq((void*) (aD+pp*2),(void*) (aD+pp*2+2))!=0 ) {//next pair is not equal to the current one
                    //un-mark duplicates
                    uint32* bamPb=(uint32*) aD[bP*2+1];//pointer to the 2nd mate of the pair
                    bamPb[4] ^= (0x400<<16);
                    bamPb=(uint32*) aD[bP*2];//pointer to the 1st mate of the pair
                    bamPb[4] ^= (0x400<<16);
                    //cout << ((char*)(bamPb+9)) <<"\n";
                    bScore=-999;//reset best score
                };
            };


            //reset for the next group
            if (bamFileEnd) break; //exit the main cycle over blocks
            rightMax=0;
            bamS=bamE;
            grN=0;
        };

        if (nMult==1) {//record this alignment in the current group, unique mappers only. Multi-mappers will not be considered for collapsing, and will remain marked as duplicates
            if (grN>=grNmax) {//reallocate
                grNmax=grN*2;
                uint *aD1=new uint[grNmax];
                memcpy((char*) aD1, (char*) aD, grN*sizeof(uint));
                delete [] aD;
                aD=aD1;
                cerr << "reallocated array "<<grNmax<<endl;
            };
            aD[grN]=(uint) bamRaw+bamE;
            ++grN;
            if (rightE>leftE) {//left mate, record coordinate of its right mate
                rightMax=max(rightMax, rightE);
            };
        };

        bamE=bamE1;//shift to the next record
        bamE1=bamE+*(uint32*)(bamRaw+bamE)+4;//next alignment

    };

    bgzf_write(bgzfOut,bamRaw,bamLength);
    bgzf_flush(bgzfOut);
    bgzf_close(bgzfOut);
};
