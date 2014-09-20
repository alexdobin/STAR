#include <unordered_map>
#include "bamRemoveDuplicates.h"
#include <iostream>
#include "htslib/htslib/sam.h"
#include "IncludeDefine.h"
#include SAMTOOLS_BGZF_H

#define compareReturn(a,b) if(a>b) {return 1;} else if (a<b) {return -1;}

uint g_bamRemoveDuplicatesMate2basesN;

int funCompareCoordName(const void *a, const void *b) {
    uint ca=*(uint*) a;
    uint cb=*(uint*) b;
    
    compareReturn(ca,cb) else {//compare read names
        uint32* pa=(uint32*) *(uint32**) (a+8);
        uint32* pb=(uint32*) *(uint32**) (b+8);
        
        uint32 la=(pa[3]<<24)>>24;
        uint32 lb=(pb[3]<<24)>>24;
        
        
        compareReturn(la,lb) else {
            char* ca=(char*) (pa+9);
            char* cb=(char*) (pb+9);
            for (int ii=0;ii<la;ii++) {
                compareReturn(ca[ii],cb[ii]);
            };
            uint32 fa=pa[4]>>16;
            uint32 fb=pb[4]>>16;

            compareReturn((fa&0x80), (fb&0x80));
            return 0;
        };
        
    };
};

int funCompareCoordFlagCigarSeq(const void *a, const void *b) {
    uint32* pa1=(uint32*) *(uint32**) a;
    uint32* pa2=(uint32*) *(uint32**) (a+8);
    
    uint32* pb1=(uint32*) *(uint32**) b;
    uint32* pb2=(uint32*) *(uint32**) (b+8);
    
    compareReturn(pa1[2],pb1[2]);//position match
    compareReturn(pa1[7],pb1[7]);//2nd mate position match
    compareReturn(pa1[8],pb1[8]);//insert length match
    compareReturn(pa1[4],pb1[4]);//FLAG and nCigar match
    compareReturn(pa2[4],pb2[4]);//FLAG and nCigar match - 2nd mate
    //compare CIGARs
    uint32 nc=(pa1[4]<<16)>>16;
    uint32* ciga=(uint32*) (((char*) pa1)+9*4+((pa1[3]<<24)>>24));
    uint32* cigb=(uint32*) (((char*) pb1)+9*4+((pb1[3]<<24)>>24));                        
    for (uint32 ii=0;ii<nc;ii++) {
        compareReturn(ciga[ii],cigb[ii]);
    };
    //compare CIGARs mate2
    nc=(pa2[4]<<16)>>16;
    ciga=(uint32*) (((char*) pa2)+9*4+((pa2[3]<<24)>>24));
    cigb=(uint32*) (((char*) pb2)+9*4+((pb2[3]<<24)>>24));                        
    for (uint32 ii=0;ii<nc;ii++) {
        compareReturn(ciga[ii],cigb[ii]);
    };                        
    //compare sequences
    uint8_t* sa=((uint8_t*) pa2)+9*4+((pa2[3]<<24)>>24)+nc*4;
    uint8_t* sb=((uint8_t*) pb2)+9*4+((pb2[3]<<24)>>24)+nc*4;
    if (((pa2[4]>>16) & 0x10) == 0) {//not reverse complemented
        uint ii=1;
        for (; ii<g_bamRemoveDuplicatesMate2basesN; ii+=2) {
            compareReturn(sa[ii/2],sb[ii/2]);
        };
        if (g_bamRemoveDuplicatesMate2basesN%2>0) {
            compareReturn((sa[ii/2]>>4),(sb[ii/2]>>4));
        };
    } else {
        int ii=pa2[5]-g_bamRemoveDuplicatesMate2basesN;
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

//int funCompareCoordFlagCigarSeqScore(const void *a, const void *b) {
//    int ret=funCompareCoordFlagCigarSeq(a,b);
//    if (ret==0) {//compare MAPQ
//        uint32 ma= ((*( ( (uint32*) *(uint32**) a ) + 3))<<16)>>24;
//        uint32 mb= ((*( ( (uint32*) *(uint32**) b ) + 3))<<16)>>24;
//
//
//     compareReturn(mb,ma);
//      return 0;
//    } else {
//        return ret;
//    };
//};

void bamRemoveDuplicates(const string bamFileName, const string bamFileNameOut, Parameters* const P) {
    g_bamRemoveDuplicatesMate2basesN=P->bamRemoveDuplicatesMate2basesN;

    bam1_t *bamA;
    bamA=bam_init1();
    
   
    BGZF *bamIn=bgzf_open(bamFileName.c_str(),"r");
    bam_hdr_t *bamHeader=bam_hdr_read(bamIn);
    
    BGZF *bgzfOut;
    bgzfOut=bgzf_open(bamFileNameOut.c_str(),("w"+to_string((long long) P->outBAMcompression)).c_str());
    bam_hdr_write(bgzfOut, bamHeader);
    
    uint bamLengthMax=P->limitBAMsortRAM/2; //max length to load
    if (bamLengthMax==0) bamLengthMax=100000000;
    uint grNmax=1000000;//max number of alignments

    char *bamRaw=new char[bamLengthMax];
    uint *aD=new uint[grNmax*2];
       
    uint bamLength=0;
    uint bamS=0, bamE=0, bamE1=0; //start/end/next-end position for read group search
    uint32 rightMax=0;
    uint grN=0;//number of reads in group
    bool bamFileEnd=false;//when the last alignment of rhe file was reached
    while (true) {
        if (bamE1>=bamLength) {//reached end of loaded BAM block, add BAM data
            if (bamLength<bamLengthMax && bamLength>0) {//reached end of BAM file, cannot load more
                bamFileEnd=true;
            } else {
                if (bamS==0 && bamLength>0) {//TODO
                    cerr << "need to increase bamLengthMax" <<endl;
                    exit(1);
                };

                //write out processed block
                bgzf_write(bgzfOut,bamRaw,bamS);
   
                bamLength-=bamS;
                memmove(bamRaw, bamRaw+bamS,bamLength); //move the non-processed part of the block to the beginning of bamRaw
                bamLength+=bgzf_read(bamIn, bamRaw+bamLength, bamLengthMax-bamLength);//marks the end of the BAM block that has been read
    //             cout <<"loaded block "<<bamLength<<endl;

    //             bamE-=bamS;
    //             bamE1-=bamS;            
                //restart search for the group
                bamS=0;
                bamE=0;
                bamE1=bamE+*(uint32*)(bamRaw+bamE)+4;//next alignment
                rightMax=0;
                grN=0;
            };
        };
   
        uint32* bamP=(uint32*) (bamRaw+bamE);//pointer to the 1st mate of the pair
        bamP[4] |= 0x400<<16;//mark all aligns as duplicate, will unmark          
        uint32 chrE=bamP[1];
        uint32 leftE=bamP[2];
        uint32 rightE=bamP[7];
        
        uint32 chrS=*(uint32*)(bamRaw+bamS+4*1);
                
        if ( chrE !=chrS ||  (rightMax>0 && leftE>rightMax) || bamFileEnd ) {//found new group of reads to be processed, call collapsing procedure

            qsort((void*) aD, grN, 2*sizeof(uint), funCompareCoordName);
//             cout <<"grN="<<grN<<endl;
            uint nP=0;
            for (uint iP=0; iP<grN/2; iP++) {
                //record pairs
                uint sP=aD[iP*4];//start of the pair, to be compared with the next one
                aD[nP*2]=aD[iP*4+1];
                aD[nP*2+1]=aD[iP*4+3];
                ++nP;                
                if (iP+1==grN/2 || aD[iP*4+4]!=sP) {//next position is different from the previous, run collapsing
                    if (nP>0) {//otherwise only one pair at this position, no duplicates
                        qsort((void*) aD, nP, 2*sizeof(uint), funCompareCoordFlagCigarSeq);
                        //go through the list and select non-duplicates
                        int bScore=-999, bP=0;
                        for (uint pp=0; pp<nP; pp++) {
                            uint32* bamP1=(uint32*) aD[pp*2];//pointer to the 1st mate of the pair                           
                            bamA->data=((uint8_t*) bamP1)+9*4+((bamP1[3]<<24)>>24)+((bamP1[4]<<16)>>16)*4+(bamP1[5]+1)/2+bamP1[5];//add length for: core, name, cigar, seq, qual
                            bamA->l_data=((uint8_t*) bamP1)+bamP1[0]+1-bamA->data;
                            int score1=bam_aux2i(bam_aux_get(bamA,"AS"));
                            if (score1>bScore) {
                                bScore=score1;
                                bP=pp;
                            };
                            
                            if ( pp+1==nP || funCompareCoordFlagCigarSeq((void*) (aD+pp*2),(void*) (aD+pp*2+2))!=0 ) {//next pair is not equal to the current one
                                //un-mark duplicates
                                uint32* bamPb=(uint32*) aD[bP*2+1];//pointer to the 2nd mate of the pair
                                bamPb[4] ^= 0x400<<16;
                                bamPb=(uint32*) aD[bP*2];//pointer to the 1st mate of the pair
                                bamPb[4] ^= 0x400<<16;   
                                cout << ((char*)(bamPb+9)) <<"\n";
                                bScore=-999;//reset best score
                            };                            
                        };
                    };
                    nP=0;
                };
            };
            
            //reset for the next group
            if (bamFileEnd) break; //exit the main cycle over blocks
            rightMax=0;
            bamS=bamE;
            grN=0;
        };
        
        bamA->data=((uint8_t*) bamP)+9*4+((bamP[3]<<24)>>24)+((bamP[4]<<16)>>16)*4+(bamP[5]+1)/2+bamP[5];//add length for: core, name, cigar, seq, qual
        bamA->l_data=((uint8_t*) bamP)+bamP[0]+1-bamA->data;
        int nMult=bam_aux2i(bam_aux_get(bamA,"NH"));
        
        if (nMult==1) {//record this alignment in the current group, unique mappers only. Multi-mappers will not be considered for collapsing, and will remain marked as duplicates
            if (rightE>leftE) {//left mate, record coordinate of its right mate
                rightMax=max(rightMax, rightE);
                aD[grN*2]=leftE;
            } else {
                aD[grN*2]=rightE;
            };        
            aD[grN*2+1]=(uint) bamRaw+bamE; 
            ++grN;
            if (grN>=grNmax) {//reallocate
                grNmax=grN*2;
                uint *aD1=new uint[grNmax*2];     
                memcpy((char*) aD1, (char*) aD, grN*2*sizeof(uint));
                delete [] aD;
                aD=aD1;
                cerr << "reallocated array "<<grNmax<<endl;                
            };            
        };
        
        bamE=bamE1;//shift to the next record
        bamE1=bamE+*(uint32*)(bamRaw+bamE)+4;//next alignment

    };

    bgzf_write(bgzfOut,bamRaw,bamLength);
    bgzf_flush(bgzfOut);
    bgzf_close(bgzfOut);
};
