#include "Transcriptome.h"
#include "ReadAlign.h"
#include "serviceFuns.cpp"
#include "AlignVsTranscript.h"
#include "ReadAnnotations.h"
#include <bitset>

int32 alignBlocksOverlapExons(Transcript &aG, uint16 exN1, uint32 *exSE1, uint64 trStart1, bool &sjConcord);

void Transcriptome::alignExonOverlap(uint nA, Transcript **aAll, int32 strandType, ReadAnnotations &readAnnot)
{
    readAnnot.geneFull_CR={};
    readAnnot.geneFull_CR_Tr={};

    struct GeneStrOverlapAlign {
        uint32 g;
        int32 ov;
        uint32 ia;
        bool str;
    };

    vector<GeneStrOverlapAlign> vGeStrOvAl;
    vGeStrOvAl.reserve(256); //TODO: check if this affects speed

    for (uint32 iag=0; iag<nA; iag++) {
        
        Transcript &aG=*aAll[iag];

        //TODO: this only works if PE mates do not protrude. It's better to find the min-start max-end genomic coordinates for all exons.
        uint64 aGstart=aG.exons[0][EX_G];
        uint64 aGend=aG.exons[aG.nExons-1][EX_G]+aG.exons[aG.nExons-1][EX_L]-1;

        //binary search through transcript starts
        uint32 tr1=binarySearch1a<uint>(aGstart, trS, nTr);//tr1 has the maximum transcript start such that it is still <= align start: aGstart>=trS[tr1]
        if (tr1==(uint32) -1) 
            continue; //this alignment is outside of range of all transcripts

        ++tr1;
        do {//cycle back through all the transcripts
            --tr1;
            if ( aGend>trE[tr1] ) //alignment end is outside of this transcript
                continue;
                 
            bool sjConcord;
            int32 nOverlap = alignBlocksOverlapExons(aG, trExN[tr1], exSE+2*trExI[tr1], trS[tr1], sjConcord);
            if ( nOverlap>=0 ) {
                bool str1 = (strandType==0 ? aG.Str : 1-aG.Str) == (trStr[tr1]-1);
                str1 = str1 || (strandType==-1);
                //vGeStrOvAl.push_back({trGene[tr1], nOverlap, iag, str1});


                uint64 exl=0;
                for (uint64 iex=0; iex<aG.nExons;iex++)
                    exl += aG.exons[iex][EX_L];
                bool sjann=true;
                for (uint64 iex=0; iex<aG.nExons-1;iex++)
                    sjann = sjann & (aG.canonSJ[iex]<0 || aG.sjAnnot[iex]!=0);

                if ( str1 && nOverlap==exl && sjConcord) {//this should match geneConcordant
                    readAnnot.geneFull_CR.insert(trGene[tr1]);
                };

            };

            //cout << trGene[tr1] <<" "<< nOverlap << " " <<flush;
        } while (trEmax[tr1]>=aGend && tr1>0);
    };

    if (readAnnot.geneFull_CR != readAnnot.geneConcordant ){
        cout << 0;
    };
};

int32 alignBlocksOverlapExons(Transcript &aG, uint16 exN1, uint32 *exSE1, uint64 trStart1, bool &sjConcord)
{//calculate overlap between blocks of align and transcript

    uint32 i1=0, i2=0;
    int32 nOverlap = 0;
    sjConcord = true;
    //trStart1 += exSE1[0]; //not really needed since exSE1[0]=0
    uint64 trEnd1 = trStart1 + exSE1[2*exN1-1] + 1;//1st base after the end

    while (i1<aG.nExons && i2<exN1) {//scan through all blocks and exons
        uint64 rs1 = aG.exons[i1][EX_G];
        uint64 re1 = aG.exons[i1][EX_G]+aG.exons[i1][EX_L];//1st base after the end

        uint64 rs2 = trStart1 + exSE1[2*i2];
        uint64 re2 = trStart1 + exSE1[2*i2+1] + 1;

        if (rs1 < trStart1 || re1 > trEnd1) //this can happen for PE reads, when the 2nd mate protrudes to the left of the first or vice versa
            return -1;

        if (rs1>=re2) {//t1 block is on the right to t2, no hope of overlap
            i2++;
            if (i1>0 && aG.canonSJ[i1-1]>=0)
                sjConcord = false; //only 1st align block does not have to match exon start
        } else if (rs2>=re1) {//t2 block is on the right to t1, no hope of overlap
            i1++;
            sjConcord = false;
        } else {//overlap
            nOverlap += min(re1,re2) - max(rs1,rs2);

            if (i1>0 && rs1!=rs2 && aG.canonSJ[i1-1]>=0 )
                sjConcord = false;
            
            if (i1<aG.nExons-1 && re1!=re2 && aG.canonSJ[i1]>=0 )
                sjConcord = false;                

            if (re1>=re2) 
                i2++;//1 is on the right of 2
            if (re2>=re1) 
                i1++;//2 is on the right of 1
        };
    };


    return nOverlap;
};