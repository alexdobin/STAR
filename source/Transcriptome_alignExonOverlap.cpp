#include "Transcriptome.h"
#include "ReadAlign.h"
#include "serviceFuns.cpp"
#include "AlignVsTranscript.h"
#include "ReadAnnotations.h"
#include <bitset>

int32 alignBlocksOverlapExons(Transcript &aG, uint16 exN1, uint32 *exSE1, uint64 trStart1, bool &sjConcord);

void Transcriptome::alignExonOverlap(uint nA, Transcript **aAll, int32 strandType, ReadAnnotFeature &annFeat)
{
    struct GeneStrOverlapAlign {
        uint32 g;
        int32 ov;
        uint32 ia;
        uint32 exl;
        bool str;
        bool sjc;
    };
    //vector<GeneStrOverlapAlign> vGeStrOvAl;
    //vGeStrOvAl.reserve(256); //TODO: check if this affects speed


    typedef array<bool,6> OverlapTypes;
    struct GeneInfo1 {
        uint32 g;
        uint32 ia;
        OverlapTypes ot; //overlap types, prioritized
    };
    vector<GeneInfo1> vGeneInfo1;
    vGeneInfo1.reserve(256); //TODO: check if this affects speed
    OverlapTypes otAS={false,true,false,true,false,true}; //which OverlapTypes is antisense, it will not be counted

    /*//GeneFullClosest3p
    uint64 minDist3p=(uint64)-1;
    uint32 minDist3pGene=(uint32)-1;
    */

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
                 
            bool str1 = int(strandType==0 ? aG.Str : 1-aG.Str) == (trStr[tr1]-1);
            str1 = str1 || (strandType==-1);

            /*
            //GeneFullClosest3p
            if ( str1 ) {
                uint64 dist3p = ( aG.Str==0 ? trE[tr1]-aGend : aGstart-trS[tr1] );
                if ( dist3p < minDist3p ) {
                    minDist3p = dist3p;
                    minDist3pGene = trGene[tr1];
                    continue;
                };
            };
            */

            bool sjConcord;
            int32 nOverlap = alignBlocksOverlapExons(aG, trExN[tr1], exSE+2*trExI[tr1], trS[tr1], sjConcord);
            if ( nOverlap>=0 ) {

                int exl=0;
                for (uint64 iex=0; iex<aG.nExons;iex++)
                    exl += aG.exons[iex][EX_L];

                //vGeStrOvAl.push_back({trGene[tr1], nOverlap, iag, exl, str1, sjConcord});
                
                vGeneInfo1.push_back({trGene[tr1], iag, {
                                                            str1 && nOverlap==exl && sjConcord,
                                                            !str1 && nOverlap==exl && sjConcord,
                                                            str1 && nOverlap>exl/2, 
                                                            !str1 && nOverlap>exl/2,
                                                            str1,
                                                            !str1
                                                        } 
                                     });

                /* matches Gene
                vGeneInfo1.push_back({trGene[tr1], iag, {
                                                            str1 && nOverlap==exl && sjConcord, 
                                                            false, 
                                                            false, 
                                                            false
                                                        } 
                                     });
                */

                /* matches GeneFull_ExonOverIntron
                vGeneInfo1.push_back({trGene[tr1], iag, {
                                                            str1 && nOverlap==exl && sjConcord, 
                                                            false, 
                                                            false, 
                                                            str1
                                                        } 
                                     });                                     
                */

                /*

                bool sjann=true;
                for (uint64 iex=0; iex<aG.nExons-1;iex++)
                    sjann = sjann & (aG.canonSJ[iex]<0 || aG.sjAnnot[iex]!=0);

                if ( str1 && nOverlap==exl && sjConcord ) {//this should match geneConcordant
                    annFeat.fSet.insert(trGene[tr1]);
                };
                */               
            };

            //cout << trGene[tr1] <<" "<< nOverlap << " " <<flush;
        } while (trEmax[tr1]>=aGend && tr1>0);
    };

    /*
    {//GeneFullClosest3p
        annFeat.fSet={};
        annFeat.fAlign = {};
        annFeat.fAlign.resize(nA);
        if (nA==1 && minDist3pGene!=(uint32)-1) {
            annFeat.fSet.insert(minDist3pGene);
            annFeat.fAlign[0].insert(minDist3pGene);
        };
        return;
    };
    */

    {//prune geneFull_Ex50pAS according to priorities
        /*
        bool bTranscr = false; //transcriptomic
        bool bExonic = false; //sense, exonic: >=50% exonic
        bool bExonicAS = false; 
        bool bIntronic = false; //sense, intronic
        for (auto &v1: vGeStrOvAl) {
            if ( v1.str ) {//sense
                if (v1.ov==v1.exl && v1.sjc) {
                    bTranscr = true;
                    break; //this is the highest priority
                } else if (v1.ov>=v1.exl/2) {
                    bExonic = true;
                } else {
                    bIntronic = true;
                };
            } else if (v1.ov>=v1.exl/2) {
                bExonicAS = true;
            };
        };
        

        for (auto &v1: vGeneInfo1) {
            if ( v1.tr ) {//sense
                bTranscr = true;
                break; //this is the highest priority
            } else if ( v1.ex ) {
                bExonic = true;
            } else if ( v1.as ) {
                bExonicAS = true;
            } else if (v1.in) {
                bIntronic = true;
            };
        };
        if (bTranscr) {
            for (auto &v1: vGeneInfo1)
                if ( v1.tr )
                    annFeat.fSet.insert(v1.g);
        } else if (bExonic) {
            for (auto &v1: vGeneInfo1)
                if ( v1.ex )
                    annFeat.fSet.insert(v1.g);
        } else if (bIntronic) {
            for (auto &v1: vGeneInfo1)
                if ( v1.in )
                    annFeat.fSet.insert(v1.g);
        };
        */

        OverlapTypes otFinal={};
        for ( auto &v1: vGeneInfo1 ) {
            for ( uint32 it=0; it<otFinal.size(); it++ ) {
                if ( v1.ot[it] ) {
                    otFinal[it] = true;
                    break; //lower it are prioritized over higher
                };
           };
        };

        if (otFinal[0]) {
            annFeat.ovType = ReadAnnotFeature::overlapTypes::exonic;
        } else if (otFinal[1]) {
            annFeat.ovType = ReadAnnotFeature::overlapTypes::exonicAS;
        } else if (otFinal[2]) {
            annFeat.ovType = ReadAnnotFeature::overlapTypes::exonic50p;
        } else if (otFinal[3]) {
            annFeat.ovType = ReadAnnotFeature::overlapTypes::exonic50pAS;
        } else if (otFinal[4]) {
            annFeat.ovType = ReadAnnotFeature::overlapTypes::intronic;
        } else if (otFinal[5]) {
            annFeat.ovType = ReadAnnotFeature::overlapTypes::intronicAS;
        } else {
            annFeat.ovType = ReadAnnotFeature::overlapTypes::intergenic; //intergenic. i.e. no overlap with genes
        };

        // annFeat.fSet={};
        // annFeat.fAlign = {};
        annFeat.fAlign.resize(nA);
        for ( uint32 it=0; it<otFinal.size(); it++ ) {
            if ( otFinal[it] ) {
                if (otAS[it])
                    return; //AS reads are not counted TODO: add to stats
                for ( auto &v1: vGeneInfo1 ) {
                    if ( v1.ot[it] ) {
                        annFeat.fSet.insert(v1.g);
                        annFeat.fAlign[v1.ia].insert(v1.g);
                    };
                };
                break;//first otFinal wins
            };
        };      
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