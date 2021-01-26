#include "Transcriptome.h"
#include "ReadAlign.h"
#include "Transcript.h"
#include "serviceFuns.cpp"
#include <random>

uint ReadAlign::quantTranscriptome (Transcriptome *Tr, uint nAlignG, Transcript **alignG, Transcript *alignT) {
    uint nAlignT=0;
    for (uint iag=0; iag<nAlignG; iag++) {//transform all alignments

        Transcript *align1=alignG[iag];

        if (!P.quant.trSAM.indel && (align1->nDel>0 || align1->nIns>0) ) {
            //prevent indels if requested
            continue;
        };
        if (!P.quant.trSAM.singleEnd && (P.readNmates==2 && align1->exons[0][EX_iFrag]==align1->exons[align1->nExons-1][EX_iFrag]) ) {//not readNends: this is alignment
        //prevent single end alignments
            continue;
        };

        if (!P.quant.trSAM.softClip) {
            //soft clipping not allowed, extend them if possible
            uint nMM1=0;
            char* R=Read1[align1->roStr==0 ? 0:2];
            Transcript align2=*align1; //copy this transcript to avoid changing the original one
            
            for (uint32 iab=0; iab<align2.nExons; iab++) {
                uint left1=0,right1=0;//how many bases to move left or right
                if (iab==0) {
                    left1=align2.exons[iab][EX_R];
                } else if (align2.canonSJ[iab-1]==-3) {
                    left1=align2.exons[iab][EX_R]-readLength[align2.exons[iab-1][EX_iFrag]]-1;
                };
                if (iab==align2.nExons-1) {//last block of left mates
                    right1=Lread-align2.exons[iab][EX_R]-align2.exons[iab][EX_L];

                } else if (align2.canonSJ[iab]==-3) {//last block of the right mate (i.e. whole read)
                    right1=readLength[align2.exons[iab][EX_iFrag]]-align2.exons[iab][EX_R]-align2.exons[iab][EX_L];
                };

                for (uint b=1; b<=left1 ; b++) {//extend to the left
                    char r1=R[align2.exons[iab][EX_R]-b];
                    char g1=mapGen.G[align2.exons[iab][EX_G]-b];
                    if ( r1!=g1 && r1<4 && g1<4) ++nMM1;
                };
                for (uint b=0; b<right1 ; b++) {//extend to the left
                    char r1=R[align2.exons[iab][EX_R]+align2.exons[iab][EX_L]+b];
                    char g1=mapGen.G[align2.exons[iab][EX_G]+align2.exons[iab][EX_L]+b];
                    if ( r1!=g1 && r1<4 && g1<4) ++nMM1;
                };
                align2.exons[iab][EX_R] -= left1;
                align2.exons[iab][EX_G] -= left1;
                align2.exons[iab][EX_L] += left1+right1;
            };

            if ( (align2.nMM + nMM1) > min(outFilterMismatchNmaxTotal, (uint) (P.outFilterMismatchNoverLmax*(Lread-1)) ) ) {
                //extension of soft clips yielded too many mismatches, no output
                continue;
            };
            
            align1 = &align2;
        };

        nAlignT += Tr->quantAlign(*align1,alignT+nAlignT);
    };

    if (P.quant.trSAM.bamYes) {//output Aligned.toTranscriptome.bam
        alignT[int(rngUniformReal0to1(rngMultOrder)*nAlignT)].primaryFlag=true;

        for (uint iatr=0;iatr<nAlignT;iatr++) {//write all transcripts
            alignBAM(alignT[iatr], nAlignT, iatr, 0, (uint) -1, (uint) -1, 0, -1, NULL, P.outSAMattrOrderQuant, outBAMoneAlign, outBAMoneAlignNbytes);
            for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
                outBAMquant->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], (imate>0 || iatr>0) ? 0 : (outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1])*2*nAlignT);
            };
        };
    };

    //not used anymore per Colin Dewey's request
    //     if (nAlignT==0 && P.outSAMunmapped=="Within") {//read could be mapped to genome, but not transcriptome - output as unmapped
    //         uint unmapType=5;
    //         bool mateMapped[2]={false,false};
    //         alignBAM(*alignG[0], 0, 0, mapGen.chrStart[alignG[0]->Chr], (uint) -1, (uint) -1, 0,  unmapType, mateMapped, P.outSAMattrOrder);
    //             for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
    //                 outBAMquant->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1]);
    //             };
    //
    //     };

    return nAlignT;
};
