#include "Transcriptome.h"
#include "ReadAlign.h"
#include "Transcript.h"
#include "serviceFuns.cpp"
#include <random>

uint ReadAlign::quantTranscriptome (Transcriptome *Tr, uint nAlignG, Transcript **alignG, Transcript *alignT) {
    uint nAlignT=0;
    for (uint iag=0; iag<nAlignG; iag++) {//transform all alignments

        if (!P.quant.trSAM.indel && (alignG[iag]->nDel>0 || alignG[iag]->nIns>0) )
        {//prevent indels if requested
            continue;
        };
        if (!P.quant.trSAM.singleEnd && (P.readNmates==2 && alignG[iag]->exons[0][EX_iFrag]==alignG[iag]->exons[alignG[iag]->nExons-1][EX_iFrag]) )
        {//prevent single end alignments
            continue;
        };

        uint nMM1=0;
        char* R=Read1[alignG[iag]->roStr==0 ? 0:2];
        if (!P.quant.trSAM.softClip)
        {//soft clipping not allowed, extend them if possible
            for (uint32 iab=0; iab<alignG[iag]->nExons; iab++) {
                uint left1=0,right1=0;//how many bases to move left or right
                if (iab==0) {
                    left1=alignG[iag]->exons[iab][EX_R];
                } else if (alignG[iag]->canonSJ[iab-1]==-3) {
                    left1=alignG[iag]->exons[iab][EX_R]-readLength[alignG[iag]->exons[iab-1][EX_iFrag]]-1;
                };
                if (iab==alignG[iag]->nExons-1)
                {//last block of left mates
                    right1=Lread-alignG[iag]->exons[iab][EX_R]-alignG[iag]->exons[iab][EX_L];

                } else if (alignG[iag]->canonSJ[iab]==-3)
                {//last block of the right mate (i.e. whole read)
                    right1=readLength[alignG[iag]->exons[iab][EX_iFrag]]-alignG[iag]->exons[iab][EX_R]-alignG[iag]->exons[iab][EX_L];
                };

                for (uint b=1; b<=left1 ; b++) {//extend to the left
                    char r1=R[alignG[iag]->exons[iab][EX_R]-b];
                    char g1=mapGen.G[alignG[iag]->exons[iab][EX_G]-b];
                    if ( r1!=g1 && r1<4 && g1<4) ++nMM1;
                };
                for (uint b=0; b<right1 ; b++) {//extend to the left
                    char r1=R[alignG[iag]->exons[iab][EX_R]+alignG[iag]->exons[iab][EX_L]+b];
                    char g1=mapGen.G[alignG[iag]->exons[iab][EX_G]+alignG[iag]->exons[iab][EX_L]+b];
                    if ( r1!=g1 && r1<4 && g1<4) ++nMM1;
                };
                alignG[iag]->exons[iab][EX_R] -= left1;
                alignG[iag]->exons[iab][EX_G] -= left1;
                alignG[iag]->exons[iab][EX_L] += left1+right1;
            };

            if ( (alignG[iag]->nMM + nMM1) > min(outFilterMismatchNmaxTotal, (uint) (P.outFilterMismatchNoverLmax*(Lread-1)) ) )
            {//extension of soft clips yielded too many mismatches, no output
                continue;
            };
        };

//         if (alignG[iag]->mappedLength==(readLength[0]+readLength[1]) && alignG[iag]->nDel==0) {//remove transcripts that contain indels of soft-clipping //TODO make this optional
        nAlignT += Tr->quantAlign(*alignG[iag],alignT+nAlignT);
//         };
    };

    alignT[int(rngUniformReal0to1(rngMultOrder)*nAlignT)].primaryFlag=true;

    for (uint iatr=0;iatr<nAlignT;iatr++) {//write all transcripts
//         alignBAM(alignT[iatr], nAlignT, iatr, 0, (uint) -1, (uint) -1, 0, -1, NULL, outBAMoneAlign, outBAMoneAlignNbytes);
//         outBAMarray1+=bamBytes;
//         outBAMbytes1+=bamBytes;
        alignBAM(alignT[iatr], nAlignT, iatr, 0, (uint) -1, (uint) -1, 0, -1, NULL, P.outSAMattrOrderQuant, outBAMoneAlign, outBAMoneAlignNbytes);
        for (uint imate=0; imate<P.readNmates; imate++) {//output each mate
            outBAMquant->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], (imate>0 || iatr>0) ? 0 : (outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1])*2*nAlignT);
        };
    };

    //not used anymore, at Colin Dewey's request
//     if (nAlignT==0 && P.outSAMunmapped=="Within") {//read could be mapped to genome, but not transcriptome - output as unmapped
//         uint unmapType=5;
//         bool mateMapped[2]={false,false};
//         alignBAM(*alignG[0], 0, 0, mapGen.chrStart[alignG[0]->Chr], (uint) -1, (uint) -1, 0,  unmapType, mateMapped, P.outSAMattrOrder);
//             for (uint imate=0; imate<P.readNmates; imate++) {//output each mate
//                 outBAMquant->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1]);
//             };
//
//     };

    return nAlignT;
};
