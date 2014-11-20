#include "Transcriptome.h"
#include "ReadAlign.h"
#include "Transcript.h"
#include "serviceFuns.cpp"

uint ReadAlign::quantTranscriptome (Transcriptome *Tr, uint nAlignG, Transcript **alignG, Transcript *alignT) {
    uint nAlignT=0;    
    for (uint iag=0; iag<nAlignG; iag++) {//transform all alignments
        if (alignG[iag]->nDel>0 || alignG[iag]->nIns>0) continue; //prevent indels
        uint nMM1=0;
        char* R=Read1[alignG[iag]->roStr==0 ? 0:2];
        for (uint32 iab=0; iab<alignG[iag]->nExons; iab++) {//check for soft-clips and indels
            uint left1=0,right1=0;//how many bases to move left or right
            if (iab==0) {
                left1=alignG[iag]->exons[iab][EX_R];
            } else if (alignG[iag]->canonSJ[iab-1]==-3) {
                left1=alignG[iag]->exons[iab][EX_R]-readLength[alignG[iag]->exons[iab-1][EX_iFrag]]-1;
            };
            if (alignG[iag]->canonSJ[iab]==-3) {
                right1=readLength[alignG[iag]->exons[iab][EX_iFrag]]-alignG[iag]->exons[iab][EX_R]-alignG[iag]->exons[iab][EX_L];
            } else if (iab==alignG[iag]->nExons-1){
                right1=Lread-alignG[iag]->exons[iab][EX_R]-alignG[iag]->exons[iab][EX_L];
            };

            for (uint b=1; b<=left1 ; b++) {//extend to the left
                char r1=R[alignG[iag]->exons[iab][EX_R]-b];
                char g1=G[alignG[iag]->exons[iab][EX_G]-b];
                if ( r1!=g1 && r1<4 && g1<4) ++nMM1;
            };
            for (uint b=0; b<right1 ; b++) {//extend to the left
                char r1=R[alignG[iag]->exons[iab][EX_R]+alignG[iag]->exons[iab][EX_L]+b];
                char g1=G[alignG[iag]->exons[iab][EX_G]+alignG[iag]->exons[iab][EX_L]+b];
                if ( r1!=g1 && r1<4 && g1<4) ++nMM1;
            };            
            alignG[iag]->exons[iab][EX_R] -= left1;
            alignG[iag]->exons[iab][EX_G] -= left1;
            alignG[iag]->exons[iab][EX_L] += left1+right1;
        };
        if ( (alignG[iag]->nMM + nMM1) > min(outFilterMismatchNmaxTotal, (uint) (P->outFilterMismatchNoverLmax*(Lread-1)) ) ) {
            continue;
        };

        
//         if (alignG[iag]->mappedLength==(readLength[0]+readLength[1]) && alignG[iag]->nDel==0) {//remove transcripts that contain indels of soft-clipping //TODO make this optional
        nAlignT += Tr->quantAlign(*alignG[iag],alignT+nAlignT);
//         };
    };
    
    for (uint iatr=0;iatr<nAlignT;iatr++) {//write all transcripts
//         alignBAM(alignT[iatr], nAlignT, iatr, 0, (uint) -1, (uint) -1, 0, -1, NULL);
//         outBAMarray1+=bamBytes;
//         outBAMbytes1+=bamBytes;
        alignBAM(alignT[iatr], nAlignT, iatr, 0, (uint) -1, (uint) -1, 0, -1, NULL, P->outSAMattrOrderQuant);        
        for (uint imate=0; imate<P->readNmates; imate++) {//output each mate
            outBAMquant->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1]);
        };        
    };
    
    if (nAlignT==0 && P->outSAMunmapped=="Within") {//read could not be mapped to transcriptome
        uint unmapType=5;
        bool mateMapped[2]={false,false};
        alignBAM(*alignG[0], 0, 0, P->chrStart[alignG[0]->Chr], (uint) -1, (uint) -1, 0,  unmapType, mateMapped, P->outSAMattrOrder);
            for (uint imate=0; imate<P->readNmates; imate++) {//output each mate
                outBAMquant->unsortedOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1]);
            };

    };
    return nAlignT;    
};
