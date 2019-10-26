#include "ReadAlign.h"
#include "ErrorWarning.h"

void ReadAlign::transformGenome() 
{//convert to new genome
    if (!mapGen.genomeOut.convYes || mapGen.pGe.transform.type==0 || nTr > P.outFilterMultimapNmax || nTr==0)
        return;
    
    uint32 nTr1=0;
    for (uint32 iTr=0; iTr<nTr; iTr++) {//convert output transcripts into new genome
        
        trMult[iTr]->haploType = ( trMult[iTr]->Chr >= mapGen.nChrReal/2 ? 1 : 2 );
        
        *trMultOut[nTr1]=*trMult[iTr];//copy information before conversion
        if (trMult[iTr]->transformGenome(*mapGen.genomeOut.g, *trMultOut[nTr1])) {            
            ++nTr1;
        };
    };
    
    if (mapGen.pGe.transform.type==2) {//remove identical transcripts
        bool keepTr[nTr1];
        for (uint32 ia1=0; ia1<nTr1; ia1++)
            keepTr[ia1]=true;
        
        for (uint32 ia1=0; ia1<nTr1; ia1++) {
            if (!keepTr[ia1])
                continue;
            for (uint32 ia2=ia1+1; ia2<nTr1; ia2++) {//compare to all previous ones
                
                if (!keepTr[ia1])
                    continue;
                
                Transcript &a1 = *trMultOut[ia1];
                Transcript &a2 = *trMultOut[ia2];
                
                if ( a1.Chr==a2.Chr && a1.Str==a2.Str 
                      && a1.exons[0][EX_G]-a1.exons[0][EX_R] == a2.exons[0][EX_G]-a2.exons[0][EX_R]
                      && a1.exons[a1.nExons][EX_G]+a1.exons[a1.nExons][EX_L]-a1.exons[a1.nExons][EX_R] == a2.exons[a2.nExons][EX_G]+a2.exons[a2.nExons][EX_L]-a2.exons[a2.nExons][EX_R] 
                    ) {//matching ends, remove one of the transcripts
                    trMultOut[ia1]->haploType = 0; //undefined haplotype
                    trMultOut[ia2]->haploType = 0;
                    if (a1.maxScore>a2.maxScore) {
                        keepTr[ia2]=false;
                    } else {
                        keepTr[ia1]=false;
                    };
                };
            };
        };
        
        nTr=0;
        for (uint32 ia1=0; ia1<nTr1; ia1++) {
            if (keepTr[ia1]) {
                trMult[nTr] = trMultOut[ia1];
                nTr++;
            };
        };
            
    } else {
        
        for (uint32 iTr=0; iTr<nTr1; iTr++)
            trMult[iTr] = trMultOut[iTr]; //point to new transcsript

        nTr=nTr1;
    };
};
