#include "ReadAlign.h"
#include "BAMfunctions.h"
#include "ChimericAlign.h"

void ReadAlign::chimericDetectionOldOutput() {

    if (!chimRecord) {
        return;
    };

    //re-calculate the score for chimeric transcripts
    trChim[0].alignScore(Read1, mapGen.G, P);
    trChim[1].alignScore(Read1, mapGen.G, P);

    if (P.pCh.out.bam) //BAM output
        ChimericAlign::chimericBAMoutput(&trChim[0], &trChim[1], this, 0, 1, true, P);

    if (P.pCh.out.samOld) {

        chimN=2; //this  is hard-coded for now
        if (trChim[0].exons[0][EX_iFrag]!=trChim[0].exons[trChim[0].nExons-1][EX_iFrag]) {//tr0 has both mates
            trChim[0].primaryFlag=true;//paired portion is primary
            trChim[1].primaryFlag=false;
        } else if (trChim[1].exons[0][EX_iFrag]!=trChim[1].exons[trChim[1].nExons-1][EX_iFrag]) {//tr1 has both mates
            trChim[1].primaryFlag=true;//paired portion is primary
            trChim[0].primaryFlag=false;
        } else if (trChim[0].exons[0][EX_iFrag]!=trChim[1].exons[0][EX_iFrag]) {//tr0 and tr1 are single different mates
            trChim[0].primaryFlag=true;
            trChim[1].primaryFlag=true;
        } else  {//two chimeric segments are on the same mate - this can only happen for single-end reads
            int chimRepresent = (trChim[0].maxScore > trChim[1].maxScore) ? 0 : 1;
            trChim[chimRepresent].primaryFlag=true;
            trChim[1-chimRepresent].primaryFlag=false;
        };

        for (uint iTr=0;iTr<chimN;iTr++)
        {//write all chimeric pieces to Chimeric.out.sam/junction
            if (P.readNmates==2) {//PE: need mate info //not readNends: this is alignment
                uint iex=0;
                if ( trChim[1-iTr].exons[0][EX_iFrag] != trChim[1-iTr].exons[trChim[1-iTr].nExons-1][EX_iFrag] )
                {//the other segment has 2 mates, need to find the opposite mate
                    for (;iex<trChim[1-iTr].nExons;iex++) {
                        if (trChim[1-iTr].exons[iex][EX_iFrag]!=trChim[iTr].exons[0][EX_iFrag]) {
                            break;
                        };
                    };
                };

                uint mateChr=trChim[1-iTr].Chr;
                uint mateStart=trChim[1-iTr].exons[iex][EX_G];
                char mateStrand=(char) (trChim[1-iTr].Str!=trChim[1-iTr].exons[iex][EX_iFrag]);

                outputTranscriptSAM(trChim[iTr], chimN, iTr, mateChr, mateStart, mateStrand, -1, NULL, &chunkOutChimSAM);
            } else
            {
                outputTranscriptSAM(trChim[iTr], chimN, iTr, -1, -1, -1, -1, NULL, &chunkOutChimSAM);
            };
        };
    };

    if (P.pCh.out.junctions) {
        //junction + SAMp
        *chunkOutChimJunction << mapGen.chrName[trChim[0].Chr] <<"\t"<< chimJ0 - mapGen.chrStart[trChim[0].Chr]+1 <<"\t"<< (trChim[0].Str==0 ? "+":"-") \
                <<"\t"<< mapGen.chrName[trChim[1].Chr] <<"\t"<< chimJ1 - mapGen.chrStart[trChim[1].Chr]+1 <<"\t"<< (trChim[1].Str==0 ? "+":"-") \
                <<"\t"<< chimMotif <<"\t"<< chimRepeat0  <<"\t"<< chimRepeat1 <<"\t"<< readName+1 \
                <<"\t"<< trChim[0].exons[0][EX_G] - mapGen.chrStart[trChim[0].Chr]+1 <<"\t"<< outputTranscriptCIGARp(trChim[0]) \
                <<"\t"<< trChim[1].exons[0][EX_G] - mapGen.chrStart[trChim[1].Chr]+1 <<"\t"<<  outputTranscriptCIGARp(trChim[1]);
        if (P.outSAMattrPresent.RG)
            *chunkOutChimJunction <<"\t"<< P.outSAMattrRG.at(readFilesIndex);
        *chunkOutChimJunction <<"\n"; //<<"\t"<< trChim[0].exons[0][EX_iFrag]+1 --- no need for that, since trChim[0] is always on the first mate
    };

    return;
};
