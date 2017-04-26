#include "ReadAlign.h"
#include "OutSJ.h"

void ReadAlign::outputTranscriptSJ(Transcript const &trOut, uint nTrOut, OutSJ *chunkOutSJ, uint sjReadStartN ) {//record junctions in chunkOutSJ array

    //TODO: make sure that a junction is recorded onyl once from one read.
    //For a multimapper, several alignments may contain the same junctions - now it's recorded several time.
//     if (nTrOut>1) return; //junctions from multi-mappers are not recorded

//     if (P.outSAMmode=="None") return; //no SAM output

    for (uint iex=0;iex<trOut.nExons-1;iex++) {//record all junctions
        if (trOut.canonSJ[iex]>=0) {//only record junctions, not indels or mate gap
            chunkOutSJ->oneSJ.junctionPointer(chunkOutSJ->data, chunkOutSJ->N);//get pointer to an empty junction in the data array
            *chunkOutSJ->oneSJ.start=trOut.exons[iex][EX_G]+trOut.exons[iex][EX_L]; //start of the intron
            *chunkOutSJ->oneSJ.gap=trOut.exons[iex+1][EX_G]-*chunkOutSJ->oneSJ.start;
            //overhangs: basic method
            //*chunkOutSJ->oneSJ.overhangLeft  = (uint32) trOut.exons[iex][EX_L];//TODO calculate the lengh of overhangs taking into account indels
            //*chunkOutSJ->oneSJ.overhangRight = (uint32) trOut.exons[iex+1][EX_L];
            //overhangs: min method
            *chunkOutSJ->oneSJ.overhangLeft = min ( (uint32) trOut.exons[iex][EX_L],(uint32) trOut.exons[iex+1][EX_L] );
            *chunkOutSJ->oneSJ.overhangRight = *chunkOutSJ->oneSJ.overhangLeft;

            //check if this junction has been recorded from this read - this happens when the mates overlap and cross the same junctions
            bool duplicateSJ(false);
            for (uint ii=sjReadStartN; ii<chunkOutSJ->N; ii++) {//TODO if there are many junctions, need to make more efficient
                if ( *chunkOutSJ->oneSJ.start == *((uint*) (chunkOutSJ->data+ii*Junction::dataSize+Junction::startP)) \
                  && *chunkOutSJ->oneSJ.gap   == *((uint32*) (chunkOutSJ->data+ii*Junction::dataSize+Junction::gapP)) ) {
                    duplicateSJ=true;
                    uint16* overhang1=(uint16*) (chunkOutSJ->data+ii*Junction::dataSize+Junction::overhangLeftP);
                    if (*overhang1<*chunkOutSJ->oneSJ.overhangLeft) {
                        *overhang1=*chunkOutSJ->oneSJ.overhangLeft;
                        * ((uint16*) (chunkOutSJ->data+ii*Junction::dataSize+Junction::overhangRightP))=*overhang1;
                    };
                    break;
                };
            };
            if (duplicateSJ) continue; //do not record this junctions

            *chunkOutSJ->oneSJ.motif=trOut.canonSJ[iex];
            *chunkOutSJ->oneSJ.strand=(char) (trOut.canonSJ[iex]==0 ? 0 : (trOut.canonSJ[iex]+1)%2+1);
            *chunkOutSJ->oneSJ.annot=trOut.sjAnnot[iex];
            if (nTrOut==1) {
                *chunkOutSJ->oneSJ.countUnique=1;
                *chunkOutSJ->oneSJ.countMultiple=0;
            } else {
                *chunkOutSJ->oneSJ.countMultiple=1; //TODO: 1/nTrOut?
                *chunkOutSJ->oneSJ.countUnique=0; //TODO: 1/nTrOut?
            };

            chunkOutSJ->N++;//increment the number of recorded junctions
        };
    };
};
