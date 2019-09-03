#include "ChimericAlign.h"
#include "ReadAlign.h"

void ChimericAlign::chimericJunctionOutput(fstream &outStream, uint chimN, int maxNonChimAlignScore, bool PEmerged_flag, int chimScoreBest, int maxPossibleAlignScore)
{
    outStream << mapGen.chrName[al1->Chr] <<"\t"<< chimJ1 - mapGen.chrStart[al1->Chr]+1 <<"\t"<< (al1->Str==0 ? "+":"-") \
        <<"\t"<< mapGen.chrName[al2->Chr] <<"\t"<< chimJ2 - mapGen.chrStart[al2->Chr]+1 <<"\t"<< (al2->Str==0 ? "+":"-") \
        <<"\t"<< chimMotif <<"\t"<< chimRepeat1  <<"\t"<< chimRepeat2 <<"\t"<< al1->readName+1 \
        <<"\t"<< al1->exons[0][EX_G] - mapGen.chrStart[al1->Chr]+1 <<"\t"<< al1->generateCigarP() \
        <<"\t"<< al2->exons[0][EX_G] - mapGen.chrStart[al2->Chr]+1 <<"\t"<< al2->generateCigarP()
        <<"\t"<< chimN // number of multimapping chimeric alignments for this read.
        << "\t" << maxPossibleAlignScore  // the maximum possible alignment score (currently the sum of the (paired) read lengths)
        << "\t" << maxNonChimAlignScore // trBest - the best alignment score from a non-chimeric alignment of this read to the ref genome.
        << "\t" << chimScore    // current chimeric alignment score
        << "\t" << chimScoreBest  // best chimeric score among multimapping chimeric alignments.
        << "\t" << PEmerged_flag;  // boolean indicating paired reads were merged into a single read before alignment & chimer detection.

        if (P.outSAMattrPresent.RG)
            outStream <<"\t"<< P.outSAMattrRG.at(RA->readFilesIndex);
        if (P.pSolo.type>0)
            outStream <<"\t"<< RA->soloRead->readBar->cbSeq <<"\t"<< RA->soloRead->readBar->umiSeq;
        outStream <<"\n"; //<<"\t"<< trChim[0].exons[0][EX_iFrag]+1 --- no need for that, since trChim[0] is always on the first mate
};
