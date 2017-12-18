#include "ChimericAlign.h"
#include "ReadAlign.h"

void ChimericAlign::chimericJunctionOutput(fstream &outStream, uint chimN)
{
    outStream << mapGen.chrName[al1->Chr] <<"\t"<< chimJ1 - mapGen.chrStart[al1->Chr]+1 <<"\t"<< (al1->Str==0 ? "+":"-") \
        <<"\t"<< mapGen.chrName[al2->Chr] <<"\t"<< chimJ2 - mapGen.chrStart[al2->Chr]+1 <<"\t"<< (al2->Str==0 ? "+":"-") \
        <<"\t"<< chimMotif <<"\t"<< chimRepeat1  <<"\t"<< chimRepeat2 <<"\t"<< al1->readName+1 \
        <<"\t"<< al1->exons[0][EX_G] - mapGen.chrStart[al1->Chr]+1 <<"\t"<< al1->generateCigarP() \
        <<"\t"<< al2->exons[0][EX_G] - mapGen.chrStart[al2->Chr]+1 <<"\t"<< al2->generateCigarP() <<"\t"<< chimN; 
        if (P.outSAMattrPresent.RG)
            outStream <<"\t"<< P.outSAMattrRG.at(RA->readFilesIndex);
        outStream <<"\n"; //<<"\t"<< trChim[0].exons[0][EX_iFrag]+1 --- no need for that, since trChim[0] is always on the first mate
};
        