#include "ChimericAlign.h"

void ChimericAlign::chimericJunctionOutput(fstream &outStream, uint chimN)
{
    outStream << P.chrName[al1->Chr] <<"\t"<< chimJ1 - P.chrStart[al1->Chr]+1 <<"\t"<< (al1->Str==0 ? "+":"-") \
        <<"\t"<< P.chrName[al2->Chr] <<"\t"<< chimJ2 - P.chrStart[al2->Chr]+1 <<"\t"<< (al2->Str==0 ? "+":"-") \
        <<"\t"<< chimMotif <<"\t"<< chimRepeat1  <<"\t"<< chimRepeat2 <<"\t"<< al1->readName+1 \
        <<"\t"<< al1->exons[0][EX_G] - P.chrStart[al1->Chr]+1 <<"\t"<< al1->generateCigarP() \
        <<"\t"<< al2->exons[0][EX_G] - P.chrStart[al2->Chr]+1 <<"\t"<< al2->generateCigarP() <<"\t"<< chimN <<"\n"; 
};
        