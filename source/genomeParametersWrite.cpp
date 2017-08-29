#include "genomeParametersWrite.h"
#include "streamFuns.h"

void genomeParametersWrite(string fileName, Parameters& P, string errorOut)
{//write the genome information into the genomePar stream
    ofstream & genomePar = ofstrOpen(fileName, errorOut, P);

    genomePar << "### "<<P.commandLineFull <<"\n";

    genomePar << "versionGenome\t" << P.versionSTAR <<"\n";
    genomePar << "pGe.gFastaFiles\t";
    for (uint ii=0;ii<P.pGe.gFastaFiles.size();ii++) genomePar << P.pGe.gFastaFiles.at(ii) << " ";
    genomePar << "\n";
    genomePar << "pGe.gSAindexNbases\t" << P.pGe.gSAindexNbases << "\n";
    genomePar << "pGe.gChrBinNbits\t" << P.pGe.gChrBinNbits << "\n";
    genomePar << "pGe.gSAsparseD\t" << P.pGe.gSAsparseD <<"\n";
    genomePar << "pGe.sjdbOverhang\t" << P.pGe.sjdbOverhang <<"\n";

    genomePar << "pGe.sjdbFileChrStartEnd\t";
    for (uint ii=0;ii<P.pGe.sjdbFileChrStartEnd.size();ii++) genomePar<< P.pGe.sjdbFileChrStartEnd.at(ii) << " ";
    genomePar<<"\n";

    genomePar << "pGe.sjdbGTFfile\t" << P.pGe.sjdbGTFfile <<"\n";
    genomePar << "pGe.sjdbGTFchrPrefix\t" << P.pGe.sjdbGTFchrPrefix <<"\n";
    genomePar << "pGe.sjdbGTFfeatureExon\t" << P.pGe.sjdbGTFfeatureExon <<"\n";
    genomePar << "pGe.sjdbGTFtagExonParentTranscript\t" << P.pGe.sjdbGTFtagExonParentTranscript <<"\n";
    genomePar << "pGe.sjdbGTFtagExonParentGene\t" << P.pGe.sjdbGTFtagExonParentGene <<"\n";

    genomePar << "pGe.sjdbInsertSave\t" << P.pGe.sjdbInsertSave <<"\n";
    
    genomePar << "pGe.gFileSizes\t" << P.pGe.gFileSizes.at(0);
    for (uint ii=1;ii<P.pGe.gFileSizes.size();ii++) 
        genomePar << " " << P.pGe.gFileSizes.at(ii) ;
    genomePar << "\n";

    genomePar.close();
};
