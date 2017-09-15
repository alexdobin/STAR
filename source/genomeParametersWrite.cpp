#include "genomeParametersWrite.h"
#include "streamFuns.h"

void genomeParametersWrite(string fileName, Parameters& P, string errorOut)
{//write the genome information into the genomePar stream
    ofstream & genomePar = ofstrOpen(fileName, errorOut, P);

    genomePar << "### "<<P.commandLineFull <<"\n";

    genomePar << "versionGenome\t" << P.versionSTAR <<"\n";
    genomePar << "genomeFastaFiles\t";
    for (uint ii=0;ii<P.pGe.gFastaFiles.size();ii++) genomePar << P.pGe.gFastaFiles.at(ii) << " ";
    genomePar << "\n";
    genomePar << "genomeSAindexNbases\t" << P.pGe.gSAindexNbases << "\n";
    genomePar << "genomeChrBinNbits\t" << P.pGe.gChrBinNbits << "\n";
    genomePar << "genomeSAsparseD\t" << P.pGe.gSAsparseD <<"\n";
    genomePar << "sjdbOverhang\t" << mapGen.sjdbOverhang <<"\n";

    genomePar << "sjdbFileChrStartEnd\t";
    for (uint ii=0;ii<P.pGe.sjdbFileChrStartEnd.size();ii++) genomePar<< P.pGe.sjdbFileChrStartEnd.at(ii) << " ";
    genomePar<<"\n";

    genomePar << "sjdbGTFfile\t" << P.pGe.sjdbGTFfile <<"\n";
    genomePar << "sjdbGTFchrPrefix\t" << P.pGe.sjdbGTFchrPrefix <<"\n";
    genomePar << "sjdbGTFfeatureExon\t" << P.pGe.sjdbGTFfeatureExon <<"\n";
    genomePar << "sjdbGTFtagExonParentTranscript\t" << P.pGe.sjdbGTFtagExonParentTranscript <<"\n";
    genomePar << "sjdbGTFtagExonParentGene\t" << P.pGe.sjdbGTFtagExonParentGene <<"\n";

    genomePar << "sjdbInsertSave\t" << P.pGe.sjdbInsertSave <<"\n";
    
    genomePar << "genomeFileSizes\t" << P.pGe.gFileSizes.at(0);
    for (uint ii=1;ii<P.pGe.gFileSizes.size();ii++) 
        genomePar << " " << P.pGe.gFileSizes.at(ii) ;
    genomePar << "\n";

    genomePar.close();
};
