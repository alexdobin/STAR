#include "genomeParametersWrite.h"
#include "streamFuns.h"

void genomeParametersWrite(string fileName, Parameters& P, string errorOut, Genome &mapGen)
{//write the genome information into the genomePar stream
    ofstream & genomePar = ofstrOpen(fileName, errorOut, P);

    genomePar << "### "<<P.commandLineFull <<"\n";
    genomePar << "### GstrandBit"<< mapGen.GstrandBit <<"\n";

    genomePar << "versionGenome\t" << P.versionSTAR <<"\n";
    genomePar << "genomeFastaFiles\t";
    for (uint ii=0;ii<mapGen.pGe.gFastaFiles.size();ii++) genomePar << mapGen.pGe.gFastaFiles.at(ii) << " ";
    genomePar << "\n";
    genomePar << "genomeSAindexNbases\t" << mapGen.pGe.gSAindexNbases << "\n";
    genomePar << "genomeChrBinNbits\t" << mapGen.pGe.gChrBinNbits << "\n";
    genomePar << "genomeSAsparseD\t" << mapGen.pGe.gSAsparseD <<"\n";
    genomePar << "sjdbOverhang\t" << mapGen.sjdbOverhang <<"\n";

    genomePar << "sjdbFileChrStartEnd\t";
    for (uint ii=0;ii<mapGen.pGe.sjdbFileChrStartEnd.size();ii++) genomePar<< mapGen.pGe.sjdbFileChrStartEnd.at(ii) << " ";
    genomePar<<"\n";

    genomePar << "sjdbGTFfile\t" << mapGen.pGe.sjdbGTFfile <<"\n";
    genomePar << "sjdbGTFchrPrefix\t" << mapGen.pGe.sjdbGTFchrPrefix <<"\n";
    genomePar << "sjdbGTFfeatureExon\t" << mapGen.pGe.sjdbGTFfeatureExon <<"\n";
    genomePar << "sjdbGTFtagExonParentTranscript\t" << mapGen.pGe.sjdbGTFtagExonParentTranscript <<"\n";
    genomePar << "sjdbGTFtagExonParentGene\t" << mapGen.pGe.sjdbGTFtagExonParentGene <<"\n";

    genomePar << "sjdbInsertSave\t" << mapGen.pGe.sjdbInsertSave <<"\n";
    
    genomePar << "genomeFileSizes\t" << mapGen.pGe.gFileSizes.at(0);
    for (uint ii=1;ii<mapGen.pGe.gFileSizes.size();ii++) 
        genomePar << " " << mapGen.pGe.gFileSizes.at(ii) ;
    genomePar << "\n";

    genomePar.close();
};
