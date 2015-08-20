#include "genomeParametersWrite.h"
#include "streamFuns.h"

void genomeParametersWrite(string fileName, Parameters* P, string errorOut)
{//write the genome information into the genomePar stream
    ofstream & genomePar = ofstrOpen(fileName, errorOut, P);
    
    genomePar << "### "<<P->commandLineFull <<"\n";
    
    genomePar << "versionGenome\t" << P->versionSTAR <<"\n";
    genomePar << "genomeFastaFiles\t";
    for (uint ii=0;ii<P->genomeFastaFiles.size();ii++) genomePar << P->genomeFastaFiles.at(ii) << " ";
    genomePar << "\n";
    genomePar << "genomeSAindexNbases\t" << P->genomeSAindexNbases << "\n";
    genomePar << "genomeChrBinNbits\t" << P->genomeChrBinNbits << "\n";
    genomePar << "genomeSAsparseD\t" << P->genomeSAsparseD <<"\n";
    genomePar << "sjdbOverhang\t" << P->sjdbOverhang <<"\n";

    genomePar << "sjdbFileChrStartEnd\t";
    for (uint ii=0;ii<P->sjdbFileChrStartEnd.size();ii++) genomePar<< P->sjdbFileChrStartEnd.at(ii) << " ";
    genomePar<<"\n";
    
    genomePar << "sjdbGTFfile\t" << P->sjdbGTFfile <<"\n";
    genomePar << "sjdbGTFchrPrefix\t" << P->sjdbGTFchrPrefix <<"\n";
    genomePar << "sjdbGTFfeatureExon\t" << P->sjdbGTFfeatureExon <<"\n";
    genomePar << "sjdbGTFtagExonParentTranscript\t" << P->sjdbGTFtagExonParentTranscript <<"\n";
    genomePar << "sjdbGTFtagExonParentGene\t" << P->sjdbGTFtagExonParentGene <<"\n";
    
    genomePar << "sjdbInsertSave\t" << P->sjdbInsert.save <<"\n";
    genomePar.close();      
};   
