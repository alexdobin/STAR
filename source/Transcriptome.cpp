#include "Transcriptome.h"
#include "ErrorWarning.h"

Transcriptome::Transcriptome (Parameters* Pin) {
    
    P=Pin;

    //load tr and ex info
    ifstream trinfo((P->genomeDir+"/transcriptInfo.tab").c_str());
    if (trinfo.fail()) {//could not open file
        ostringstream errOut;
        errOut << "ERROR_01101: exiting because of *INPUT FILE* error: could not open for reading "<< (P->genomeDir+"/transcriptInfo.tab") <<"\n";
        errOut << "SOLUTION: if this file is missing from the genome directory, you will need to *re-generate the genome*, \n";
        errOut << "          if this file is present, check its read permissions\n";
        exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_FILE_OPEN, *P);
    };

    trinfo >> nTr;
    trS=new uint [nTr];
    trE=new uint [nTr];
    trEmax=new uint [nTr];
    trExI=new uint32 [nTr];
    trExN=new uint16 [nTr];
    trStr=new uint8 [nTr];
    trID.resize(nTr);
    
    for (uint32 itr=0; itr<nTr; itr++) {
        uint16 str1;
        trinfo >> trID[itr] >> trS[itr] >> trE[itr] >> trEmax[itr] >> str1 >> trExN[itr] >> trExI[itr];
        trStr[itr]=str1;
    };
    P->inOut->logMain << "Loaded transcript database, nTr="<<nTr<<endl;

    trinfo.close();
    
    ifstream exinfo((P->genomeDir+"/exonInfo.tab").c_str());

    exinfo >> nEx;
    exSE = new uint32 [2*nEx];
    exLenCum = new uint32 [nEx];
    for (uint32 iex=0; iex<nEx; iex++) {
        exinfo >> exSE[2*iex] >> exSE[2*iex+1] >> exLenCum[iex]; //reading all elements one after another
    };
    P->inOut->logMain << "Loaded exon database, nEx="<<nEx<<endl;

    exinfo.close();
    
    //
    Ptr=Pin;

};