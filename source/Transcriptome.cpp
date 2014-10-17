#include "Transcriptome.h"

Transcriptome::Transcriptome (Parameters* Pin) {
    
    P=Pin;

    //load tr ans ex info
    ifstream trinfo((P->genomeDir+"/transcriptInfo.tab").c_str());

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