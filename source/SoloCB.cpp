#include "SoloCB.h"
#include "streamFuns.h"

SoloCB::SoloCB(Parameters &Pin, int iChunk) : P(Pin), pSolo(P.pSolo) {
    if (pSolo.type==0)
        return;
    
    for (uint32 ii=0; ii<stats.nStats; ii++)
        stats.V[ii]=0;
    
    cbReadCount = new uint32[pSolo.cbWL.size()];
    for (uint32 ii=0; ii<pSolo.cbWL.size(); ii++)
        cbReadCount[ii]=0;
    
    if (iChunk>=0)
        strU_0 = &fstrOpen(P.outFileTmp+"/soloCB_U_0_"+std::to_string(iChunk),ERROR_OUT, P);
    //strU_1 = &ofstrOpen(P.outFileTmp+"/soloCB_U_1_"+std::to_string(iChunk),ERROR_OUT, P);
    
    for (uint32 jj=0;jj<4;jj++) {
        homoPolymer[jj]=0;
        for (uint32 ii=0; ii<pSolo.umiL;ii++) {
            homoPolymer[jj]=(homoPolymer[jj]<<2)+jj;
        };
    };
};

void SoloCB::addSoloCB(const SoloCB &soloCBin) {
    for (uint32 ii=0; ii<stats.nStats; ii++)
        stats.V[ii] += soloCBin.stats.V[ii];
    
    for (uint32 ii=0; ii<pSolo.cbWL.size(); ii++)
        cbReadCount[ii] += soloCBin.cbReadCount[ii];
};

void SoloCB::statsOut(ofstream &streamOut) {
    for (uint32 ii=0; ii<stats.nStats; ii++) {
        streamOut << setw(25) << stats.names.at(ii) << setw(15) << stats.V[ii] << '\n';
    };
    streamOut.flush();
};

void SoloCB::readCBgeneUMIfromFiles(uint32 ** cbP) {
    strU_0->flush();
    strU_0->seekg(0,ios::beg);
     
    uint32 cb1, g1, umi1;
    for (uint64 ii=0; ii<stats.V[stats.nMatch]; ii++) {
        *strU_0 >> cb1 >> g1 >> umi1;
        //strU_0->ignore(1000000000,'\n'); //in case more fields were output, may remove this
        cbP[cb1][0]=g1;
        cbP[cb1][1]=umi1;
        cbP[cb1]+=2;
    };
};
