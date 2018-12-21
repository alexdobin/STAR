#include "SoloCB.h"
#include "streamFuns.h"

SoloCB::SoloCB(Parameters &Pin, int iChunk) : P(Pin), pSolo(P.pSolo) {
    if (pSolo.type==0)
        return;
    
    for (uint32 ii=0; ii<stats.nStats; ii++)
        stats.V[ii]=0;
    
    cbReadCount = new uint32[pSolo.cbWL.size()];
    cbReadCountExact = new uint32[pSolo.cbWL.size()];    
    for (uint32 ii=0; ii<pSolo.cbWL.size(); ii++) {
        cbReadCount[ii]=0;
        cbReadCountExact[ii]=0;        
    };
    
    if (iChunk>=0) {
        strU_0 = &fstrOpen(P.outFileTmp+"/soloCB_U_0_"+std::to_string(iChunk),ERROR_OUT, P);
        strU_1 = &fstrOpen(P.outFileTmp+"/soloCB_U_1_"+std::to_string(iChunk),ERROR_OUT, P);
        strU_2 = &fstrOpen(P.outFileTmp+"/soloCB_U_2_"+std::to_string(iChunk),ERROR_OUT, P);        
    };
    
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
    
    for (uint32 ii=0; ii<pSolo.cbWL.size(); ii++) {
        cbReadCount[ii] += soloCBin.cbReadCount[ii];
        cbReadCountExact[ii] += soloCBin.cbReadCountExact[ii];
    };
};

void SoloCB::statsOut(ofstream &streamOut) {
    for (uint32 ii=0; ii<stats.nStats; ii++) {
        streamOut << setw(25) << stats.names.at(ii) << setw(15) << stats.V[ii] << '\n';
    };
    streamOut.flush();
};

void SoloCB::readCBgeneUMIfromFiles(uint32 ** cbP, uint32 *cbReadCountExactTotal) {
    
    {//load exact matches
        strU_0->flush();
        strU_0->seekg(0,ios::beg);
        uint32 cb1, g1, umi1;
        while (*strU_0 >> cb1 >> g1 >> umi1) {
            //strU_0->ignore(1000000000,'\n'); //in case more fields were output
            cbP[cb1][0]=g1;
            cbP[cb1][1]=umi1;
            cbP[cb1]+=2;
        };
    };
    
    {//load 1MM
        strU_1->flush();
        strU_1->seekg(0,ios::beg);
        uint32 cb1, g1, umi1;
        while (*strU_1 >> cb1 >> g1 >> umi1) {
            if (cbReadCountExactTotal[cb1]>0) {
                cbP[cb1][0]=g1;
                cbP[cb1][1]=umi1;
                cbP[cb1]+=2;
            };
        };
    };
    
    {//load 2MM
        strU_2->flush();
        strU_2->seekg(0,ios::beg);
        uint32 ncb1, cb1, g1, umi1;
        while (*strU_2 >> g1 >> umi1 >> ncb1) {
            float ptot=0.0,pmax=0.0;
            for (uint32 ii=0; ii<ncb1; ii++) {
                uint32 cbin;
                char  qin;
                float pin;
                *strU_2 >> cbin >> qin;
                if (cbReadCountExactTotal[cbin]>0) {//otherwise this cbin does not work
                    qin -= pSolo.QSbase;
                    qin = qin < pSolo.QSbase ? qin : pSolo.QSbase;
                    pin=cbReadCountExactTotal[cbin]*pow(10.0,-qin/10.0);
                    ptot+=pin;
                    if (pin>pmax) {
                        cb1=cbin;
                        pmax=pin;
                    };
                };
            };
            if (ptot>0.0 && pmax>=pSolo.cbMinP*ptot) {
                cbP[cb1][0]=g1;
                cbP[cb1][1]=umi1;
                cbP[cb1]+=2;
            };        
        };
    };
};
