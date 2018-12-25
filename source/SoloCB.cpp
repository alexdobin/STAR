#include "SoloCB.h"
#include "streamFuns.h"

SoloCB::SoloCB(int32 feTy, Parameters &Pin, int iChunk) 
             : featureType(feTy), P(Pin), pSolo(P.pSolo)
{
    if (pSolo.type==0 || !pSolo.featureYes[featureType])
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
//         vector<string> prefix={"soloGene","soloSJ"};
//         constexpr string prefix1=prefix[soloFeature];
        string prefix1="solo";
//         if (soloFeature==1) {
//             prefix=;
//         } else if (soloFeature==1) {
//             prefix=;
//         };
        strU_0 = &fstrOpen(P.outFileTmp+"/"+prefix1+"_U_0_"+std::to_string(iChunk),ERROR_OUT, P);
        strU_1 = &fstrOpen(P.outFileTmp+"/"+prefix1+"_U_1_"+std::to_string(iChunk),ERROR_OUT, P);
        strU_2 = &fstrOpen(P.outFileTmp+"/"+prefix1+"_U_2_"+std::to_string(iChunk),ERROR_OUT, P);        
    };
    
    for (uint32 jj=0;jj<4;jj++) {
        homoPolymer[jj]=0;
        for (uint32 ii=0; ii<pSolo.umiL;ii++) {
            homoPolymer[jj]=(homoPolymer[jj]<<2)+jj;
        };
    };
};

void SoloCB::addSoloCBcounts(const SoloCB &soloCBin)
{   
    for (uint32 ii=0; ii<pSolo.cbWL.size(); ii++) {
        cbReadCount[ii] += soloCBin.cbReadCount[ii];
        cbReadCountExact[ii] += soloCBin.cbReadCountExact[ii];
    };
};

void SoloCB::addSoloCBstats(const SoloCB &soloCBin)
{
    for (uint32 ii=0; ii<stats.nStats; ii++)
        stats.V[ii] += soloCBin.stats.V[ii];
};

void SoloCB::statsOut(ofstream &streamOut)
{
    //streamOut << setw(50) << "CELL BARCODES IN READS:\n"
    for (uint32 ii=0; ii<stats.nStats; ii++) {
        streamOut << setw(50) << stats.names.at(ii) << setw(15) << stats.V[ii] << '\n';
    };
    streamOut.flush();
};
