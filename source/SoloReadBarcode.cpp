#include "SoloReadBarcode.h"
#include "streamFuns.h"

SoloReadBarcode::SoloReadBarcode(Parameters &Pin) : P(Pin), pSolo(P.pSolo)
{
    if (pSolo.type==0)
        return;
    
    for (uint32 ii=0; ii<stats.nStats; ii++)
        stats.V[ii]=0;    
    
    cbReadCountExact = new uint32[pSolo.cbWL.size()];    
    for (uint32 ii=0; ii<pSolo.cbWL.size(); ii++) {
        cbReadCountExact[ii]=0;        
    };
    
    for (uint32 jj=0;jj<4;jj++) {
        homoPolymer[jj]=0;
        for (uint32 ii=0; ii<pSolo.umiL;ii++) {
            homoPolymer[jj]=(homoPolymer[jj]<<2)+jj;
        };
    };
};

void SoloReadBarcode::addCounts(const SoloReadBarcode &rfIn)
{   
    for (uint32 ii=0; ii<pSolo.cbWL.size(); ii++) {
        cbReadCountExact[ii] += rfIn.cbReadCountExact[ii];
    };
};

void SoloReadBarcode::addStats(const SoloReadBarcode &rfIn)
{
    for (uint32 ii=0; ii<stats.nStats; ii++)
        stats.V[ii] += rfIn.stats.V[ii];
};

void SoloReadBarcode::statsOut(ofstream &streamOut)
{
    //streamOut << setw(50) << "CELL BARCODES IN READS:\n"
    for (uint32 ii=0; ii<stats.nStats; ii++) {
        streamOut << setw(50) << stats.names.at(ii) << setw(15) << stats.V[ii] << '\n';
    };
    streamOut.flush();
};
