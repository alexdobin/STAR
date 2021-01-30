#include "SoloReadBarcode.h"
#include "streamFuns.h"

SoloReadBarcode::SoloReadBarcode(Parameters &P) : P(P), pSolo(P.pSolo)
{
    if (pSolo.type==0)
        return;

    for (uint32 ii=0; ii<stats.nStats; ii++)
        stats.V[ii]=0;

    if (pSolo.cbWLyes) {
        cbReadCountExact.resize(pSolo.cbWLsize,0);
    };

    for (uint32 jj=0;jj<4;jj++) {
        homoPolymer[jj]=0;
        for (uint32 ii=0; ii<pSolo.umiL;ii++) {
            homoPolymer[jj]=(homoPolymer[jj]<<2)+jj;
        };
    };
    
    qualHist.fill(0);
};

void SoloReadBarcode::addCounts(const SoloReadBarcode &rfIn)
{
    if (pSolo.cbWLyes) {
        for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
            cbReadCountExact[ii] += rfIn.cbReadCountExact[ii];
        };
    };
        
    for (uint32 ii=0; ii<qualHist.size(); ii++)
        qualHist[ii] += rfIn.qualHist[ii];
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
