#include "SoloCB.h"
#include "streamFuns.h"

SoloBarcodes::SoloBarcodes(Parameters &Pin) 
             : featureType(feTy), P(Pin), pSolo(P.pSolo)
{
    if (pSolo.type==0)
        return;
    
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
