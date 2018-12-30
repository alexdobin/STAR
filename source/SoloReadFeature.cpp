#include "SoloReadFeature.h"
#include "streamFuns.h"

SoloReadFeature::SoloReadFeature(int32 feTy, Parameters &Pin, int iChunk) 
             : featureType(feTy), P(Pin), pSolo(P.pSolo)
{
    if (pSolo.type==0)
        return;
    
    for (uint32 ii=0; ii<stats.nStats; ii++)
        stats.V[ii]=0;
    
    cbReadCount = new uint32[pSolo.cbWL.size()];
    for (uint32 ii=0; ii<pSolo.cbWL.size(); ii++) {
        cbReadCount[ii]=0;
    };
    
    if (iChunk>=0) {
        strU_0 = &fstrOpen(P.outFileTmp+"/solo"+pSolo.featureNames[featureType]+"_0_"+std::to_string(iChunk),ERROR_OUT, P);
        strU_1 = &fstrOpen(P.outFileTmp+"/solo"+pSolo.featureNames[featureType]+"_1_"+std::to_string(iChunk),ERROR_OUT, P);
        strU_2 = &fstrOpen(P.outFileTmp+"/solo"+pSolo.featureNames[featureType]+"_2_"+std::to_string(iChunk),ERROR_OUT, P);        
    };
};

void SoloReadFeature::addCounts(const SoloReadFeature &rfIn)
{   
    for (uint32 ii=0; ii<pSolo.cbWL.size(); ii++) {
        cbReadCount[ii] += rfIn.cbReadCount[ii];
    };
};

void SoloReadFeature::addStats(const SoloReadFeature &rfIn)
{
    for (uint32 ii=0; ii<stats.nStats; ii++)
        stats.V[ii] += rfIn.stats.V[ii];
};

void SoloReadFeature::statsOut(ofstream &streamOut)
{
    //streamOut << setw(50) << "CELL BARCODES IN READS:\n"
    for (uint32 ii=0; ii<stats.nStats; ii++) {
        streamOut << setw(50) << stats.names.at(ii) << setw(15) << stats.V[ii] << '\n';
    };
    streamOut.flush();
};
