#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"

void SoloFeature::processRecords(ReadAlignChunk **RAchunk)
{
    if (pSolo.type==0)
        return;

    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Starting Solo post-map for " <<SoloFeatureTypes::Names[featureType] <<endl;
  
    ///////////////////////////// collect RAchunk->RA->soloRead->readFeat
    nReadsInput=0;
    for (int ii=0; ii<P.runThreadN; ii++) {//point to
        readFeatAll[ii]= RAchunk[ii]->RA->soloRead->readFeat[pSolo.featureInd[featureType]];
        nReadsInput = max(nReadsInput,RAchunk[ii]->RA->iReadAll);
    };
    ++nReadsInput;
    
    countCBgeneUMI();

    ofstream *statsStream = &ofstrOpen(outputPrefix+"Features.stats",ERROR_OUT, P);
    readFeatSum->statsOut(*statsStream);
    statsStream->close();
    
    //output nU per gene per CB
    outputResults(false); //unfiltered
    
    if (pSolo.cellFilter.type[0]!="None" && (featureType==SoloFeatureTypes::Gene || featureType==SoloFeatureTypes::GeneFull)) {
        cellFiltering();
        outputResults(true);
    };
    
    //summary stats output
    statsOutput();
    
    //delete big arrays allocated in the previous functions
    delete[] rGeneUMI;
    delete[] rCBp;
    delete[] indCB;

};
