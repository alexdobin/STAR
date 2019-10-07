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
     
    sumThreads(RAchunk);
    
    if (featureType==SoloFeatureTypes::Velocyto) {
        countVelocyto();
    } else {//all others, standard processing
        countCBgeneUMI();
    };

    //output
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
    //delete[] indCB;
};
