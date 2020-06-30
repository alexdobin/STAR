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
     
    //prepare for feature-specific counting:
    if (featureType==SoloFeatureTypes::SJ && P.sjAll[0].empty()) {
        ifstream &sjIn = ifstrOpen(P.outFileTmp+"SJ.start_gap.tsv",  ERROR_OUT, "SOLUTION: re-run STAR", P);
        P.sjAll[0].reserve(10000000);
        P.sjAll[1].reserve(10000000);
        uint64 start1, gap1;
        while ( sjIn >> start1 >> gap1 ) {            
            P.sjAll[0].emplace_back(start1);
            P.sjAll[1].emplace_back(gap1);
        };
        sjIn.close();
        P.inOut->logMain <<"Read splice junctions for Solo SJ feature: "<< P.sjAll[0].size() <<endl;
    };
    
    //number of features
    switch (featureType) {
    	case SoloFeatureTypes::Gene :
    	case SoloFeatureTypes::GeneFull :
    	case SoloFeatureTypes::Velocyto :
    		featuresNumber=Trans.nGe;
    		break;
    	case SoloFeatureTypes::SJ :
    		featuresNumber=P.sjAll[0].size();
	};

    sumThreads(RAchunk);
    
    
    //call counting method
    if (featureType==SoloFeatureTypes::Velocyto) {
        countVelocyto();
    } else if (featureType==SoloFeatureTypes::Transcript3p) {
        quantTranscript();
        return;
    } else {//all others, standard processing
        if (pSolo.type==pSolo.SoloTypes::SmartSeq) {
            countSmartSeq();
        } else {
            countCBgeneUMI();
        };
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
