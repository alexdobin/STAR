#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"
#include "systemFunctions.h"

void SoloFeature::processRecords()
{
    if (pSolo.type==0)
        return;

    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Starting Solo post-map for " <<SoloFeatureTypes::Names[featureType] <<endl;
    
    outputPrefix= P.outFileNamePrefix+pSolo.outFileNames[0];
    outputPrefix+= SoloFeatureTypes::Names[featureType] +'/';
    outputPrefixFiltered= outputPrefix + "filtered/";
    
    if (mkdir(outputPrefix.c_str(),P.runDirPerm)!=0 && errno!=EEXIST) {//create directory
        ostringstream errOut;
        errOut << "EXITING because of fatal OUTPUT FILE error: could not create Solo output directory"<<outputPrefix<<"\n";
        errOut << "SOLUTION: check the path and permisssions";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    };    
     
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

    SoloFeature::sumThreads();
    
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
    
    // no need to include multi-gene in yesWLmatch:
    // for multi-gene options, it's already done in _collapseUMI_ function
    //if (!(pSolo.multiMap.yes.multi && (featureType==SoloFeatureTypes::Gene || featureType==SoloFeatureTypes::GeneFull)))
    //    readFeatSum->stats.V[readFeatSum->stats.yesWLmatch] += readFeatSum->stats.V[readFeatSum->stats.MultiFeature];   

    //output
    ofstream *statsStream = &ofstrOpen(outputPrefix+"Features.stats",ERROR_OUT, P);
    readFeatSum->statsOut(*statsStream);
    statsStream->close();
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Solo: writing raw matrix" <<endl;

    //output nU per gene per CB
    outputResults(false,  outputPrefix + "/raw/"); //unfiltered
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Solo: cell filtering" <<endl;    
    cellFiltering();
    
    //summary stats output
    statsOutput();
    
    //delete big arrays allocated in the previous functions
    clearLarge();
    //delete[] indCB;

    P.inOut->logMain << "RAM after completing solo:\n"
                     <<  linuxProcMemory() << flush;   
};
