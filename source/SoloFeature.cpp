#include "SoloFeature.h"
#include "streamFuns.h"
#include "ErrorWarning.h"

SoloFeature::SoloFeature(int32 feTy, Parameters &Pin, Transcriptome &inTrans, SoloReadBarcode *readBarSumIn, SoloFeature **soloFeatAll)
            :featureType(feTy), P(Pin), Trans(inTrans), soloFeatAll(soloFeatAll), pSolo(P.pSolo), readBarSum(readBarSumIn)
{

    readFeatSum = new SoloReadFeature(featureType,P,-1);
    readFeatAll = new SoloReadFeature*[P.runThreadN];

    if (pSolo.type==0)
        return;

    outputPrefix=P.outFileNamePrefix+pSolo.outFileNames[0];
    outputPrefix += SoloFeatureTypes::Names[featureType] +'/';
    if (mkdir(outputPrefix.c_str(),P.runDirPerm)!=0 && errno!=EEXIST) {//create directory
        ostringstream errOut;
        errOut << "EXITING because of fatal OUTPUT FILE error: could not create Solo output directory"<<outputPrefix<<"\n";
        errOut << "SOLUTION: check the path and permisssions";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    };
    
    if (featureType==SoloFeatureTypes::Transcript3p) {//for now - open output file
        streamTranscriptsOut = &ofstrOpen(outputPrefix+"transcriptCbUmiIndexDistance.txt",ERROR_OUT, P);
    };
};
