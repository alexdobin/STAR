#include "SoloFeature.h"
#include "streamFuns.h"
#include "ErrorWarning.h"

SoloFeature::SoloFeature(int feTy, Parameters &Pin, Transcriptome &inTrans)
            :featureType(feTy), P(Pin), Trans(inTrans), pSolo(P.pSolo)
{

    readFeatSum = new SoloReadFeature(featureType,P,-1);
    readBarSum = new SoloReadBarcode(P);
    readFeatAll = new SoloReadFeature*[P.runThreadN];
    readBarAll = new SoloReadBarcode*[P.runThreadN];

    if (pSolo.type==0)
        return;

    outputPrefix=P.outFileNamePrefix+pSolo.outFileNames[0];
    if (featureType!=0) {//
        outputPrefix += pSolo.featureNames[featureType] +'/';
        if (outputPrefix.back()=='/' && mkdir(outputPrefix.c_str(),P.runDirPerm)!=0 && errno!=EEXIST) {//create directory
            ostringstream errOut;
            errOut << "EXITING because of fatal OUTPUT FILE error: could not create Solo output directory"<<outputPrefix<<"\n";
            errOut << "SOLUTION: check the path and permisssions";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        };
    };
    
    if (featureType==3) {//for now - open output file
        streamTranscriptsOut = &ofstrOpen(outputPrefix+"transcriptCbUmiIndexDistance.txt",ERROR_OUT, P);
    };
};
