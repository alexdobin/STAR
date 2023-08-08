#include "ParametersGenome.h"
#include "Parameters.h"
#include "ErrorWarning.h"

void ParametersGenome::initialize(Parameters *pPin)
{
    pP=pPin;
    
    if (gDir.back()!='/') {
        gDir += '/';
    };
    
    //genome transformation
    if (transform.typeString=="None") {
        transform.type=0;
    } else if (transform.typeString=="Haploid") {
        transform.type=1;
    } else if (transform.typeString=="Diploid") {
        transform.type=2;
    } else {
        ostringstream errOut;
        errOut << "EXITING because of FATAL PARAMETER ERROR: unrecognized option in --outTransformType = " << transform.typeString << "\n";
        errOut << "SOLUTION: use one of the allowed values for --outTransformType : 'None' or 'Haploid' or 'Diploid' \n";
        exitWithError(errOut.str(), std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };
    
    transform.outYes = transform.outSAM = transform.outSJ = false;
    if (transform.output.at(0) == "None") {
        //nothing to do
    } else {
        for (auto &ot: transform.output) {
            if (ot == "SAM") {
                transform.outYes = transform.outSAM = true;
            } else if (ot == "SJ") {
                transform.outYes = transform.outSJ = true;
            } else if (ot == "Quant") {
                transform.outYes = transform.outQuant = true;                
            } else {
                exitWithError("EXITING because of FATAL PARAMETER ERROR: unrecognized option in --outTransformOutput = " + ot + '\n'
                              + "SOLUTION: use allowed values for --outTransformOutput: None or SAM and/or SJ\n"
                              ,std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
            };
        };
    };
    
    if (gTypeString!="Full" && gTypeString!="Transcriptome" && gTypeString!="SuperTranscriptome") {
        ostringstream errOut;
        errOut << "EXITING because of FATAL parameter error: --genomeType=" << gTypeString << "\n";
        errOut << "SOLUTION: use one of the allowed values for --genomeLoad : Full OR Transcriptome OR SuperTranscriptome\n" <<flush;
        exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };

    if (gLoad!="LoadAndKeep" && gLoad!="LoadAndRemove" && gLoad!="Remove" && gLoad!="LoadAndExit" && gLoad!="NoSharedMemory") {// find shared memory fragment
        ostringstream errOut;
        errOut << "EXITING because of FATAL INPUT ERROR: --genomeLoad=" << gLoad << "\n" <<flush;
        errOut << "SOLUTION: use one of the allowed values for --genomeLoad : NoSharedMemory,LoadAndKeep,LoadAndRemove,LoadAndExit,Remove.\n" <<flush;
        exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    };    
};