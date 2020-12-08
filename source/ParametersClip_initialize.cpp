#include "ParametersClip.h"
#include "Parameters.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"

void ParametersClip::initialize(Parameters *pPin)
{
    pP = pPin;
    //yes = false;
    
    if (adapterType[0]!="Hamming" && adapterType[0]!="CellRanger4") {
        exitWithError("EXITING because of fatal PARAMETER error: --clipAdapterType = " + adapterType[0] + " is not a valid option\n" +
                      "SOLUTION: use valid --adapterType options: Hamming OR CellRanger4\n", std::cerr, pPin->inOut->logMain, EXIT_CODE_PARAMETER, *pPin);
    };
       
	for (uint32 im=0; im<in[0].adSeq.size(); im++) {//TODO implement standard 5p clipping
        if (in[0].adSeq[im]!="-") {
            exitWithError("EXITING because of fatal PARAMETER error: --clip5pAdapterSeq is not supported yet.\
                           \nSOLUTION: Do not use --clip5pAdapter* options.\n", std::cerr, pPin->inOut->logMain, EXIT_CODE_PARAMETER, *pPin);
        };      
    };
    
    if (adapterType[0]=="CellRanger4") {
        
        if (in[1].adSeq.size()>1 || in[1].adSeq[0]!="-") {
            exitWithError("EXITING because of fatal PARAMETER error: --clipAdapterType CellRanger4 uses fixed sequences for 5' (TSO) and 3' (polyA) adapters.\
                           \nSOLUTION: Do not use --clip3pAdapter* or --clip5pAdapter* options.\n", std::cerr, pPin->inOut->logMain, EXIT_CODE_PARAMETER, *pPin);
        };
        
        in[0].adSeq[0] = "AAGCAGTGGTATCAACGCAGAGTACATGGG";
        in[1].adSeq[0] = "A";
    };

};

void ParametersClip::initializeClipMates(vector<vector<ClipMate>> &clipMates)
{        
    clipMates.resize(min(2LLU, pP->readNmates));//3rd read can only be barcode - do not trim it (for now).

    for (uint32 im=0; im<clipMates.size(); im++) {
        
        clipMates[im].resize(2);
        
        for (int ip=0; ip<2; ip++) {//fill in the ip
            clipMates[im][ip].type=ip;
            
            if (adapterType[0]=="CellRanger4") {
                clipMates[im][ip].type += 10;
            };            
            
            if (im==0) {
                clipMates[0][ip].initialize(in[ip].N[0], in[ip].adSeq[0], in[ip].NafterAd[0], in[ip].adMMp[0]);
            } else if (im==1) {
                //back() effectively duplicates the values for 2nd mate - if only one value was given in paramerter input
                clipMates[1][ip].initialize(in[ip].N.back(), in[ip].adSeq.back(), in[ip].NafterAd.back(), in[ip].adMMp.back()); 
            };
        };
    }; 
};
