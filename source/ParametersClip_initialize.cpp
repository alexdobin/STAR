#include "ParametersClip.h"
#include "Parameters.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"

void ParametersClip::initialize(Parameters *pPin, vector<vector<ClipMate>> &clipMates)
{
    pP = pPin;
    yes = false;
    
    if (adapterType[0]!="Hamming" && adapterType[0]!="CellRanger4") {
        exitWithError("EXITING because of fatal PARAMETER error: --adapterType = " + adapterType[0] + " is not a valid option\n" +
                      "SOLUTION: use valid --adapterType options: Hamming OR CellRanger4\n", std::cerr, pPin->inOut->logMain, EXIT_CODE_PARAMETER, *pPin);
    };
    
    clipMates.resize(min(2LLU, pP->readNmates));//3rd read can only be barcode - do not trim it (for now).

	for (uint32 im=0; im<clipMates.size(); im++) {
		for (int ip=0; ip<2; ip++) {//fill in the ip
			clipMates[im].resize(2);
			clipMates[im][ip].type=ip;
            
            if (adapterType[0]=="CellRanger4") {
                clipMates[im][ip].type += 10;
                if (ip==0 && in[ip].adSeq[im]=="-") {
                    in[ip].adSeq[im] = "AAGCAGTGGTATCAACGCAGAGTACATGGG"; //standard TSO adapter
                } else if (ip==1) {
                    in[ip].adSeq[im] = "A"; //always polyA
                };
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
