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
                      "SOLUTION: use valid --clipAdapterType options: Hamming OR CellRanger4\n", std::cerr, pPin->inOut->logMain, EXIT_CODE_PARAMETER, *pPin);
    };
    
    if (adapterType[0]=="CellRanger4") {
        
        if (in[1].adSeq.size()>1 || in[1].adSeq[0]!="-") {
            exitWithError("EXITING because of fatal PARAMETER error: --clipAdapterType CellRanger4 uses fixed sequences for 3' polyA adapters."
                          "\nSOLUTION: Do not use --clip5pAdapter* options.\n", std::cerr, pPin->inOut->logMain, EXIT_CODE_PARAMETER, *pPin);
        };
        in[1].adSeq[0] = "A";
        
        if (in[0].adSeq.size()>1 && in[0].adSeq[1]!="-") {
            exitWithError("EXITING because of fatal PARAMETER error: when using --clipAdapterType CellRanger4, only 5' adapter for the 1st mate can be specified."
                          "\nSOLUTION: Use only one sequence in --clip5pAdapterSeq (or '-' instead of 2nd sequence).\n", std::cerr, pPin->inOut->logMain, EXIT_CODE_PARAMETER, *pPin);
        };

        if (in[0].adSeq[0]=="-") {//standard 10X TSO
            in[0].adSeq[0] = "AAGCAGTGGTATCAACGCAGAGTACATGGG";
        };

    } else {
        for (uint32 im=0; im<in[0].adSeq.size(); im++) {//TODO implement standard 5p clipping
            if (in[0].adSeq[im]!="-") {
                exitWithError("EXITING because of fatal PARAMETER error: --clip5pAdapterSeq is not supported yet, except for --clipAdapterType CellRanger4.\
                            \nSOLUTION: Do not use --clip5pAdapter* options without --clipAdapterType CellRanger4.\n", std::cerr, pPin->inOut->logMain, EXIT_CODE_PARAMETER, *pPin);
            };      
        };
    };
   
    for (int ip=0; ip<2; ip++) {//if no clipping, repeat for all mates
        if (in[ip].adSeq[0]=="-") {
            in[ip].adSeq.resize(pP->readNmates, "-");
            in[ip].adMMp.resize(pP->readNmates, 0);
        };
            
        if (in[ip].N[0]==0)
            in[ip].N.resize(pP->readNmates, 0);
            
        if (in[ip].NafterAd[0]==0)
            in[ip].NafterAd.resize(pP->readNmates, 0);        
    };
    
    array<string,2> p53={"5","3"};

    for (int ip=0; ip<2; ip++) {//check that readNmates values are provided
        if (in[ip].adSeq.size() != pP->readNmates) {
            exitWithError("EXITING because of fatal PARAMETER error: --clip" +p53[ip]+ "pAdapterSeq has to contain " +to_string(pP->readNmates)+ " values to match the number of mates."
                            "\nSOLUTION: specify " +to_string(pP->readNmates)+ "values in --clip" +p53[ip]+ "pAdapterSeq , for no clipping use -"
                           , std::cerr, pPin->inOut->logMain, EXIT_CODE_PARAMETER, *pPin);
        };
            
        if (in[ip].adMMp.size() != pP->readNmates) {
            exitWithError("EXITING because of fatal PARAMETER error: --clip" +p53[ip]+ "pAdapterMMp has to contain " +to_string(pP->readNmates)+ " values to match the number of mates."
                            "\nSOLUTION: specify " +to_string(pP->readNmates)+ "values in --clip" +p53[ip]+ "pAdapterMMp"
                            , std::cerr, pPin->inOut->logMain, EXIT_CODE_PARAMETER, *pPin);
        };
            
        if (in[ip].NafterAd.size() != pP->readNmates) {
            exitWithError("EXITING because of fatal PARAMETER error: --clip" +p53[ip]+ "pAfterAdapterNbases has to contain " +to_string(pP->readNmates)+ " values to match the number of mates."
                            "\nSOLUTION: specify " +to_string(pP->readNmates)+ "values in --clip" +p53[ip]+ "pAfterAdapterNbases , for no clipping use 0"
                            , std::cerr, pPin->inOut->logMain, EXIT_CODE_PARAMETER, *pPin);
        };
        
        if (in[ip].N.size() != pP->readNmates) {
            exitWithError("EXITING because of fatal PARAMETER error: --clip" +p53[ip]+ "pNbases has to contain " +to_string(pP->readNmates)+ " values to match the number of mates."
                            "\nSOLUTION: specify " +to_string(pP->readNmates)+ "values in --clip" +p53[ip]+ "pNbases , for no clipping use 0"
                            , std::cerr, pPin->inOut->logMain, EXIT_CODE_PARAMETER, *pPin);
        };  
    };
};

void ParametersClip::initializeClipMates(vector<vector<ClipMate>> &clipMates)
{        
    clipMates.resize(pP->readNends);

    for (uint32 im=0; im<clipMates.size(); im++) {
        
        clipMates[im].resize(2);
        
        for (int ip=0; ip<2; ip++) {//fill in the ip
            clipMates[im][ip].type=ip;
            
            if (adapterType[0]=="CellRanger4") {
                clipMates[im][ip].type += 10;
            };            
            
            if (im<pP->readNmates) {//true mates
                clipMates[im][ip].initialize(in[ip].N[im], in[ip].adSeq[im], in[ip].NafterAd[im], in[ip].adMMp[im]);
            } else {//barcode read, no clipping
                //back() effectively duplicates the values for 2nd mate - if only one value was given in paramerter input
                clipMates[im][ip].initialize(0, "-", 0, 0); 
            };
        };
    };
};
