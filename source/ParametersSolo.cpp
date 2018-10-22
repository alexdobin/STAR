#include "ParametersSolo.h"
#include "Parameters.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include "SequenceFuns.h"
#include "serviceFuns.cpp"

#include <stdlib.h>

void ParametersSolo::initialize(Parameters *pPin) {
    pP=pPin;
    
    if (typeStr=="None") {
        type=0;
    } else if (typeStr=="10XchromiumV2") {
        type=1;
        bL=cbL+umiL;
        pP->readNmates=1; //output mates TODO: check that readNmatesIn==2
    } else  {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in of --soloType="<<typeStr<<"\n";
        errOut << "SOLUTION: use allowed option: None or 10XchromiumV2";
        exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_PARAMETER, *pP);
    }; 
    
    //load the CB whitlist and create unordered_map
    ifstream & cbWlStream = ifstrOpen(soloCBwhitelist, ERROR_OUT, "SOLUTION: check the path and permissions of the CB whitelist file: " + soloCBwhitelist, *pP);
    string seq1;
    while (!cbWlStream.eof()) {
        cbWlStream >> seq1;
        if (seq1.size() != cbL) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in input CB whitelist file: "<< soloCBwhitelist <<" the total length of barcode sequence is "  << seq1.size() << " not equal to expected " <<bL <<"\n"  ;
            errOut << "SOLUTION: make sure that the barcode read is the second in --readFilesIn and check that is has the correct formatting\n";
            exitWithError(errOut.str(),std::cerr, pP->inOut->logMain, EXIT_CODE_INPUT_FILES, *pP);
        };
        uint32 cb1;
        //convert to 2-bit format
        if (convertNuclStrToInt32(seq1,cb1)) {
            //cbWL.insert(cb1);
            cbWL.push_back(cb1);
        } else {
            pP->inOut->logMain << "WARNING: CB whitelist sequence contains non-ACGT and is ignored: " << seq1 <<endl;
        };
    };
    //int comp1 = [] (const void *a, const void *b) {uint32 aa=*(uint32*) a; uint32 bb=*(uint32*) b; if (a 
    qsort(cbWL.data(),cbWL.size(),sizeof(uint32),funCompareNumbers<uint32>);
    
    time_t rawTime;
    time(&rawTime);
    pP->inOut->logMain << timeMonthDayTime(rawTime) << "Finished reading CB whitelist sequences: " << cbWL.size() <<endl;
    *pP->inOut->logStdOut << timeMonthDayTime(rawTime) << " Finished reading CB whitelist sequences: " << cbWL.size() <<endl;
    
    uint64 nn=0;
    for (uint ii=0; ii<100000000; ii++){
        uint32 a=ii;//(uint32) rand(); 
//         nn+=cbWL.count(a);
    };
    time(&rawTime);
    *pP->inOut->logStdOut << timeMonthDayTime(rawTime) <<" "<< nn <<endl;
    exit(1);
};