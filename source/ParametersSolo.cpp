#include "ParametersSolo.h"

ParametersSolo::initialize(Parameters *pPin) : pP (pPin) {
    if (typeStr=="None") {
        type=0;
    } else if (typeStr=="10XchromiumV2") {
        type=1;
        bL=cbL+umiL;
        readNmates=1; //output mates TODO: check that readNmatesIn==2
    } else  {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETERS error: unrecognized option in of --soloType="<<typeStr<<"\n";
        errOut << "SOLUTION: use allowed option: None or 10XchromiumV2";
        exitWithError(errOut.str(),std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };    
}