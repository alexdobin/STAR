#ifndef CODE_ErrorWarning
#define CODE_ErrorWarning

#include "IncludeDefine.h"
#include "Parameters.h"

void exitWithError(string messageOut, ostream &streamOut1, ostream &streamOut2, int errorInt, const Parameters &P);
void warningMessage(string messageOut, ostream &streamOut1, ostream &streamOut2, Parameters &P);
#endif
