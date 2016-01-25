/*
functions that handle errors and warnings
*/
#include "ErrorWarning.h"
#include "TimeFunctions.h"

void exitWithError(string messageOut, ostream &streamOut1, ostream &streamOut2, int errorInt, Parameters &P) {
    time_t timeCurrent;
    time( &timeCurrent);
    if (streamOut1.good()) {
        streamOut1 << "\n" << messageOut << endl << timeMonthDayTime(timeCurrent) <<" ...... FATAL ERROR, exiting\n"  <<flush;
    };
    if (streamOut2.good()) {
        streamOut2 << "\n" << messageOut << endl << timeMonthDayTime(timeCurrent) <<" ...... FATAL ERROR, exiting\n"  <<flush;
    };
    delete P.inOut; //to close files

    exit(errorInt);
};
