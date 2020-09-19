/*
    functions that handle errors and warnings
*/
#include "ErrorWarning.h"
#include "TimeFunctions.h"
#include "GlobalVariables.h"

void exitWithError(string messageOut, ostream &streamOut1, ostream &streamOut2, int errorInt, Parameters &P)
{
    if (P.runThreadN>1) 
        pthread_mutex_lock(&g_threadChunks.mutexError);
    time_t timeCurrent;
    time( &timeCurrent);
    if (streamOut1.good()) {
        streamOut1 << "\n" << messageOut << endl << timeMonthDayTime(timeCurrent) <<" ...... FATAL ERROR, exiting\n"  <<flush;
    };
    if (streamOut2.good()) {
        streamOut2 << "\n" << messageOut << endl << timeMonthDayTime(timeCurrent) <<" ...... FATAL ERROR, exiting\n"  <<flush;
    };
    delete P.inOut; //to close files
//     if (P.runThreadN>1) pthread_mutex_unlock(&g_threadChunks.mutexError);
    exit(errorInt);
};

void warningMessage(string messageOut, ostream &streamOut1, ostream &streamOut2, Parameters &P)
{
    if (P.runThreadN>1) 
        pthread_mutex_lock(&g_threadChunks.mutexError);
    if (streamOut1.good()) {
        streamOut1 << "!!!!! WARNING: " << messageOut <<endl;
    };
    if (streamOut2.good()) {
        streamOut2 << "!!!!! WARNING: " << messageOut <<endl;
    };
    if (P.runThreadN>1) 
        pthread_mutex_unlock(&g_threadChunks.mutexError);    
};
