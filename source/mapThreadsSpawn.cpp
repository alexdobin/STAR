#include "mapThreadsSpawn.h"
#include "ThreadControl.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"

void mapThreadsSpawn (Parameters &P, ReadAlignChunk** RAchunk) {
    for (int ithread=1;ithread<P.runThreadN;ithread++) {//spawn threads
        int threadStatus=pthread_create(&g_threadChunks.threadArray[ithread], NULL, &g_threadChunks.threadRAprocessChunks, (void *) RAchunk[ithread]);
        if (threadStatus>0) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while creating thread # " << ithread <<", error code: "<<threadStatus ;
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, 1, P);
        };
        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
        P.inOut->logMain << "Created thread # " <<ithread <<"\n"<<flush;
        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    };

    RAchunk[0]->processChunks(); //start main thread

    for (int ithread=1;ithread<P.runThreadN;ithread++) {//wait for all threads to complete
        int threadStatus = pthread_join(g_threadChunks.threadArray[ithread], NULL);
        if (threadStatus>0) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while joining thread # " << ithread <<", error code: "<<threadStatus ;
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, 1, P);
        };
        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
        P.inOut->logMain << "Joined thread # " <<ithread <<"\n"<<flush;
        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    };
};

