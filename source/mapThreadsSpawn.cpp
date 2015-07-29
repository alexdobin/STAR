#include <memory>
#include "mapThreadsSpawn.h"
#include "ThreadControl.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"


#if !defined(_WIN32) && defined(USE_PTHREAD)
#include <pthread.h>
#else
#include <mutex>
#include <thread>
#endif


#if !defined(_WIN32) && defined(USE_PTHREAD)

void mapThreadsSpawn (Parameters *P, const std::vector<ReadAlignChunk*> &RAchunk)) {
    for (int ithread=1;ithread<P->runThreadN;ithread++) {//spawn threads
        int threadStatus=pthread_create(&g_threadChunks.threadArray[ithread], NULL, &g_threadChunks.threadRAprocessChunks, (void *) RAchunk[ithread]);
        if (threadStatus>0) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while creating thread # " << ithread <<", error code: "<<threadStatus ;
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, 1, *P);
        };
        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
        P->inOut->logMain << "Created thread # " <<ithread <<"\n"<<flush;
        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    };
    
    RAchunk[0]->processChunks(); //start main thread
    
    for (int ithread=1;ithread<P->runThreadN;ithread++) {//wait for all threads to complete
        int threadStatus = pthread_join(g_threadChunks.threadArray[ithread], NULL);
        if (threadStatus>0) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while joining thread # " << ithread <<", error code: "<<threadStatus ;
                exitWithError(errOut.str(),std::cerr, P->inOut->logMain, 1, *P);
        };
        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
        P->inOut->logMain << "Joined thread # " <<ithread <<"\n"<<flush;        
        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    };
};

#else // ~ #if !defined(_WIN32) && defined(USE_PTHREAD)

void mapThreadsSpawn(Parameters *P, const std::vector<ReadAlignChunk*> &RAchunk) {
	for (int ithread = 1; ithread<P->runThreadN; ithread++) {//spawn threads
		try{
			std::thread tempThread(g_threadChunks.threadRAprocessChunks, RAchunk[ithread]);
			g_threadChunks.threadArray.push_back(std::move(tempThread));
		} 
		catch(std::exception &e){
			ostringstream errOut;
			errOut << "EXITING because of FATAL ERROR: phtread error while creating thread # " << e.what();
			exitWithError(errOut.str(), std::cerr, P->inOut->logMain, 1, *P);
		}
		g_threadChunks.mutexLogMain.lock();
		P->inOut->logMain << "Created thread # " << ithread << "\n" << flush;
		g_threadChunks.mutexLogMain.unlock();
	};

	RAchunk[0]->processChunks(); //start main thread

	for (int ithread = 1; ithread<P->runThreadN; ithread++) {//wait for all threads to complete
		try{
			g_threadChunks.threadArray[ithread].join();
		}
		catch (std::exception &e){
			ostringstream errOut;
			errOut << "EXITING because of FATAL ERROR: phtread error while joining thread # " << ithread << ", error code: " << e.what();
			exitWithError(errOut.str(), std::cerr, P->inOut->logMain, 1, *P);
		}
		g_threadChunks.mutexLogMain.lock();
		P->inOut->logMain << "Joined thread # " << ithread << "\n" << flush;
		g_threadChunks.mutexLogMain.unlock();
	};
};

#endif // ~ #if !defined(_WIN32) && defined(USE_PTHREAD)