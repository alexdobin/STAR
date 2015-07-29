#ifndef THREAD_CONTROL_DEF
#define THREAD_CONTROL_DEF

#include "ReadAlignChunk.h"
#include "CrossPlatform.h"

#if !defined(_WIN32) && defined(USE_PTHREAD)
#include <pthread.h>
#else
#include <thread>
#include <mutex>
#endif

#define MAX_chunkOutBAMposition 100000

class ThreadControl {
public:
    bool threadBool;
#ifdef _WIN32
	std::vector<std::thread> threadArray;
	std::mutex mutexInRead, mutexOutSAM, mutexOutBAM1, mutexOutChimSAM, mutexOutChimJunction, mutexOutUnmappedFastx, mutexOutFilterBySJout;
	std::mutex mutexStats, mutexLogMain, mutexBAMsortBins;
#else
	pthread_t *threadArray;
	pthread_mutex_t mutexInRead, mutexOutSAM, mutexOutBAM1, mutexOutChimSAM, mutexOutChimJunction, mutexOutUnmappedFastx, mutexOutFilterBySJout;
	pthread_mutex_t mutexStats, mutexLogMain, mutexBAMsortBins;
#endif

   
    
    uint chunkInN,chunkOutN;
    
    ThreadControl();
    
    static void* threadRAprocessChunks(void *RAchunk) {
        ( (ReadAlignChunk*) RAchunk )->processChunks();

#if !defined(_WIN32) && defined(USE_PTHREAD)
        pthread_exit(0);
#else
		// The standard method to exit std::thread is just exit the thread function. 
		return nullptr; 
#endif

        return NULL;
    };
};

#endif

