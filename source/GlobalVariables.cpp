#include "GlobalVariables.h"
Stats g_statsAll;//global mapping statistics
ThreadControl g_threadChunks;

void LockMutexOutSAM()
{
#if !defined(_WIN32) && defined(USE_PTHREAD)
	pthread_mutex_lock(&g_threadChunks.mutexOutSAM) ; 
#else
	g_threadChunks.mutexOutSAM.lock();
#endif
}

void UnlockMutexOutSAM()
{
#if !defined(_WIN32) && defined(USE_PTHREAD)
	pthread_mutex_unlock(&g_threadChunks.mutexOutSAM) ; 
#else
	g_threadChunks.mutexOutSAM.unlock();
#endif
}

void LockMutexBAMsortBins()
{
#if !defined(_WIN32) && defined(USE_PTHREAD)
	pthread_mutex_lock(&g_threadChunks.mutexBAMsortBins);
#else
	g_threadChunks.mutexBAMsortBins.lock();
#endif
}

void UnlockMutexBAMsortBins()
{
#if !defined(_WIN32) && defined(USE_PTHREAD)
	pthread_mutex_unlock(&g_threadChunks.mutexBAMsortBins);
#else
	g_threadChunks.mutexBAMsortBins.unlock();
#endif
}

void LockMutexStats()
{
#if !defined(_WIN32) && defined(USE_PTHREAD)
	pthread_mutex_lock(&g_threadChunks.mutexStats);
#else
	g_threadChunks.mutexStats.lock();
#endif
}

void UnlockMutexStats()
{
#if !defined(_WIN32) && defined(USE_PTHREAD)
	pthread_mutex_unlock(&g_threadChunks.mutexStats);
#else
	g_threadChunks.mutexStats.unlock();
#endif
}

void LockMutexLogMain()
{
#if !defined(_WIN32) && defined(USE_PTHREAD)
	pthread_mutex_lock(&g_threadChunks.mutexLogMain);
#else
	g_threadChunks.mutexLogMain.lock();
#endif
}

void UnlockMutexLogMain()
{
#if !defined(_WIN32) && defined(USE_PTHREAD)
	pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
#else
	g_threadChunks.mutexLogMain.unlock();
#endif
}

void LockMutexInRead()
{
#if !defined(_WIN32) && defined(USE_PTHREAD)
	pthread_mutex_lock(&g_threadChunks.mutexInRead);
#else
	g_threadChunks.mutexInRead.lock();
#endif
}

void UnlockMutexInRead()
{
#if !defined(_WIN32) && defined(USE_PTHREAD)
	pthread_mutex_unlock(&g_threadChunks.mutexInRead);
#else
	g_threadChunks.mutexInRead.unlock();
#endif
}