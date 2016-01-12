#ifndef GLOBAL_VARIABLES_DEF
#define GLOBAL_VARIABLES_DEF

#include "ThreadControl.h"
extern Stats g_statsAll;
extern ThreadControl g_threadChunks;

void LockMutexOutSAM(); 
void UnlockMutexOutSAM();

void LockMutexBAMsortBins();
void UnlockMutexBAMsortBins();

void LockMutexStats();
void UnlockMutexStats();

void LockMutexLogMain();
void UnlockMutexLogMain();

void LockMutexInRead();
void UnlockMutexInRead();

#endif

