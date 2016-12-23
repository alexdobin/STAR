// SharedMemory.cpp
// Gery Vessere - gvessere@illumina.com, gery@vessere.com
// An abstraction over both SysV and POSIX shared memory APIs

#include "SharedMemory.h"
#include <sstream>
#include <sys/mman.h>
#include <sys/shm.h>
#include <sys/stat.h>        /* For mode constants */
#include <fcntl.h>           /* For O_* constants */
#include <semaphore.h>
#include <errno.h>

#ifdef COMPILE_FOR_MAC
  //some Mac's idiosyncrasies: standard SHM libraries are very old and missing some definitions
  #define SHM_NORESERVE 0
#endif

using namespace std;

SharedMemory::SharedMemory(key_t key, bool unloadLast): _key(key), _counterKey(key+1), _unloadLast(unloadLast), _err(&cerr)
{
    _shmID = -1;
    _sharedCounterID = -1;
    _counterMem = 0;
    _mapped=NULL;
    _length = NULL;
    _sem=NULL;
    _isAllocator = false;
    _needsAllocation = true;

    EnsureCounter();
    OpenIfExists();
}

SharedMemory::~SharedMemory()
{
    try
    {
        int inUse = SharedObjectsUseCount()-1;
        Close();

        if (_unloadLast)
        {
            if (inUse > 0)
            {
                (*_err) << inUse << " other job(s) are attached to the shared memory segment, will not remove it." <<endl;
            }
            else
            {
                (*_err) << "No other jobs are attached to the shared memory segment, removing it."<<endl;
                Clean();
            }
        }
    }
    catch (const SharedMemoryException & exc)
    {
        try
        {
           Clean();
        }
        catch (...)
        {}
    }
}

void SharedMemory::Allocate(size_t shmSize)
{
    _exception.ClearError();

    if (!_needsAllocation)
        ThrowError(EALREADYALLOCATED);

    CreateAndInitSharedObject(shmSize);

    if (_exception.HasError() && _exception.GetErrorCode() != EEXISTS)
        throw _exception;

    _exception.ClearError(); // someone else came in first so retry open

    OpenIfExists();

    _isAllocator = true;
}

string SharedMemory::GetPosixObjectKey()
{
    ostringstream key;
    key << "/" << _key;
    return key.str();
}

string SharedMemory::CounterName()
{
    ostringstream counterName;
    counterName << "/shared_use_counter" << _key;
    return counterName.str();
}


void SharedMemory::CreateAndInitSharedObject(size_t shmSize)
{
    unsigned long long toReserve = (unsigned long long) shmSize + sizeof(unsigned long long);

#ifdef POSIX_SHARED_MEM
    _shmID=shm_open(GetPosixObjectKey().c_str(), O_CREAT | O_RDWR | O_EXCL, 0666);
#else
    _shmID=shmget(_key, toReserve, IPC_CREAT | IPC_EXCL | SHM_NORESERVE | 0666); //        _shmID = shmget(shmKey, shmSize, IPC_CREAT | SHM_NORESERVE | SHM_HUGETLB | 0666);
#endif

    if (_shmID == -1)
    {
        switch (errno)
        {
            case EEXIST:
                _exception.SetError(EEXISTS, 0);
                break;
            default:
                ThrowError(EOPENFAILED, errno);
        }
        return;
    }

#ifdef POSIX_SHARED_MEM
    int err = ftruncate(_shmID, toReserve);
    if (err == -1)
    {
        ThrowError(EFTRUNCATE);
    }
#endif
}

void SharedMemory::OpenIfExists()
{
    errno=0;
    if (_shmID < 0){
#ifdef POSIX_SHARED_MEM
        _shmID=shm_open(GetPosixObjectKey().c_str(), O_RDWR, 0);
#else
        _shmID=shmget(_key,0,0);
#endif
}
    bool exists=_shmID>=0;
    if (! (exists || errno == ENOENT))
        ThrowError(EOPENFAILED, errno); // it's there but we couldn't get a handle

    if (exists)
    {
        MapSharedObjectToMemory();

        _needsAllocation = false;
    }
}

#ifdef POSIX_SHARED_MEM
struct stat SharedMemory::GetSharedObjectInfo()
{
    struct stat buf;
    int err = fstat(_shmID, &buf);
    if (err == -1)
        ThrowError(EOPENFAILED, errno);

    return buf;
}
#endif

void SharedMemory::MapSharedObjectToMemory()
{
#ifdef POSIX_SHARED_MEM
    size_t size=0;
    struct stat buf = SharedMemory::GetSharedObjectInfo();
    size = (size_t) buf.st_size;
    _mapped = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_NORESERVE, _shmID, (off_t) 0);

    if (_mapped==((void *) -1))
        ThrowError(EMAPFAILED, errno);

    _length = (size_t *) _mapped;
    *_length = size;
#else
    _mapped= shmat(_shmID, NULL, 0);

    if (_mapped==((void *) -1))
        ThrowError(EMAPFAILED, errno);

    _length = (size_t *) _mapped;
#endif
}

void SharedMemory::Close()
{
    #ifdef POSIX_SHARED_MEM
    if (_mapped != NULL)
    {
        int ret = munmap(_mapped, (size_t) *_length);
        if (ret == -1)
            ThrowError(EMAPFAILED, errno);
        _mapped = NULL;
    }

    if (_shmID != -1)
    {
        int err = close(_shmID);
        _shmID=-1;
        if (err == -1)
            ThrowError(ECLOSE, errno);
    }

    #else
    if (_mapped != NULL)
    {
        shmdt(_mapped);
        _mapped = NULL;
    }
    #endif
}

void SharedMemory::Unlink()
{
    if (!_needsAllocation)
    {
        int shmStatus=-1;
    #ifdef POSIX_SHARED_MEM
        shmStatus = shm_unlink(GetPosixObjectKey().c_str());
    #else
        struct shmid_ds buf;
        shmStatus=shmctl(_shmID,IPC_RMID,&buf);
    #endif
        if (shmStatus == -1)
            ThrowError(EUNLINK, errno);

        _needsAllocation = true;
    }
}

void SharedMemory::Clean()
{
    Close();
    Unlink();
    RemoveSharedCounter();
}

void SharedMemory::EnsureCounter()
{
    if (_sharedCounterID < 0)
        _sharedCounterID=shmget(_counterKey,0,0);

    bool exists=_sharedCounterID>=0;

    if (!exists)
    {
        errno=0;
        _sharedCounterID=shmget(_counterKey, 1, IPC_CREAT | IPC_EXCL | SHM_NORESERVE | 0666);

        if (_sharedCounterID < 0)
            ThrowError(ECOUNTERCREATE, errno);
    }

    if (_counterMem == 0)
    {
        _counterMem = shmat(_sharedCounterID, NULL, 0);

        if (_counterMem==((void *) -1))
            ThrowError(EMAPFAILED, errno);
    }
}

void SharedMemory::RemoveSharedCounter()
{
    struct shmid_ds buf;
    int shmStatus=shmctl(_sharedCounterID,IPC_RMID,&buf);
    if (shmStatus == -1)
        ThrowError(ECOUNTERREMOVE, errno);
}

int SharedMemory::SharedObjectsUseCount()
{
    EnsureCounter();
    if (_sharedCounterID != -1)
    {
        struct shmid_ds shmStat;
        int shmStatus=shmctl(_sharedCounterID,IPC_STAT,&shmStat);
        if (shmStatus == -1)
            ThrowError(ECOUNTERUSE, errno);

        return shmStat.shm_nattch;
    }
    else
        return -1;
}
