#include "SharedMemory.h"

#include <iostream>
#include <sstream>
#include <sys/mman.h>
#include <sys/shm.h>
#include <sys/stat.h>        /* For mode constants */
#include <fcntl.h>           /* For O_* constants */
#include <semaphore.h>

using namespace std;

SharedMemory::SharedMemory(key_t key, bool unloadLast): _key(key), _unloadLast(unloadLast)
{
    _needsAllocation = false;

    OpenIfExists();
}

SharedMemory::~SharedMemory()
{
    try
    {
        int inUse = SharedObjectsUseCount() - 1;
        SharedUseDecrement();
        Close();

        if (inUse > 0)
        {
            cerr << inUse << " other job(s) are attached to the shared memory segment, will not remove it." <<endl;
        } 
        else if (_unloadLast)
        {
            cerr << "No other jobs are attached to the shared memory segment, removing it."<<endl;
            Unlink();
            RemoveSharedCounter();
        };
    }
    catch (const SharedMemoryException & exc)
    {
        Unlink();
        RemoveSharedCounter();
    }
}

void SharedMemory::Allocate(size_t shmSize)
{
    ClearError();

    if (!_needsAllocation)
        ThrowError(ErrorState::EALREADYALLOCATED);

    try
    {
        CreateAndInitSharedObject(shmSize);
    }
    catch (SharedMemoryException & exc)
    {
        if (exc.GetErrorCode() != ErrorState::EEXISTS)
            rethrow_exception(current_exception());
        
        ClearError(); // someone else came in first
    }

    OpenIfExists();
    
    _isAllocator = true;

    //TODO: Error validation
}

const char * SharedMemory::GetPosixObjectKey()
{
    ostringstream key;
    key << "/" << _key;
    return key.str().c_str();
}

string SharedMemory::CounterName()
{
    ostringstream counterName;
    counterName << "/shared_use_counter" << _shmID;
    return counterName.str();
}


void SharedMemory::CreateAndInitSharedObject(size_t shmSize)
{    
    unsigned long long toReserve = (unsigned long long) shmSize + sizeof(unsigned long long);

#ifdef POSIX_SHARED_MEM
    _shmID=shm_open(GetPosixObjectKey(), O_CREAT | O_RDWR | O_EXCL  , 0666);
#else
    _shmID = shmget(_key, toReserve, IPC_CREAT | IPC_EXCL | SHM_NORESERVE | 0666); //        _shmID = shmget(shmKey, shmSize, IPC_CREAT | SHM_NORESERVE | SHM_HUGETLB | 0666);
#endif

    if (_shmID == -1 && errno == EEXIST)
        ThrowError(ErrorState::EEXISTS);

#ifdef POSIX_SHARED_MEM
    int err = ftruncate(_shmID, toReserve);
    if (err == -1)
    {
        Destroy();
        ThrowError(ErrorState::EFTRUNCATE);
    }
#endif
}

void SharedMemory::OpenIfExists()
{
#ifdef POSIX_SHARED_MEM
    _shmID=shm_open(SharedMemory::GetPosixObjectKey(), O_RDWR, 0);
#else
    _shmID=shmget(_key,0,0);
#endif

    bool exists=(_shmID!=-1);

    if (!exists && errno !=ENOENT)
        ThrowError(EOPENFAILED);

    if (exists)
    {
        // mandatory since on first re-open we don't have the size of object
        struct stat buf = SharedMemory::GetSharedObjectInfo(); 
        size_t size = (size_t) buf.st_size;

        MapSharedObjectToMemory(size);
        
        _needsAllocation = false;
        SharedUseIncrement();
    }
}

struct stat SharedMemory::GetSharedObjectInfo()
{
    struct stat buf;
    int err = fstat(_shmID, &buf);
    if (err == -1)
        ThrowError(EOPENFAILED);

    return buf;
}

void SharedMemory::MapSharedObjectToMemory(size_t size)
{
#ifdef POSIX_SHARED_MEM
    _mapped = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_NORESERVE, _shmID, (off_t) 0);
#else
    _mapped= shmat(_shmID, NULL, 0);
#endif

    if (_mapped==((void *) -1)) 
        ThrowError(EMAPFAILED);
    // set the length field
    _length = (size_t *) _mapped;
    *_length = size;
}

void SharedMemory::Close()
{
    if (_mapped != NULL)
    {
        int ret = munmap(_mapped, (size_t) _length);
        _mapped = NULL;
        if (ret == -1)
            ThrowError(EMAPFAILED);
    }

    if (_shmID != 0)
    {
        int err = close(_shmID);
        if (err == -1)
            ThrowError(ECLOSE);
    }
}

void SharedMemory::Unlink()
{
#ifdef POSIX_SHARED_MEM
    int ret = shm_unlink(SharedMemory::GetPosixObjectKey());
    if (ret == -1)
        ThrowError(EUNLINK);

#else
    struct shmid_ds *buf=NULL;
    int shmStatus=shmctl(_shmID,IPC_RMID,buf);
    if (shmStatus==-1) 
        ThrowError(EUNLINK);

#endif
}

void SharedMemory::Destroy()
{
    Close();
    Unlink();
    RemoveSharedCounter();
}

void SharedMemory::EnsureCounter()
{
    if (_sem!= NULL)
        return;

    const char * counterName = SharedMemory::CounterName().c_str();
    _sem = sem_open(counterName, 0);
    if ((errno & ENOENT) == ENOENT)
    {
        _sem = sem_open(counterName, O_CREAT, 0666, 0);
        if (errno != 0)
            ThrowError(ECOUNTERCREATE);
    }
}

void SharedMemory::RemoveSharedCounter()
{
    const char * counterName = SharedMemory::CounterName().c_str();
    int ret = sem_unlink(counterName);
    if (ret == -1)
            ThrowError(ECOUNTERREMOVE);
}

void SharedMemory::SharedUseIncrement()
{
    SharedMemory::EnsureCounter();
    int ret = sem_post(_sem);
    if (ret == -1)
        ThrowError(ECOUNTERINC);
}

void SharedMemory::SharedUseDecrement()
{
    SharedMemory::EnsureCounter();
    int ret = sem_trywait(_sem);
    if (ret == -1)
        ThrowError(ECOUNTERDEC);
}

int SharedMemory::SharedObjectsUseCount()
{
    SharedMemory::EnsureCounter();
    int sval=-1;
    int ret = sem_getvalue(_sem, &sval);

    if (ret == -1 || sval == -1)
        ThrowError(ECOUNTERUSE);

    return sval;
}
