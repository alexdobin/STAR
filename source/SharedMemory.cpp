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
    _shmID = -1;
    _mapped=NULL;
    _length = NULL;
    _sem=NULL;
    _isAllocator = false;
    _needsAllocation = true;

    OpenIfExists();
}

SharedMemory::~SharedMemory()
{
    try
    {
        int inUse = SharedObjectsUseCount() - 1;
        Close();

        if (inUse > 0 && _unloadLast)
        {
            cerr << inUse << " other job(s) are attached to the shared memory segment, will not remove it." <<endl;
        } 
        else if (_unloadLast)
        {
            cerr << "No other jobs are attached to the shared memory segment, removing it."<<endl;
            Clean();
        };
    }
    catch (const SharedMemoryException & exc)
    {
        //cerr << exc.GetErrorCode() << " " << exc.GetErrorDetail() << endl;
        Clean();
    }
}

void SharedMemory::Allocate(size_t shmSize)
{
    _exception.ClearError();

    if (!_needsAllocation)
        ThrowError(ErrorState::EALREADYALLOCATED);
    
    CreateAndInitSharedObject(shmSize);

    if (_exception.HasError() && _exception.GetErrorCode() != ErrorState::EEXISTS)
        throw _exception;

    _exception.ClearError(); // someone else came in first so retry open

    OpenIfExists();
    
    _isAllocator = true;
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
    counterName << "/shared_use_counter" << _key;
    return counterName.str();
}


void SharedMemory::CreateAndInitSharedObject(size_t shmSize)
{    
    unsigned long long toReserve = (unsigned long long) shmSize + sizeof(unsigned long long);

#ifdef POSIX_SHARED_MEM
    _shmID=shm_open(GetPosixObjectKey(), O_CREAT | O_RDWR | O_EXCL  , 0666);
#else
    _shmID=shmget(_key, toReserve, IPC_CREAT | IPC_EXCL | SHM_NORESERVE | 0666); //        _shmID = shmget(shmKey, shmSize, IPC_CREAT | SHM_NORESERVE | SHM_HUGETLB | 0666);
#endif

    if (_shmID == -1)
    {
        switch (errno)
        {
            case EEXIST:
                _exception.SetError(ErrorState::EEXISTS, 0);
                break;
            default:
                ThrowError(ErrorState::EOPENFAILED, errno);
        }
        return;
    }

#ifdef POSIX_SHARED_MEM
    int err = ftruncate(_shmID, toReserve);
    if (err == -1)
    {
        ThrowError(ErrorState::EFTRUNCATE);
    }
#endif
}

void SharedMemory::OpenIfExists()
{
    errno=0;
    if (_shmID < 0){
#ifdef POSIX_SHARED_MEM
        _shmID=shm_open(SharedMemory::GetPosixObjectKey(), O_RDWR, 0);
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
#else

    // size is not needed
    _mapped= shmat(_shmID, NULL, 0);
#endif

    if (_mapped==((void *) -1)) 
        ThrowError(EMAPFAILED, errno);
    
    SharedUseIncrement();

    // set the length field
    _length = (size_t *) _mapped;
#ifdef POSIX_SHARED_MEM
    *_length = size;
#endif
}

void SharedMemory::Close()
{
    if (_mapped != NULL)
    {
    #ifdef POSIX_SHARED_MEM
        int ret = munmap(_mapped, (size_t) *_length);
        if (ret == -1)
            ThrowError(EMAPFAILED, errno);
    #endif
        _mapped = NULL;
        SharedUseDecrement();
    }

    #ifdef POSIX_SHARED_MEM
    if (_shmID != -1)
    {
        int err = close(_shmID);
        _shmID=-1;
        if (err == -1)
            ThrowError(ECLOSE, errno);
    }
    #endif
}

void SharedMemory::Unlink()
{
    if (!_needsAllocation)
    {
        int shmStatus=-1;
    #ifdef POSIX_SHARED_MEM
        shmStatus = shm_unlink(SharedMemory::GetPosixObjectKey());

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
    if (_sem!= NULL)
        return;

    const char * counterName = SharedMemory::CounterName().c_str();

    _sem = sem_open(counterName, O_CREAT, 0666, 0);
    if (_sem == SEM_FAILED)
        ThrowError(ECOUNTERCREATE, errno);
}

void SharedMemory::RemoveSharedCounter()
{
    const char * counterName = SharedMemory::CounterName().c_str();

    if (_sem != NULL)
    {
        int ret = sem_close(_sem);
        if (ret == -1)
            ThrowError(ECLOSE, errno);
        
        ret = sem_unlink(counterName);
        if (ret == -1)
            ThrowError(ECOUNTERREMOVE, errno);
        _sem = NULL;
    }
}

void SharedMemory::SharedUseIncrement()
{
    SharedMemory::EnsureCounter();
    int ret = sem_post(_sem);
    if (ret == -1)
        ThrowError(ECOUNTERINC, errno);
    cerr << "incremented shared memory fragment usage to " << SharedObjectsUseCount() << endl;
}

void SharedMemory::SharedUseDecrement()
{
    SharedMemory::EnsureCounter();
    int ret = sem_trywait(_sem);
    if (ret == -1 && errno != EAGAIN)
        ThrowError(ECOUNTERDEC, errno);
    cerr << "decremented shared memory fragment usage to " << SharedObjectsUseCount() << endl;
}

int SharedMemory::SharedObjectsUseCount()
{
    SharedMemory::EnsureCounter();
    int sval=-1;
    int ret = sem_getvalue(_sem, &sval);

    if (ret == -1 || sval == -1)
        ThrowError(ECOUNTERUSE, errno);

    return sval;
}
