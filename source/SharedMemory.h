// SharedMemory.cpp
// Gery Vessere - gvessere@illumina.com, gery@vessere.com
// An abstraction over both SysV and POSIX shared memory APIs

#ifndef SHAREDMEMORY_H
#define SHAREDMEMORY_H


#include <string>
#include <semaphore.h>
#include <unistd.h>
#include <exception>
#include <iostream>


enum ErrorState {
ENONE,
ENOTALLOCATED,
ETRYAGAIN,
EALREADYALLOCATED,
EOPENFAILED,
EEXISTS,
EFTRUNCATE,
EMAPFAILED,
ECLOSE,
EUNLINK,
ECOUNTERCREATE,
ECOUNTERREMOVE,
ECOUNTERUSE
};


class SharedMemoryException: public std::exception
{
private:
    bool _hasError;
    ErrorState _error;
    int _errorDetail;

public:

    SharedMemoryException()
    {
        _hasError = false;
        _error = ENONE;
        _errorDetail = 0;
    };

    SharedMemoryException(ErrorState error): _error(error)
    {};

    ErrorState GetErrorCode() const
    {
        return _error;
    };

    int GetErrorDetail() const
    {
        return _errorDetail;
    }

    bool HasError() const
    {
        return _hasError;
    };

    void SetError(ErrorState error, int detail)
    {
        if (!_hasError)
        {
            _hasError = true;
            _error = error;
            _errorDetail = detail;
        }
    }

    void ClearError()
    {
        _hasError = false;
        _error = ENONE;
        _errorDetail = 0;
    };
};

class SharedMemory
{
public:
        void * GetMapped()
        {
            return (void *) ((char*) _mapped + sizeof(size_t));
        };

        size_t GetSize()
        {
            if (!_needsAllocation)
                return *_length - sizeof(size_t);

            _exception.SetError(ENOTALLOCATED, 0);
            return -1;
        };

        bool NeedsAllocation() const
        {
            return _needsAllocation;
        };

        // the owner is the first one that created the named shared memory segment
        bool IsAllocator() const
        {
            return _isAllocator;
        };

        int GetId() const
        {
            return _shmID;
        };

        bool HasError() const
        {
            return _exception.HasError();
        }

        void ThrowError(ErrorState error, int detail)
        {
            if (!_exception.HasError())
            {
                _exception.SetError(error, detail);
            }
            throw _exception;
        };

        void ThrowError(ErrorState error)
        {
            ThrowError(error, 0);
        };

        void SetErrorStream(std::ostream * err)
        {
            _err = err;
        };

        SharedMemory(key_t key, bool unloadLast);
        ~SharedMemory();
        void Allocate(size_t shmSize);
        void Clean();

private:
        SharedMemoryException _exception;

        int _shmID;
        int _sharedCounterID;
        void * _counterMem;

        void * _mapped;
        size_t * _length;
        sem_t * _sem;
        bool _isAllocator;
        bool _needsAllocation;

        key_t _key;
        key_t _counterKey;
        bool _unloadLast;
        std::ostream * _err;

        int SharedObjectsUseCount();
        void OpenIfExists();
        void CreateAndInitSharedObject(size_t shmSize);
        void MapSharedObjectToMemory();
        std::string GetPosixObjectKey();
        struct stat GetSharedObjectInfo();
        void Close();
        void Unlink();

        std::string CounterName();

        void EnsureCounter();
        void RemoveSharedCounter();
        void SharedUseIncrement();
        void SharedUseDecrement();
};

#endif
