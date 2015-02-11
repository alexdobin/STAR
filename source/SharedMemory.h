
#ifndef SHAREDMEMORY_H
#define SHAREDMEMORY_H


#include <string>
#include <semaphore.h>
#include <unistd.h>
#include <exception>



enum ErrorState { 
ENONE,
ENOTALLOCATED,
ETRYAGAIN,
EALREADYALLOCATED, 

EOPENFAILED,        	 
// errOut << "EXITING because of FATAL ERROR: problems with shared memory: error from shmget() or shm_open(): " << strerror(errno) << "\n" << flush;
// errOut << "SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n" <<flush;     
// exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_SHM, *P);

EEXISTS,
// ostringstream errOut;
// errOut <<"EXITING: fatal error from shmget() trying to allocate shared memory piece: error type: " << strerror(errno) <<"\n";
// errOut <<"Possible cause 1: not enough RAM. Check if you have enough RAM of at least " << shmSize+2000000000 << " bytes\n";
// errOut <<"Possible cause 2: not enough virtual memory allowed with ulimit. SOLUTION: run ulimit -v " <<  shmSize+2000000000 <<"\n";
// errOut <<"Possible cause 3: allowed shared memory size is not large enough. SOLUTIONS: (i) consult STAR manual on how to increase shared memory allocation; "
// "(ii) ask your system administrator to increase shared memory allocation; (iii) run STAR with --genomeLoad NoSharedMemory\n"<<flush;
// exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_MEMORY_ALLOCATION, *P);  

EFTRUNCATE,


// ostringstream errOut;
// errOut <<"EXITING: fatal error from ftruncate() error shared memory: error type: " << strerror(errno) << endl;
// errOut <<"Possible cause 1: not enough RAM. Check if you have enough RAM of at least " << shmSize+2000000000 << " bytes\n";
// void * ptr = NULL;

// exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_MEMORY_ALLOCATION, *P); 

EMAPFAILED,


// ostringstream errOut;
// errOut << "EXITING because of FATAL ERROR: problems with shared memory: error from shmat() while trying to get address of the shared memory piece: " << strerror(errno) << "\n" <<flush;
// errOut << "SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n" <<flush;     
// exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_SHM, *P);

// ostringstream errOut;
// errOut <<"EXITING because of FATAL ERROR:  could not delete the shared object: " << strerror(errno) <<flush;
// exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_SHM, *P);


ECLOSE,

// ostringstream errOut;
// errOut << "EXITING because of FATAL ERROR: could not close the shared memory object: " << strerror(errno) << "\n" <<flush;     
// exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_SHM, *P);

EUNLINK,

// ostringstream errOut;
// errOut <<"EXITING because of FATAL ERROR:  could not delete the shared object: " << strerror(errno) <<flush;
// exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_SHM, *P);


// ostringstream errOut;
// errOut << "EXITING because of FATAL ERROR: problems with shared memory: error from shmctl() while trying to remove shared memory piece: " << strerror(errno) << "\n" <<flush;
// errOut << "SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory\n" <<flush;     
// exitWithError(errOut.str(),std::cerr, P->inOut->logMain, EXIT_CODE_SHM, *P);    


ECOUNTERCREATE,
ECOUNTERREMOVE,
ECOUNTERINC,
ECOUNTERDEC,
ECOUNTERUSE


};


class SharedMemoryException: public std::exception
{
private:
	ErrorState _error;

public:

	SharedMemoryException(int error): _error((ErrorState)error)
	{};

	ErrorState GetErrorCode()
	{
		return _error;
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
			
			SetError(ENOTALLOCATED);
			return -1;
		};

		bool NeedsAllocation()
		{
			return _needsAllocation;
		};

		// the owner is the first one that created the named shared memory segment
		bool IsAllocator()
		{
			return _isAllocator;
		};

		bool HasError()
		{
			return _hasError;
		};

		void SetError(ErrorState error)
		{
			if (!_hasError)
			{
			    _hasError = true;
			    _error = error;
			}
		}

		void ThrowError(ErrorState error)
		{
			SetError(error);
			throw new SharedMemoryException(_error);
		};

		ErrorState GetErrorCode()
		{
			return _error;
		};

		void ClearError()
		{
		    _hasError = false;
		    _error = ErrorState::ENONE;
		};

		int GetId()
		{
			return _shmID;
		};

		SharedMemory(key_t key, bool unloadLast);
		~SharedMemory();
		void Allocate(size_t shmSize);
		void Destroy();

private:
		ErrorState _error;

        int _shmID = 0;
        void * _mapped=NULL;
		size_t * _length = NULL;
        sem_t * _sem=NULL;
        bool _isAllocator = false;
        bool _needsAllocation = true;
		bool _hasError = false;

		key_t _key;
        bool _unloadLast = true;
		
		int SharedObjectsUseCount();
		void OpenIfExists();
        void CreateAndInitSharedObject(size_t shmSize);
        void MapSharedObjectToMemory(size_t size);
        const char * GetPosixObjectKey();
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