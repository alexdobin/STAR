#ifndef SHAREDMEMORYSEGMENT_H
#define SHAREDMEMORYSEGMENT_H

#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <iostream>
#include <string>
#include <memory>

enum ErrorState 
{
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


class SharedMemoryException : public std::exception
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

	SharedMemoryException(ErrorState error) : _error(error)
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


class SharedMemorySegment
{
public:
	SharedMemorySegment(const key_t key, const bool unloadLast);
	~SharedMemorySegment();

	// public methods
	void Allocate(const size_t shmSize);
	void Clean();
	void* GetMapped() const
	{
		// Start address of usable shared memory, this equals mapped memory address incremented with room for use count. 
		return (void*)((char*)_pMapped + sizeof(size_t));
	}

	bool NeedsAllocation() const
	{
		return _needsAllocation; 
	}
	bool IsAllocator() const
	{
		return _isAllocator;
	}
	void SetErrorStream(std::ostream *err)
	{
		_pErr = err;
	};
	void ThrowError(ErrorState error, int detail = 0)
	{
		if (!_exception.HasError())
		{
			_exception.SetError(error, detail);
		}
		throw _exception;
	};

private :
	// private fields
	bool _needsAllocation; 
	void* _pMapped; 
	std::string _name;
	bool _isAllocator;
	bool _unloadLast;
	std::ostream* _pErr;
	SharedMemoryException _exception;
	size_t* _pUseCount; 

	boost::interprocess::shared_memory_object* _pShmMemObj; 
	boost::interprocess::mapped_region* _pMappedRegion;
	
	// private methods
	void CreateAndInitSharedObject(const size_t shmSize);
	void OpenIfExists();
	void MapSharedObjectToMemory(); 
};

#endif // SHAREDMEMORYSEGMENT_H
