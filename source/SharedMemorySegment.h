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
	SharedMemorySegment(key_t key, bool unloadLast);
	~SharedMemorySegment();

	// public methods
	void Allocate(size_t shmSize);
	void Clean();
	void* GetMapped() const
	{
		return (void*)((char*)_mapped + sizeof(size_t));
	}

	bool NeedsAllocation() const
	{
		return _needsAllocation; 
	}
	bool IsAllocator() const
	{
		return _isAllocator;
	}
	void SetErrorStream(std::ostream * err)
	{
		_err = err;
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
	key_t _key;
	bool _needsAllocation; 
	void* _mapped; 
	std::string _name;
	bool _isAllocator;
	bool _unloadLast;
	std::ostream* _err;
	SharedMemoryException _exception;
	size_t* _length; 
	std::unique_ptr<boost::interprocess::shared_memory_object> _shm_obj_ptr;
	std::unique_ptr<boost::interprocess::mapped_region> _mapped_region_ptr;

	// private methods
	void CreateAndInitSharedObject(size_t shmSize);
	void OpenIfExists();
	void MapSharedObjectToMemory(); 
};

#endif // SHAREDMEMORYSEGMENT_H
