#ifndef SHAREDMEMORYSEGMENT_H
#define SHAREDMEMORYSEGMENT_H

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

private :

	// private fields
	key_t _key;
	bool _needsAllocation; 
	void* _mapped; 
	size_t _size;
	std::string _name;
	bool _isAllocator;
	bool _unloadLast;
	std::ostream* _err;

	// private methods
	void CreateAndInitSharedObject(size_t shmSize);
	void SharedMemorySegment::OpenIfExists();
	void MapSharedObjectToMemory(); 
};

#endif // SHAREDMEMORYSEGMENT_H
