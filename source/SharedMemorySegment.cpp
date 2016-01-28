#include "CrossPlatform.h"
#include "SharedMemorySegment.h"

using namespace boost::interprocess; 

SharedMemorySegment::SharedMemorySegment(const key_t key, const bool unloadLast) : _unloadLast(unloadLast),
														_needsAllocation(true),
														_pMapped(nullptr),
														_isAllocator(false),
														_pShmMemObj(nullptr),
														_pMappedRegion(nullptr),
														_pUseCount(nullptr)
{
	_name = std::to_string(key);
	OpenIfExists(); 
}

SharedMemorySegment::~SharedMemorySegment()
{
	if (_pMapped != nullptr) // _pMapped will be nullptr if Clean() was called before. 
	{
		if (_pUseCount != nullptr)
		{
			// Decrement use count
			size_t inUse = (*_pUseCount);
			(*_pUseCount) = --inUse; 

			if (_unloadLast)
			{
				if (inUse == 0)
				{
					(*_pErr) << "No other jobs are attached to the shared memory segment, removing it." << std::endl;
					shared_memory_object::remove(_name.c_str()); // this will fail only if shared memory does not exist, not expected.
				}
				else
				{
					(*_pErr) << inUse << " other job(s) are attached to the shared memory segment, will not remove it." << std::endl;
				}
			}
			_pUseCount = nullptr; 
		}
		_pMapped = nullptr; 
	}

	delete _pMappedRegion;
	delete _pShmMemObj;
}

void SharedMemorySegment::Allocate(const size_t shmSize)
{
	_exception.ClearError(); 
	
	if (!_needsAllocation)
		ThrowError(EALREADYALLOCATED);

	CreateAndInitSharedObject(shmSize);

	if (_exception.HasError() && _exception.GetErrorCode() != EEXISTS)
		throw _exception;

	_exception.ClearError();

	OpenIfExists();

	_isAllocator = true;
}

void SharedMemorySegment::Clean()
{
	if (!_needsAllocation)
	{
		try
		{
			if (_pMapped != nullptr)
			{
				if (_pUseCount != nullptr)
				{
					// Decrement use count
					size_t inUse = (*_pUseCount);
					(*_pUseCount) = --inUse;

					// try to remove shared memory
					if (!shared_memory_object::remove(_name.c_str()))
					{
						// Shared Memory does not exist. 
						ThrowError(EUNLINK, errno);
					}
					_pUseCount = nullptr;
				}
				_pMapped = nullptr;
			}
		}
		catch (interprocess_exception &ex)
		{
			ThrowError(EUNLINK, ex.get_error_code());
		}
		// successfully removed, set flag.
		_needsAllocation = true;
	}
}

void SharedMemorySegment::CreateAndInitSharedObject(const size_t shmSize)
{
	shared_memory_object* shm_obj_ptr = nullptr;
	try
	{
		shm_obj_ptr = new shared_memory_object(
				create_only						//only create
				, _name.c_str()					//name
				, read_write					//read-write mode
				);
	}
	catch (interprocess_exception &ex)
	{	
		delete shm_obj_ptr; // delete pointer, we are going to return.
		if (ex.get_error_code() != already_exists_error)
		{
			ThrowError(EOPENFAILED, ex.get_error_code());
		}
		// Already exist, set error code and return.
		_exception.SetError(EEXISTS, 0);
		return;
	}

	try
	{
		// Total size of Shared Memory = shmsize + sizeof(size_t) for use count.
		size_t toReserve = shmSize + sizeof(size_t);

		shm_obj_ptr->truncate(toReserve);
	}
	catch (interprocess_exception &ex)
	{
		delete shm_obj_ptr;
		ThrowError(EFTRUNCATE, ex.get_error_code());
	}
	delete shm_obj_ptr;
	return; 
}

void SharedMemorySegment::OpenIfExists()
{
	try
	{
		_pShmMemObj = new shared_memory_object(
				open_only						//only open
				, _name.c_str()					//name
				, read_write					//read-write mode
				);
	}
	catch (interprocess_exception &ex)
	{
		if (ex.get_error_code() != not_found_error) //check if error is something other than not found, we need to throw it. 
		{
			ThrowError(EOPENFAILED, ex.get_error_code());
		}
		// Shared Memory not found, return. 
		return; 

	}

	// Shared Memory exist, map it and get address.
	try
	{
		_pMappedRegion = new mapped_region(
			*_pShmMemObj                         //Memory-mappable object
			, read_write					     //read-write mode
			);
		
		_pMapped = _pMappedRegion->get_address();
		if (_pMapped == nullptr)
		{
			ThrowError(EMAPFAILED);
		}

		// Get pointer to use count
		_pUseCount = (size_t*)_pMapped;

		size_t inUse = (*_pUseCount); 

		// shared memory is initialized to 0 when created, so it is safe to 
		// just increment use count even if this is first time we open it. 
		(*_pUseCount) = ++inUse; 

		_needsAllocation = false;
	}

	catch (interprocess_exception &ex)
	{
		ThrowError(EMAPFAILED, ex.get_error_code());
	}
}

