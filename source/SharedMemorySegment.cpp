#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <string>
#include <malloc.h>
#include "CrossPlatform.h"
#include "SharedMemorySegment.h"

using namespace boost::interprocess;

SharedMemorySegment::SharedMemorySegment(key_t key, bool unloadLast) : _key(key),
														_unloadLast(unloadLast),
														_needsAllocation(true),
														_mapped(nullptr),
														_isAllocator(false),
														_length(nullptr)
{
	_name = std::to_string(_key); 
	OpenIfExists(); 
}

SharedMemorySegment::~SharedMemorySegment()
{
	Clean(); 
}

void SharedMemorySegment::Allocate(size_t shmSize)
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

void SharedMemorySegment::Clean()
{
	try
	{
		if (_mapped != nullptr)
		{
			if (!shared_memory_object::remove(_name.c_str()))
			{
				ThrowError(EUNLINK, errno);
			}
		}
	}
	catch (interprocess_exception &ex)
	{
		ThrowError(EUNLINK, ex.get_error_code());
	}
	_needsAllocation = true; 
}

void SharedMemorySegment::CreateAndInitSharedObject(size_t shmSize)
{
	std::unique_ptr<shared_memory_object> shm_obj_ptr;
	try
	{
		std::string key = std::to_string(_key);
		shm_obj_ptr = std::make_unique<shared_memory_object>(
			create_only						//only create
			, _name.c_str()					//name
			, read_write					//read-write mode
			);
	}
	catch (interprocess_exception &ex)
	{
		if (ex.get_error_code() != already_exists_error)
		{
			ThrowError(EOPENFAILED, ex.get_error_code());
		}
		_exception.SetError(EEXISTS, 0);
		return;
	}

	try
	{
		// TODO : Check why this size adjustment is done.
		unsigned long long toReserve = (unsigned long long) shmSize + sizeof(unsigned long long);

		(*shm_obj_ptr).truncate(toReserve);
	}
	catch (interprocess_exception &ex)
	{
		// revisit this, can we blindly return this error ? 
		ThrowError(EFTRUNCATE, ex.get_error_code());
	}
	return; 
}

void SharedMemorySegment::OpenIfExists()
{
	errno = 0;
	std::unique_ptr<shared_memory_object> shm_obj_ptr;

	try
	{
		std::string key = std::to_string(_key);
		shm_obj_ptr = std::make_unique<shared_memory_object>(
			open_only						//only create
			, key.c_str()					//name
			, read_write					//read-write mode
			);
	}
	catch (interprocess_exception &ex)
	{
		// TODO : TCD, check if this return same error code in linux too.
		if (ex.get_error_code() != not_found_error) // Shared Memory does not exist.
		{
			ThrowError(EOPENFAILED, ex.get_error_code());
		}
		// Shared Memory not found, return. 
		return; 

	}

	// Shared Memory exist,  map that and get address.
	try
	{
		mapped_region region(
			*shm_obj_ptr                        //Memory-mappable object
			, read_write					    //Access mode
			);
		
		_mapped = region.get_address();
		if (_mapped == nullptr)
		{
			ThrowError(EMAPFAILED);
		}

		_length = (size_t*)_mapped; 
		*_length = region.get_size(); 

		_needsAllocation = false;
	}

	catch (interprocess_exception &ex)
	{
		// TODO : TCD, check if this return same error code in linux too.
		int error_code = ex.get_error_code(); 
		error_code = ex.get_native_error();
		ThrowError(EMAPFAILED, ex.get_error_code());
		
	}
	
}

