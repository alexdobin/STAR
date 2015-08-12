#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <string>
#include "CrossPlatform.h"
#include "SharedMemorySegment.h"

using namespace boost::interprocess;

SharedMemorySegment::SharedMemorySegment(key_t key, bool unloadLast) : _key(key),
														_unloadLast(unloadLast),
														_needsAllocation(true),
														_mapped(nullptr),
														_size(0),
														_name(std::to_string(_key)),
														_isAllocator(false)
{
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
	try
	{
		std::string key = std::to_string(_key);
		shared_memory_object shm_obj
			(create_only					//only create
			, _name.c_str()					//name
			, read_write					//read-write mode
			);

		// TODO : Check why this size adjustment is done.
		unsigned long long toReserve = (unsigned long long) shmSize + sizeof(unsigned long long);

		// TODO : Check if this can throw different exception that we need to handle.
		shm_obj.truncate(toReserve);
	}
	catch (interprocess_exception &ex)
	{
		// revisit this, can we blindly return this error ? 
		_exception.SetError(EEXISTS, ex.get_error_code());
	}
	return; 
}

void SharedMemorySegment::OpenIfExists()
{
	errno = 0;
	try
	{
		std::string key = std::to_string(_key);
		shared_memory_object shm_obj
			(open_only						//only create
			, key.c_str()					//name
			, read_write					//read-write mode
			);

		mapped_region region
			(shm_obj                        //Memory-mappable object
			, read_write					//Access mode
			);
		_mapped = region.get_address();

		_needsAllocation = false;
	}
	catch (interprocess_exception &ex)
	{
		ThrowError(EOPENFAILED, ex.get_error_code());
	}
}

