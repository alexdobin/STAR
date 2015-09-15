#include "CrossPlatform.h"
#include "SharedMemorySegment.h"

using namespace boost::interprocess; 

SharedMemorySegment::SharedMemorySegment(key_t key, bool unloadLast) : _key(key),
														_unloadLast(unloadLast),
														_needsAllocation(true),
														_mapped(nullptr),
														_isAllocator(false),
														_length(nullptr),
														_shm_obj_ptr(nullptr),
														_mapped_region_ptr(nullptr)
{
	_name = std::to_string(_key); 
	OpenIfExists(); 
}

SharedMemorySegment::~SharedMemorySegment()
{
	if (_unloadLast)
	{
		Clean();
	}
	delete _shm_obj_ptr;
	delete _mapped_region_ptr;
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
			_mapped = nullptr; 
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
		if (ex.get_error_code() != already_exists_error)
		{
			delete shm_obj_ptr; 
			ThrowError(EOPENFAILED, ex.get_error_code());
		}
		_exception.SetError(EEXISTS, 0);
		return;
	}

	try
	{
		unsigned long long toReserve = (unsigned long long) shmSize + sizeof(unsigned long long);

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
		_shm_obj_ptr = new shared_memory_object(
				open_only						//only open
				, _name.c_str()					//name
				, read_write					//read-write mode
				);
	}
	catch (interprocess_exception &ex)
	{
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
		_mapped_region_ptr = new mapped_region(
			*_shm_obj_ptr                        //Memory-mappable object
			, read_write					     //read-write mode
			);
		
		_mapped = _mapped_region_ptr->get_address();
		if (_mapped == nullptr)
		{
			ThrowError(EMAPFAILED);
		}

		_length = (size_t*)_mapped; 
		*_length = _mapped_region_ptr->get_size();

		_needsAllocation = false;
	}

	catch (interprocess_exception &ex)
	{
		ThrowError(EMAPFAILED, ex.get_error_code());
		
	}
	
}

