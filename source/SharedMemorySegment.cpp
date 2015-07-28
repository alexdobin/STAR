#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <string>
#include "CrossPlatform.h"
#include "SharedMemorySegment.h"

using namespace boost::interprocess;

SharedMemorySegment::SharedMemorySegment(key_t key) : _key(key),
														_needsAllocation(true),
														_mapped(nullptr),
														_size(0),
														_name(std::to_string(_key))
{

}

SharedMemorySegment::~SharedMemorySegment()
{
	Clean(); 
}

void SharedMemorySegment::Allocate(size_t shmSize)
{
	CreateAndInitSharedObject(shmSize);
}

void SharedMemorySegment::Clean()
{
	try
	{
		if (!shared_memory_object::remove(_name.c_str()))
		{
			// TODO : Throw Exception.
		}
	}
	catch (interprocess_exception &ex)
	{

	}
}

void SharedMemorySegment::CreateAndInitSharedObject(size_t shmSize)
{
	
	try
	{
		if (!_needsAllocation)
		{
			//TODO : throw exception.
			return;
		}
		std::string key = std::to_string(_key);
		shared_memory_object shm_obj
			(create_only					//only create
			, _name.c_str()					//name
			, read_write					//read-write mode
			);
		// TODO : Check adjustements done in size
		shm_obj.truncate(shmSize);
	}
	catch (interprocess_exception &ex)
	{
		//TODO:
	}

}

void SharedMemorySegment::OpenIfExists()
{
	try
	{
		std::string key = std::to_string(_key);
		shared_memory_object shm_obj
			(open_only					//only create
			, key.c_str()					//name
			, read_write					//read-write mode
			);

		// TODO : Check if exist and Map.

		mapped_region region
			(shm_obj                      //Memory-mappable object
			, read_write               //Access mode
			);
		_mapped = region.get_address();
	}
	catch (interprocess_exception &ex)
	{
		//TODO :
	}

}

