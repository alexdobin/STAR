#ifndef SHAREDMEMORYSEGMENT_H
#define SHAREDMEMORYSEGMENT_H

class SharedMemorySegment
{
public:
	SharedMemorySegment(key_t key);
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
	
private :

	// private fields
	key_t _key;
	bool _needsAllocation; 
	void* _mapped; 
	size_t _size;
	std::string _name;
	bool _isAllocator;

	// private methods
	void CreateAndInitSharedObject(size_t shmSize);
	void SharedMemorySegment::OpenIfExists();
	void MapSharedObjectToMemory(); 
};

#endif // SHAREDMEMORYSEGMENT_H
