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
private :
	key_t _key;
	bool _needsAllocation; 
	void* _mapped; 
	size_t _size;
	std::string _name; 

	// private methods
	void CreateAndInitSharedObject(size_t shmSize);
	void SharedMemorySegment::OpenIfExists();
	void MapSharedObjectToMemory(); 
};

#endif 
