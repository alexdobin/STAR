#ifndef MMAP_MGR_H
#define MMAP_MGR_H

#include <stdlib.h>

void* mmap_reserve(const size_t size);
int mmap_maptemp(void* addr, const size_t size, char* template);
int mmap_unmap(void* addr, const size_t size);

#endif