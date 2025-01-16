#ifndef _MMAP_MALLOC_H
#define _MMAP_MALLOC_H 1

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

void* mmap_malloc(size_t size);
void mmap_free(void* addr);
void* mmap_calloc(size_t num_elements, size_t element_size);
void* mmap_realloc(void *addr, size_t size);
void* mmap_reallocarray(void *addr, size_t size, size_t count);

#ifdef __cplusplus
}
#endif

#endif /* mmap_malloc.h */
