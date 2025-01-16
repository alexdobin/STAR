#define _GNU_SOURCE // To enable various non-standard GNU extensions.
#include <assert.h>
#include <dlfcn.h>
#include <errno.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "constants.h"
#include "default_config.h"
#include "heap.h"
#include "list.h"
#include "mmap_mgr.h"
#include <stdlib.h>
#include "profiling.h"

#if defined(__APPLE__)
#include <TargetConditionals.h>
#if TARGET_OS_OSX
#include <malloc/malloc.h>
#define MALLOCSIZE malloc_size
#else
#error "Unsupported Apple platform"
#endif
#elif __linux__
#include <malloc.h>
#define MALLOCSIZE malloc_usable_size
#else
#error "Unsupported compiler"
#endif

/*---------------------------------------------------------------------------*/
// Out of memory.
#define OUT_OF_MEMORY() do { \
  fprintf(stderr, "Out of space. Failed to allocate more memory.\n"); \
  errno = ENOMEM; \
} while (0)

/*---------------------------------------------------------------------------*/
// Allocator status
enum allocator_status {
  FAILED_TO_LOAD = -1,
  NOT_LOADED = 0,
  LOADED = 1
};
static enum allocator_status allocator_status = NOT_LOADED;

/*---------------------------------------------------------------------------*/
// Global lock.
static pthread_mutex_t glock = PTHREAD_MUTEX_INITIALIZER;

#define GLOBAL_LOCK_ACQUIRE() do {\
  pthread_mutex_lock(&glock);\
} while (0)

#define GLOBAL_LOCK_RELEASE() do {\
  pthread_mutex_unlock(&glock);\
} while (0)

/*---------------------------------------------------------------------------*/
// Parameters

// Fixed
static size_t page_size = 0;
static void* mmap_region_base = NULL;
static void* mmap_region_end = NULL;

// Configurable
static size_t mmap_heap_size = 0;
static size_t mmap_alloctor_min_bsize = 0;
static char naming_template[MAX_NAMING_TEMPLATE_SIZE + 1] = {0};

#define CEILING_PAGE_SIZE(size) (((size) + (page_size - 1)) & ~(page_size - 1))

/*---------------------------------------------------------------------------*/
// Heap
heap_t mmap_heap;

// Wrapper for heap check. 
// Notice that heap check should acquire the global heap lock.
// Be aware of dead lock. Don't die!
#if DEBUG_HEAP
#define HEAP_CHECK() do {\
  GLOBAL_LOCK_ACQUIRE();\
  heap_check(&mmap_heap);\
  GLOBAL_LOCK_RELEASE();\
} while (0)
#else
#define HEAP_CHECK() do {} while (0)
#endif

bool LOCAL_HELPER
config_parameters() {
  assert(page_size && "Page size should be initialized.");

  const char* env_template = getenv(env_naming_template);
  if (env_template) {
    if (strlen(env_template) > MAX_NAMING_TEMPLATE_SIZE) {
      fprintf(stderr, "Config error: %s is too long. ", env_naming_template);
      fprintf(stderr, "Max acceptable size: %d\n", MAX_NAMING_TEMPLATE_SIZE);
      return false;
    }
    strcpy(naming_template, env_template);
  } else {
    strcpy(naming_template, default_naming_template);
  }

  const char* evn_size = getenv(env_mmap_heap_size);
  if (evn_size) {
    mmap_heap_size = (size_t) strtoull(evn_size, NULL, 10);
    if (mmap_heap_size % page_size != 0) {
      fprintf(
        stderr, 
        "Config error: %s is not page aligned. Page size: %ld\n", 
        env_mmap_heap_size, 
        page_size
      );
      return false;
    }
  } else {
    mmap_heap_size = default_mmap_heap_size;
  }

  const char* evn_bsize = getenv(env_mmap_alloctor_min_bsize);
  if (evn_bsize) {
    mmap_alloctor_min_bsize = (size_t) strtoull(evn_bsize, NULL, 10);
    if (mmap_alloctor_min_bsize > mmap_heap_size) {
      fprintf(
        stderr, 
        "Config error: %s is not larger than the heap size.\n", 
        env_mmap_alloctor_min_bsize
      );
      return false;
    }
  } else {
    mmap_alloctor_min_bsize = default_mmap_alloctor_min_bsize;
  }

  const char* env_profile = getenv(env_profile_file_path);
  if (env_profile) {
    profile_file = fopen(env_profile, "a");
    if (!profile_file) {
      fprintf(
        stderr,
        "Failed tp open the profile file in path: %s\n",
        env_profile
      );
      return false;
    }
  }

  const char* env_profile_freq = getenv(env_profile_frequency);
  if (env_profile_freq) {
    profile_frequency = (size_t) strtoull(env_profile_freq, NULL, 10);
    if (profile_frequency == 0) {
      fprintf(
        stderr, 
        "Config error: %s cannot be zero.\n", 
        env_mmap_alloctor_min_bsize
      );
      return false;
    }
  }

  return true;
}

/*---------------------------------------------------------------------------*/
// Allocator API from mmap_heap

LOCAL_HELPER void mmap_allocator_init() {
  GLOBAL_LOCK_ACQUIRE();
  if (allocator_status != NOT_LOADED) {
    GLOBAL_LOCK_RELEASE();
    return;
  }

  // Initialize page size.
  page_size = sysconf(_SC_PAGE_SIZE);

  // Read environment variables
  if (!config_parameters()) {
    fprintf(stderr, "Failed to configure parameters from env.\n");
    allocator_status = FAILED_TO_LOAD;
    GLOBAL_LOCK_RELEASE();
    return;
  }

  // Configure profiling
  if (profile_file) {
    pthread_t profile_thread;
    if (pthread_create(&profile_thread, NULL, profile_thread_entry, NULL)) {
      fprintf(stderr, "Profiling thread init failed.\n");
      allocator_status = FAILED_TO_LOAD;
      GLOBAL_LOCK_RELEASE();
      return;
    }
  }

  // Reserve mmap region.
  mmap_region_base = mmap_reserve(mmap_heap_size);
  if (!mmap_region_base) {
    // Failed to reserve the mmap region.
    fprintf(stderr, "Failed to reserve mmap alloctor region.\n");
    allocator_status = FAILED_TO_LOAD;
    GLOBAL_LOCK_RELEASE();
    return;
  }
  mmap_region_end = mmap_region_base + mmap_heap_size - 1;

  // Initialize mmap heap.
  if (!heap_init(&mmap_heap, mmap_region_base, mmap_heap_size)) {
    allocator_status = FAILED_TO_LOAD;
    GLOBAL_LOCK_RELEASE();
    return;
  }

  fprintf(stderr, "MMap Allocator is successfully loaded.\n");
  allocator_status = LOADED;
  GLOBAL_LOCK_RELEASE();
}

LOCAL_HELPER void* mmap_allocate(size_t size) {
  if (size == 0) {
    // Ideally will never happen.
    return NULL;
  }

  // Round to ensure page aligned
  size = CEILING_PAGE_SIZE(size);

  GLOBAL_LOCK_ACQUIRE();
  list_node_t block = heap_allocate(&mmap_heap, size);
  GLOBAL_LOCK_RELEASE();

  if (!block) {
    return NULL;
  }

  assert(block->size == size);
  const int retval = mmap_maptemp(block->addr, block->size, naming_template);
  if (0 != retval) {
    // Failed to allocate the swap file.
    // We will return this node back to theheap.
    GLOBAL_LOCK_ACQUIRE();
    heap_free(&mmap_heap, block);
    GLOBAL_LOCK_RELEASE();
    fprintf(stderr, "Failed to create mmap region. Error code: %d\n", retval);
    return NULL;
  }

  return (void*) block->addr;
}

LOCAL_HELPER bool mmap_release(void* addr) {
  GLOBAL_LOCK_ACQUIRE();
  list_node_t block = list_find_in_use(&mmap_heap.node_list, addr);
  if (!block) {
    fprintf(stderr, "Failed to find the block with the given addr: %p\n", addr);
    GLOBAL_LOCK_RELEASE();
    return false;
  }

  // Unmap the region.
  const int retval = mmap_unmap(addr, block->size);
  if (0 != retval) {
    fprintf(stderr, "Failed to unmap the region.\n");
    GLOBAL_LOCK_RELEASE();
    return false;
  }

  if (!heap_free(&mmap_heap, block)) {
    fprintf(
      stderr, 
      "Failed to free as stdlib cannot allocate more memory for the heap.\n"
    );
    GLOBAL_LOCK_RELEASE();
    return false;
  }

  GLOBAL_LOCK_RELEASE();
  return true;
}

LOCAL_HELPER bool mmap_release_with_copy(
  list_node_t block, 
  void* dest, 
  const size_t new_size
) {
  void* addr = block->addr;
  // We will need to copy the data from the src to the new dest.
  const size_t copy_size = (block->size < new_size) ? block->size : new_size;
  memcpy(dest, addr, copy_size);

  // Unmap the region.
  if (mmap_unmap(addr, block->size)) {
    fprintf(stderr, "Failed to unmap the region.\n");
    GLOBAL_LOCK_RELEASE();
    return false;
  }

  GLOBAL_LOCK_ACQUIRE();
  bool success = heap_free(&mmap_heap, block);
  GLOBAL_LOCK_RELEASE();

  if (!success) {
    fprintf(
      stderr, 
      "Failed to free as stdlib cannot allocate more memory for the heap.\n"
    );
    return false;
  }
  return true;
}

/*---------------------------------------------------------------------------*/
// User library interface
void* mmap_malloc(size_t size);
void* mmap_calloc(size_t num_elements, size_t element_size);
void* mmap_reallocarray(void *addr, size_t size, size_t count);
void* mmap_realloc(void *addr, size_t size);

/*---------------------------------------------------------------------------*/
// malloc implementation
void* mmap_malloc(size_t size) {
  if (allocator_status == NOT_LOADED) {
    mmap_allocator_init();
  }

  if (allocator_status != LOADED || size < mmap_alloctor_min_bsize) {
    return malloc(size);
  }

  void* ret = mmap_allocate(size);
  if (!ret) {
    OUT_OF_MEMORY();
  }
  HEAP_CHECK();
  return ret;
}

/*---------------------------------------------------------------------------*/
// calloc implementation
void* mmap_calloc(size_t num_elements, size_t element_size) {
  if (allocator_status == NOT_LOADED) {
    mmap_allocator_init();
  }

  const size_t total_size = num_elements * element_size;
  if (allocator_status != LOADED || total_size < mmap_alloctor_min_bsize) {
    return calloc(num_elements, element_size);
  }

  void* ret = mmap_allocate(total_size);
  if (!ret) {
    OUT_OF_MEMORY();
  }
  // As required of calloc, memory region must be initialized to 0.
  memset(ret, 0, total_size);
  return ret;
}

/*---------------------------------------------------------------------------*/
// reallocarray implementation
void* mmap_reallocarray(void *addr, size_t size, size_t count) {
  return mmap_realloc(addr, size * count);
}

/*---------------------------------------------------------------------------*/
// realloc implementation
void* mmap_realloc(void *addr, size_t size) {
  if (!addr) return mmap_malloc(size); // It should behave like malloc.

  if (allocator_status == NOT_LOADED) {
    mmap_allocator_init();
  }

  if (allocator_status != LOADED || addr < mmap_region_base || addr > mmap_region_end) {
    if (size < mmap_alloctor_min_bsize) {
      return realloc(addr, size);
    } else {
      size_t old_size = MALLOCSIZE(addr);
	  fprintf(stderr, "Realloc Heap[%zu] to Mmap[%zu].\n",old_size,size);
      if (old_size > size) old_size = size;
      void* new_region = mmap_allocate(size);
      if (!new_region) {
        OUT_OF_MEMORY();
        return NULL;
      }
      memcpy(new_region, addr, old_size);
      free(addr);
      HEAP_CHECK();
      return new_region;
    }
  }

  // mmap allocator is properly initialized.
  if (addr >= mmap_region_base && addr <= mmap_region_end) {

    GLOBAL_LOCK_ACQUIRE();
    list_node_t block = list_find_in_use(&mmap_heap.node_list, addr);
    GLOBAL_LOCK_RELEASE();
    if (!block) {
      fprintf(stderr, 
        "Realloc failed: cannot find the given address: %p.\n", addr);
      return NULL;
    }
    
    if (block->size >= size) return addr; // Already large enough.
    
    // const size_t alloc_size = CEILING_PAGE_SIZE(size);
    void* new_region = mmap_allocate(size);
    if (!new_region) {
      OUT_OF_MEMORY();
      return NULL;
    }

    if (!mmap_release_with_copy(block, new_region, size)) {
      fprintf(stderr, "Realloc copy failed.\n");
      mmap_release(new_region);
      return NULL;
    }

    HEAP_CHECK();
    return new_region;
  }
  fprintf(stderr, "Failed in mmap_relloc.\n");
  return NULL;
}

/*---------------------------------------------------------------------------*/
// free implementation
void mmap_free(void* addr) {
  if (!addr) return;

  if (allocator_status == NOT_LOADED) {
    mmap_allocator_init();
  }

  if (allocator_status != LOADED || addr < mmap_region_base || addr > mmap_region_end) {
  /* On macOS, `addr` can be larger than `mmap_region_end` */
    free(addr);
    return;
  }

  if (!mmap_release(addr)) {
    fprintf(stderr, "Failed to free at addr: %p\n", addr);
  }
  HEAP_CHECK();
}