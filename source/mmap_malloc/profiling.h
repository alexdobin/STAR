#ifndef PROFILING_H
#define PROFILING_H

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#include "constants.h"
#include "default_config.h"

extern pthread_mutex_t profile_lock;
extern FILE* profile_file;
extern size_t profile_frequency;

extern size_t num_mmap_file;
extern size_t mmap_heap_total_size;

#define PROFILE_LOCK_ACQUIRE() do {\
  pthread_mutex_lock(&profile_lock);\
} while (0)

#define PROFILE_LOCK_RELEASE() do {\
  pthread_mutex_unlock(&profile_lock);\
} while (0)

void FORCE_INLINE
profile_allocate(const size_t size) {
  if (!profile_file) return;
  PROFILE_LOCK_ACQUIRE();
  ++num_mmap_file;
  mmap_heap_total_size += size;
  PROFILE_LOCK_RELEASE();
}

void FORCE_INLINE
profile_deallocate(const size_t size) {
  if (!profile_file) return;
  PROFILE_LOCK_ACQUIRE();
  --num_mmap_file;
  mmap_heap_total_size -= size;
  PROFILE_LOCK_RELEASE();
}

void* profile_thread_entry(void* arg);

#endif