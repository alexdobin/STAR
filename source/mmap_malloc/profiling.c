#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include "profiling.h"

pthread_mutex_t profile_lock = PTHREAD_MUTEX_INITIALIZER;
FILE* profile_file = NULL;
size_t profile_frequency = 30;

size_t num_mmap_file = 0;
size_t mmap_heap_total_size = 0;

void* profile_thread_entry(void* arg) {
  (void) arg;
  assert(profile_file);

  while (1) {
    PROFILE_LOCK_ACQUIRE();
    time_t curr = time(NULL);
    size_t curr_heap_size = mmap_heap_total_size / (1024 * 1024); // In MBytes
    size_t curr_num_block = num_mmap_file;
    PROFILE_LOCK_RELEASE();

    fprintf(profile_file, "<timestamp> %ld\n", curr);
    fprintf(profile_file, "<# mmap blocks> %ld\n", curr_num_block);
    fprintf(
      profile_file, 
      "<current size of mmap heap> %ldMBytes\n\n", 
      curr_heap_size
    );
    fflush(profile_file);

    sleep(profile_frequency);
  }
}
