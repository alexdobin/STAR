#define _GNU_SOURCE // To enable various non-standard GNU extensions.
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>

#include "constants.h"
#include "mmap_mgr.h"
#include "profiling.h"
#include <stdlib.h>

/*
This function will reserve a region for later mmap use.
It does not map to an existing file.
*/
void* mmap_reserve(const size_t size) {
  void* addr = mmap(NULL, size, PROT_NONE, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0);
  if (addr == MAP_FAILED) {
    return NULL;
  }
  return addr;
}

int mmap_maptemp(void* addr, const size_t size, char* template) {
  char* filename = (char*) malloc(strlen(template) + 1);
  if (!filename) {
    return -1;
  }
  strcpy(filename, template);
  int fd = mkstemp(filename);
  if (fd < 0) {
    free(filename);
    return -2;
  }

  // After unlinking the file it will not appear in the file system,
  // however, since we have the fd, it can still be written.
  int retval;
  retval = unlink(filename);
  free(filename);
  if (retval) {
    return -3;
  }
  // Resize the file to ensure we have enough space.
  retval = ftruncate(fd, size);
  if (retval) {
    fprintf(stderr, "ftruncate failed! %s. Size: %ld\n", strerror(errno), size);
    return -4;
  }

  void* ret_addr = mmap(
    addr, // Map to this address
    size, // With the given size.
    PROT_READ | PROT_WRITE,
    MAP_SHARED | MAP_FIXED,
    fd,
    0
  );
  if (ret_addr == MAP_FAILED) {
    return -5;
  }

  profile_allocate(size);
  
  retval = close(fd);
  if (retval) {
    return -6;
  }

  return 0;
}

int mmap_unmap(void* addr, const size_t size) {
  void* ret_addr = mmap(
    addr, // Previous mapped address
    size, // With the given size.
    PROT_NONE,
    MAP_ANONYMOUS | MAP_FIXED | MAP_PRIVATE,
    -1, // Reserved
    0
  );
  if (ret_addr == MAP_FAILED) {
    return -1;
  }

  profile_deallocate(size);
  return 0;
}