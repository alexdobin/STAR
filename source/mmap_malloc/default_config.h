#ifndef DEFAULT_CONFIG_H
#define DEFAULT_CONFIG_H

#include <stdlib.h>

#define MAX_NAMING_TEMPLATE_SIZE 256

extern char* const env_mmap_heap_size;
extern char* const env_mmap_alloctor_min_bsize;
extern char* const env_naming_template;

extern char* const env_profile_file_path;
extern char* const env_profile_frequency;

extern const size_t default_mmap_heap_size;
extern const size_t default_mmap_alloctor_min_bsize;
extern char* const  default_naming_template;

#endif