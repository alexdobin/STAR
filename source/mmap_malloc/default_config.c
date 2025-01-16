#include "default_config.h"

char* const env_mmap_heap_size = "ENV_MMAP_HEAP_SIZE";
char* const env_mmap_alloctor_min_bsize = "ENV_MMAP_ALLOCATOR_MIN_BSIZE";
char* const env_naming_template = "ENV_NAMING_TEMPLATE";

char* const env_profile_file_path = "ENV_PROFILE_FILE_PATH";
char* const env_profile_frequency = "ENV_PROFILE_FREQUENCY";

const size_t default_mmap_heap_size = (size_t) 1024 * 1024 * 1024 * 1024 * 4;
/*
$ grep -m1 "address sizes" /proc/cpuinfo
address sizes	: 46 bits physical, 57 bits virtual
*/
const size_t default_mmap_alloctor_min_bsize = (size_t) 1024 * 1024 * 1024 * 2;
char* const  default_naming_template = ".mmap_alloc.XXXXXXXXXX";

/*
48 bit is 256T, 44bit:4T, 16-1 possible.

https://www.linuxquestions.org/questions/linux-hardware-18/address-sizes-in-cpuinfo-757456/

Every x86 and x86_64 CPU has an address translation system that translates virtual addresses into physical addresses.

The virtual addresses are the ones programs actually use (they are the numerical values of pointer variables in the program).

The physical addresses are the ones the CPU sends to the memory controller to generate the signals needed to select specific memory locations.

Each process has its own translation tables, so a pointer in one process with identical value to a pointer in another process doesn't normally point to the same physical memory.

The older x86 architecture has a translation design that translates a 32 bit virtual address to a 32 bit physical address.

Newer x86 architecture also supports a translation (called PAE) that translates 32 bit virtual to 36 bit physical.

The x86_64 architecture has 48 bit virtual addresses. Pointers take 64 bits, but the top 17 bits of any 64 bit pointer are required to be identical (all 17 zeros or all 17 ones) so 16 of those 17 bits are redundant.

The x86_64 address translation takes in 48 bit virtual addresses on all models. But the number of physical addresses generated varies by model of CPU chip.

The architecture has an upper limit on the number of physical address bits. I think that is 52, but I'm not sure I remember correctly. But since no one will be using that much physical ram before current CPU models are obsolete, CPU chips don't actually support that much. They support various numbers of physical address bits from 36 up. 40, which you have, is pretty typical.
*/
