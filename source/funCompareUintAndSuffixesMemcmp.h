#ifndef CODE_funCompareUintAndSuffixesMemcmp
#define CODE_funCompareUintAndSuffixesMemcmp

#include <stdint.h>

extern char* g_funCompareUintAndSuffixesMemcmp_G;
extern uint64_t g_funCompareUintAndSuffixesMemcmp_N;
// Maximum number of bases to compare; UINT64_MAX means compare to the sentinel.
extern uint64_t g_funCompareUintAndSuffixesMemcmp_L;
int funCompareUintAndSuffixesMemcmp ( const void *a, const void *b);

#endif
