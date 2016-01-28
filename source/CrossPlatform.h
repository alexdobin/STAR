#ifdef _WIN32
#define sleep(x) Sleep(x)
typedef __int64 int64;
typedef int key_t;
typedef unsigned int mode_t;
#else
#include <inttypes.h>
typedef int64_t int64;
#endif
