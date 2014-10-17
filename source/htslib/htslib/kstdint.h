#ifndef KSTDINT_H
#define KSTDINT_H

#include <limits.h>

/* Basic assumptions: 1) "char" is 8-bit; 2) there is a 8-bit, 16-bit, 32-bit
 * and 64-bit integer type, respectively; 3) "short" is no less than "char",
 * "int" is no less than "short", "long" is no less than "int" and "long long"
 * is no less than "long"; 4) "int" is at least 16-bit, "long" at least 32-bit
 * and "long long" at least 64-bit. The last two assumptions are enforced by
 * the C99 spec.
 *
 * Following assumptions 1) and 2), we know that "signed char"=="int8_t" and
 * "short"=="int16_t" for sure. Further from the assumptions, a 32-bit integer
 * type must be either "int" or "long". We can test (UINT16_MAX==UINT_MAX) to
 * see which is the case. Similarly, a 64-bit integer must be either "long" or
 * "long long". We can test (UINT16_MAX==UINT_MAX) to get the definite answer.
 */

/* 8-bit integers */
typedef signed char          int8_t;
typedef unsigned char       uint8_t;
#define INT8_MIN     (-SCHAR_MAX-1)
#define INT8_MAX          SCHAR_MAX
#define UINT8_MAX         UCHAR_MAX

/* 16-bit integers */
typedef signed short        int16_t;
typedef unsigned short     uint16_t;
#define INT16_MIN     (-SHRT_MAX-1)
#define INT16_MAX          SHRT_MAX
#define UINT16_MAX        USHRT_MAX

/* 32-bit integers */
#if UINT16_MAX != UINT_MAX
typedef signed int          int32_t;
typedef unsigned int       uint32_t;
#define INT32_MIN      (-INT_MAX-1)
#define INT32_MAX           INT_MAX
#define UINT32_MAX         UINT_MAX
#else /* then int is 16-bit and long is 32-bit, which may happen to compilers for embedded CPUs */
typedef signed long         int32_t;
typedef unsigned long      uint32_t;
#define INT32_MIN     (-LONG_MAX-1)
#define INT32_MAX          LONG_MAX
#define UINT32_MAX        ULONG_MAX
#endif /* ~UINT16_MAX!=UINT_MAX */

/* 64-bit integers */
#if UINT32_MAX != ULONG_MAX
typedef signed long         int64_t;
typedef unsigned long      uint64_t;
#define INT64_MIN     (-LONG_MAX-1)
#define INT64_MAX          LONG_MAX
#define UINT64_MAX        ULONG_MAX
#else
typedef signed long long    int64_t;
typedef unsigned long long uint64_t;
#define INT64_MIN    (-LLONG_MAX-1)
#define INT64_MAX         LLONG_MAX
#define UINT64_MAX       ULLONG_MAX
#endif /* ~UINT32_MAX!=ULONG_MAX */

#endif /* ~defined(KSTDINT_H) */
