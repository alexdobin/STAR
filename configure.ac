#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([STAR], [2.4.2a], [https://github.com/alexdobin/STAR/issues])
AC_CONFIG_SRCDIR([source/STAR.cpp])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([foreign])

AC_LANG(C++)

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_OPENMP

# Checks for libraries.
PKG_CHECK_MODULES(HTSLIB, htslib)

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h limits.h stdint.h stdlib.h string.h sys/statvfs.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT8_T
AC_TYPE_MODE_T
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_FUNC_REALLOC
AC_FUNC_FORK
AC_CHECK_FUNCS([floor ftruncate memmove memset mkdir mkfifo munmap rmdir strchr strerror strstr])

AC_CONFIG_FILES([Makefile
                 source/Makefile])

AC_OUTPUT
