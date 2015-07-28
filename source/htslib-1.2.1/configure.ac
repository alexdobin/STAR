# Configure script for htslib, a C library for high-throughput sequencing data.
#
#    Copyright (C) 2015 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

dnl Process this file with autoconf to produce a configure script
AC_INIT([HTSlib], m4_esyscmd_s([make print-version]),
        [samtools-help@lists.sourceforge.net], [], [http://www.htslib.org/])
AC_PREREQ(2.63)  dnl This version introduced 4-argument AC_CHECK_HEADER
AC_CONFIG_SRCDIR(hts.c)

dnl Copyright notice to be copied into the generated configure script
AC_COPYRIGHT([Portions copyright (C) 2015 Genome Research Ltd.

This configure script is free software: you are free to change and
redistribute it.  There is NO WARRANTY, to the extent permitted by law.])

AC_PROG_CC
AC_PROG_RANLIB

AC_ARG_WITH([irods],
  [AS_HELP_STRING([[--with-irods[=DIR]]],
                  [use RodsAPIs library (in DIR) to support iRODS URLs])],
  [case $withval in
     no)  irods=disabled ;;
     yes) irods=enabled ;;
     *)   irods=enabled; IRODS_HOME=$withval ;;
   esac],
  [irods=disabled])

save_LIBS=$LIBS
zlib_devel=ok
dnl Set a trivial non-empty INCLUDES to avoid excess default includes tests
AC_CHECK_HEADER([zlib.h], [], [zlib_devel=missing], [;])
AC_CHECK_LIB(z, inflate,  [], [zlib_devel=missing])
LIBS=$save_LIBS

if test $zlib_devel != ok; then
  AC_MSG_ERROR([zlib development files not found

HTSlib uses compression routines from the zlib library <http://zlib.net>.
Building HTSlib requires zlib development files to be installed on the build
machine; you may need to ensure a package such as zlib1g-dev (on Debian or
Ubuntu Linux) or zlib-devel (on RPM-based Linux distributions) is installed.

FAILED.  This error must be resolved in order to build HTSlib successfully.])
fi

if test $irods = enabled; then
  # TODO Also test whether we require libgssapi_krb5 and AC_CHECK_LIB it
  save_LDFLAGS=$LDFLAGS
  LDFLAGS="$LDFLAGS -L$IRODS_HOME/lib/core/obj"
  AC_CHECK_LIB([RodsAPIs], [getRodsEnvFileName],
    [case $with_irods in
       yes) define_IRODS_HOME='# Uses $(IRODS_HOME) from the environment' ;;
       *)   define_IRODS_HOME="IRODS_HOME = $with_irods" ;;
     esac],
    [AC_MSG_ERROR([iRODS development files not found

Support for iRODS URLs requires the libRodsAPI client library and headers.
Configure with --with-irods=DIR (or just --with-irods if \$IRODS_HOME has
been exported with a suitable value), where DIR is the base of an iRODS tree
such that the library is present as DIR/lib/core/obj/libRodsAPI.* and headers
are present under DIR/lib/api/include and so on.])],
    [-lgssapi_krb5 -lpthread])
  LDFLAGS=$save_LDFLAGS
else
  define_IRODS_HOME='IRODS_HOME ?= /disabled'
fi
AC_SUBST([irods])
AC_SUBST([define_IRODS_HOME])

AC_CONFIG_FILES(config.mk)
AC_OUTPUT
