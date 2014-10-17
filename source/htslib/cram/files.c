/*
Copyright (c) 1994, 1996-1997, 2000, 2003 MEDICAL RESEARCH COUNCIL
All rights reserved

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

   1 Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

   2 Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

   3 Neither the name of the MEDICAL RESEARCH COUNCIL, THE LABORATORY OF 
MOLECULAR BIOLOGY nor the names of its contributors may be used to endorse or 
promote products derived from this software without specific prior written 
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR 
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include "cram/misc.h"

#include <sys/types.h>
#include <sys/stat.h>
/* Alliant's Concentrix <sys/stat.h> is hugely deficient */
/* Define things we require in this program              */
/* Methinks S_IFMT and S_IFDIR aren't defined in POSIX   */
#ifndef S_ISDIR
#define S_ISDIR(m)      (((m)&S_IFMT) == S_IFDIR)
#endif /*!S_ISDIR*/
#ifndef S_ISREG
#define S_ISREG(m)      (((m)&S_IFMT) == S_IFREG)
#endif /*!S_ISREG*/

int is_directory(char * fn)
{
    struct stat buf;
    if ( stat(fn,&buf) ) return 0;
    return S_ISDIR(buf.st_mode);
}

int is_file(char * fn)
{
    struct stat buf;
    if ( stat(fn,&buf) ) return 0;
    return S_ISREG(buf.st_mode);
}

int file_exists(char * fn)
{
    struct stat buf;
    return ( stat(fn,&buf) == 0);
}

int file_size(char * fn)
{
    struct stat buf;
    if ( stat(fn,&buf) != 0) return 0;
    return buf.st_size;
}

