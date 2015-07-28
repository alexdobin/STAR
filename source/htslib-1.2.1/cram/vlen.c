/*
Author: James Bonfield (jkb@sanger.ac.uk)

Copyright (c) 1995-1996 MEDICAL RESEARCH COUNCIL
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

/*
Copyright (c) 2004, 2009, 2011-2012 Genome Research Ltd.

Author: James Bonfield <jkb@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/types.h>
#include <string.h>

#include "cram/vlen.h"
#include "cram/os.h"

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef ABS
#define ABS(a) ((a)>0?(a):-(a))
#endif

/* #define DEBUG_printf(a,n) printf(a,n) */
#define DEBUG_printf(a,n)

/*
 * vlen: 27/10/95 written by James Bonfield, jkb@mrc-lmb.cam.ac.uk
 *
 * Given sprintf style of arguments this routine returns the maximum
 * size of buffer needed to allocate to use with sprintf. It errs on
 * the side of caution by being simplistic in its approach: we assume
 * all numbers are of maximum length.
 *
 * Handles the usual type conversions (%[%diuaxXcfeEgGpns]), but not
 * the 'wide' character conversions (%C and %S).
 * Precision is handled in the correct formats, including %*.*
 * notations.
 * Additionally, some of the more dubious (but probably illegal) cases
 * are supported (eg "%10%" will expand to "         %" on many
 * systems).
 *
 * We also assume that the largest integer and larger pointer are 64
 * bits, which at least covers the machines we'll need it for.
 */
int flen(char *fmt, ...)
{
    va_list args;

    va_start(args, fmt);
    return vflen(fmt, args);
}

int vflen(char *fmt, va_list ap)
{
    int len = 0;
    char *cp, c;
    long long l;
    int i;
    double d; 

    /*
     * This code modifies 'ap', but we do not know if va_list is a structure
     * or a pointer to an array so we do not know if it is a local variable
     * or not.
     * C99 gets around this by defining va_copy() to make copies of ap, but
     * this does not exist on all systems.
     * For now, I just assume that when va_list is a pointer the system also
     * provides a va_copy macro to work around this problem. The only system
     * I have seen needing this so far was Linux on AMD64.
     */
#if defined(HAVE_VA_COPY)
    va_list ap_local;
    va_copy(ap_local, ap);
#    define ap ap_local
#endif

    for(cp = fmt; *cp; cp++) {
	switch(*cp) {

	/* A format specifier */
	case '%': {
	    char *endp;
	    long conv_len1=0, conv_len2=0, conv_len=0;
	    signed int arg_size;

	    /* Firstly, strip the modifier flags (+-#0 and [space]) */
	    for(; (c=*++cp);) {
		if ('#' == c)
		    len+=2; /* Worst case of "0x" */
		else if ('-' == c || '+' == c || ' ' == c)
		    len++;
		else
		    break;
	    }

	    /* Width specifier */
	    l = strtol(cp, &endp, 10);
	    if (endp != cp) {
		cp = endp;
		conv_len = conv_len1 = l;
	    } else if (*cp == '*') {
		conv_len = conv_len1 = (int)va_arg(ap, int);
		cp++;
	    }

	    /* Precision specifier */
	    if ('.' == *cp) {
		cp++;
		conv_len2 = strtol(cp, &endp, 10);
		if (endp != cp) {
		    cp = endp;
		} else if (*cp == '*') {
		    conv_len2 = (int)va_arg(ap, int);
		    cp++;
		}
		conv_len = MAX(conv_len1, conv_len2);
	    }

	    /* Short/long identifier */
	    if ('h' == *cp) {
		arg_size = -1; /* short */
		cp++;
	    } else if ('l' == *cp) {
		arg_size = 1; /* long */
		cp++;
		if ('l' == *cp) {
		    arg_size = 2; /* long long */
		    cp++;
		}
	    } else {
		arg_size = 0; /* int */
	    }

	    /* The actual type */
	    switch (*cp) {
	    case '%':
		/*
		 * Not real ANSI I suspect, but we'll allow for the
		 * completely daft "%10%" example.
		 */
		len += MAX(conv_len1, 1);
		break;

	    case 'd':
	    case 'i':
	    case 'u':
	    case 'a':
	    case 'x':
	    case 'X':
		/* Remember: char and short are sent as int on the stack */
		if (arg_size == -1)
		    l = (long)va_arg(ap, int);
		else if (arg_size == 1)
		    l = va_arg(ap, long); 
		else if (arg_size == 2)
		    l = va_arg(ap, long long); 
		else 
		    l = (long)va_arg(ap, int);

		DEBUG_printf("%d", l);

		/*
		 * No number can be more than 24 characters so we'll take
		 * the max of conv_len and 24 (23 is len(2^64) in octal).
		 * All that work above and we then go and estimate ;-),
		 * but it's needed incase someone does %500d.
		 */
		len += MAX(conv_len, 23);
		break;

	    case 'c':
		i = va_arg(ap, int);
		DEBUG_printf("%c", i);
		/*
		 * Note that %10c and %.10c act differently.
		 * Besides, I think precision is not really allowed for %c.
		 */
		len += MAX(conv_len1, i>=0x80 ?MB_CUR_MAX :1);
		break;

	    case 'f':
		d = va_arg(ap, double);
		DEBUG_printf("%f", d);
		/*
		 * Maybe "Inf" or "NaN", but we'll not worry about that.
		 * Again, err on side of caution and take max of conv_len
		 * and max length of a double. The worst case I can
		 * think of is 317 characters (-1[308 zeros].000000)
		 * without using precision codes. That's horrid. I
		 * cheat and either use 317 or 15 depending on how
		 * large the number is as I reckon 99% of floats
		 * aren't that long.
		 */
		l = (ABS(d) > 1000000) ? 317 : 15;
		l = MAX(l, conv_len1 + 2);
		if (conv_len2) l += conv_len2 - 6;
		len += l;
		break;

	    case 'e':
	    case 'E':
	    case 'g':
	    case 'G':
		d = va_arg(ap, double);
		DEBUG_printf("%g", d);
		/*
		 * Maybe "Inf" or "NaN", but we'll not worry about that
		 * Again, err on side of caution and take max of conv_len
		 * and max length of a double (which defaults to only
		 * '-' + 6 + '.' + 'E[+-]xxx' == 13.
		 */
		len += MAX(conv_len, 13);
		break;

	    case 'p':
		l = (long)va_arg(ap, void *);
		/*
		 * Max pointer is 64bits == 16 chars (on alpha),
		 * == 20 with + "0x".
		 */
		DEBUG_printf("%p", (void *)l);
		len += MAX(conv_len, 20);
		break;

	    case 'n':
		/* produces no output */
		break;

	    case 's': {
		char *s = (char *)va_arg(ap, char *);
		DEBUG_printf("%s", s);

		if (!conv_len2) {
		    len += MAX(conv_len, (int)strlen(s));
		} else {
		    len += conv_len;
		}
		break;
	    }

	    default:
		/* wchar_t types of 'C' and 'S' aren't supported */
		DEBUG_printf("Arg is %c\n", *cp);
	    }
	    
	}

	case '\0':
	    break;

	default:
	    DEBUG_printf("%c", *cp);
	    len++;
	}
    }

    va_end(ap);

    return len+1; /* one for the null character */
}

#if 0
int main() {
    int l;
    char buf[10000];

    sprintf(buf, "d: %d\n", 500);
    l = flen("d: %d\n", 500);
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "");
    l = flen("");
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%s\n","test");
    l = flen("%s\n", "test");
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%c\n", 'a');
    l = flen("%c\n", 'a');
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%31.30f\n", -9999.99);
    l = flen("%31.30f\n", -9999.99);
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%f\n", -1e308);
    l = flen("%f\n", -1e308);
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%.9f\n", -1e308);
    l = flen("%.9f\n", -1e308);
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%10.20f\n", -1.999222333);
    l = flen("%10.20f\n", -1.999222333);
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%#g\n", -3.14159265358e-222);
    l = flen("%#g\n", -3.1415927e-222);
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%e\n", -123456789123456789.1);
    l = flen("%e\n", -123456789123456789.1);
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%c %f %d %s %c %g %ld %s\n", 'a', 3.1, 9, "one", 'b', 4.2, 9, "two");
    l = flen("%c %f %d %s %c %g %ld %s\n", 'a', 3.1, 9, "one", 'b', 4.2, 9, "two");
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%*.*e %*c\n", 10, 5, 9.0, 20, 'x');
    l = flen("%*.*e %*c\n", 10, 5, 9.0, 20, 'x');
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%10c\n", 'z');
    l = flen("%10c\n", 'z');
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%.10c\n", 'z');
    l = flen("%.10c\n", 'z');
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%10d\n", 'z');
    l = flen("%10d\n", 'z');
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%.10d\n", 'z');
    l = flen("%.10d\n", 'z');
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%10%\n");
    l = flen("%10%\n");
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%.10%\n");
    l = flen("%.10%\n");
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%s\n", "0123456789");
    l = flen("%s\n", "0123456789");
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%5s\n", "0123456789");
    l = flen("%5s\n", "0123456789");
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%50s\n", "0123456789");
    l = flen("%50s\n", "0123456789");
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%.5s\n", "0123456789");
    l = flen("%.5s\n", "0123456789");
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%.50s\n", "0123456789");
    l = flen("%.50s\n", "0123456789");
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%5.50s\n", "0123456789");
    l = flen("%5.50s\n", "0123456789");
    printf("%d %d\n\n", strlen(buf), l);

    sprintf(buf, "%50.5s\n", "0123456789");
    l = flen("%50.5s\n", "0123456789");
    printf("%d %d\n\n", strlen(buf), l);

    return 0;
}
#endif
