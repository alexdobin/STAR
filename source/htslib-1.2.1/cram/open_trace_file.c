/*
Author: James Bonfield

Copyright (c) 2000-2001 MEDICAL RESEARCH COUNCIL
All rights reserved

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

   3. Neither the name of the MEDICAL RESEARCH COUNCIL, THE LABORATORY OF 
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
Copyright (c) 2008, 2009, 2013, 2014 Genome Research Ltd.
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "cram/os.h"
#ifndef PATH_MAX
#  define PATH_MAX 1024
#endif

#include "cram/open_trace_file.h"
#include "cram/misc.h"
#include "htslib/hfile.h"

/*
 * Tokenises the search path splitting on colons (unix) or semicolons
 * (windows).
 * We also  explicitly add a "./" to the end of the search path
 *
 * Returns: A new search path with items separated by nul chars. Two nul
 *          chars in a row represent the end of the tokenised path.
 * Returns NULL for a failure.
 *
 * The returned data has been malloced. It is up to the caller to free this
 * memory.
 */
char *tokenise_search_path(char *searchpath) {
    char *newsearch;
    unsigned int i, j;
    size_t len;
#ifdef _WIN32
    char path_sep = ';';
#else
    char path_sep = ':';
#endif

    if (!searchpath)
	searchpath="";

    newsearch = (char *)malloc((len = strlen(searchpath))+5);
    if (!newsearch)
	return NULL;

    for (i = 0, j = 0; i < len; i++) {
	/* "::" => ":". Used for escaping colons in http://foo */
	if (i < len-1 && searchpath[i] == ':' && searchpath[i+1] == ':') {
	    newsearch[j++] = ':';
	    i++;
	    continue;
	}

	/* Handle http:// and ftp:// too without :: */
	if (path_sep == ':') {
	    if ((i == 0 || (i > 0 && searchpath[i-1] == ':')) &&
		(!strncmp(&searchpath[i], "http:",     5) ||
		 !strncmp(&searchpath[i], "ftp:",      4) ||
		 !strncmp(&searchpath[i], "|http:",    6) ||
		 !strncmp(&searchpath[i], "|ftp:",     5) ||
		 !strncmp(&searchpath[i], "URL=http:", 9) ||
		 !strncmp(&searchpath[i], "URL=ftp:",  8))) {
		do {
		    newsearch[j++] = searchpath[i];
		} while (i<len && searchpath[i++] != ':');
		if (searchpath[i] == ':')
		    i++;
		if (searchpath[i]=='/')
		    newsearch[j++] = searchpath[i++];
		if (searchpath[i]=='/')
		    newsearch[j++] = searchpath[i++];
		// Look for host:port
		do {
		    newsearch[j++] = searchpath[i++];
		} while (i<len && searchpath[i] != ':' && searchpath[i] != '/');
		newsearch[j++] = searchpath[i++];
		if (searchpath[i] == ':')
		    i++;
	    }
	}

	if (searchpath[i] == path_sep) {
	    /* Skip blank path components */
	    if (j && newsearch[j-1] != 0)
		newsearch[j++] = 0;
	} else {
	    newsearch[j++] = searchpath[i];
	}
    }

    if (j)
	newsearch[j++] = 0;
    newsearch[j++] = '.';
    newsearch[j++] = '/';
    newsearch[j++] = 0;
    newsearch[j++] = 0;
    
    return newsearch;
}

mFILE *find_file_url(char *file, char *url) {
    char buf[8192], *cp;
    mFILE *mf = NULL;
    int maxlen = 8190 - strlen(file), len;
    hFILE *hf;

    /* Expand %s for the trace name */
    for (cp = buf; *url && cp - buf < maxlen; url++) {
	if (*url == '%' && *(url+1) == 's') {
	    url++;
	    cp += strlen(strcpy(cp, file));
	} else {
	    *cp++ = *url;
	}
    }
    *cp++ = 0;

    if (!(hf = hopen(buf, "r")))
	return NULL;

    if (NULL == (mf = mfcreate(NULL, 0)))
	return NULL;
    while ((len = hread(hf, buf, 8192)) > 0) {
	if (mfwrite(buf, len, 1, mf) <= 0) {
	    hclose_abruptly(hf);
	    mfdestroy(mf);
	    return NULL;
	}
    }
    if (hclose(hf) < 0) {
	mfdestroy(mf);
	return NULL;
    }

    mrewind(mf);
    return mf;
}

/*
 * Searches for file in the directory 'dirname'. If it finds it, it opens
 * it. This also searches for compressed versions of the file in dirname
 * too.
 *
 * Returns mFILE pointer if found
 *         NULL if not
 */
static mFILE *find_file_dir(char *file, char *dirname) {
    char path[PATH_MAX+1];
    size_t len = strlen(dirname);
    char *cp;

    if (dirname[len-1] == '/')
	len--;

    /* Special case for "./" or absolute filenames */
    if (*file == '/' || (len==1 && *dirname == '.')) {
	sprintf(path, "%s", file);
    } else {
	/* Handle %[0-9]*s expansions, if required */
	char *path_end = path;
	*path = 0;
	while ((cp = strchr(dirname, '%'))) {
	    char *endp;
	    long l = strtol(cp+1, &endp, 10);
	    if (*endp != 's') {
		strncpy(path_end, dirname, (endp+1)-dirname);
		path_end += (endp+1)-dirname;
		dirname = endp+1;
		continue;
	    }
	    
	    strncpy(path_end, dirname, cp-dirname);
	    path_end += cp-dirname;
	    if (l) {
		strncpy(path_end, file, l);
		path_end += MIN(strlen(file), l);
		file     += MIN(strlen(file), l);
	    } else {
		strcpy(path_end, file);
		path_end += strlen(file);
		file     += strlen(file);
	    }
	    len -= (endp+1) - dirname;
	    dirname = endp+1;
	}
	strncpy(path_end, dirname, len);
	path_end += MIN(strlen(dirname), len);
	*path_end = 0;
	if (*file) {
	    *path_end++ = '/';
	    strcpy(path_end, file);
	}

	//fprintf(stderr, "*PATH=\"%s\"\n", path);
    }

    if (is_file(path)) {
	return mfopen(path, "rb");
    }

    return NULL;
}

/*
 * ------------------------------------------------------------------------
 * Public functions below.
 */

/*
 * Opens a trace file named 'file'. This is initially looked for as a
 * pathname relative to a file named "relative_to". This may (for
 * example) be the name of an experiment file referencing the trace
 * file. In this case by passing relative_to as the experiment file
 * filename the trace file will be picked up in the same directory as
 * the experiment file. Relative_to may be supplied as NULL.
 *
 * 'file' is looked for at relative_to, then the current directory, and then
 * all of the locations listed in 'path' (which is a colon separated list).
 * If 'path' is NULL it uses the RAWDATA environment variable instead.
 *
 * Returns a mFILE pointer when found.
 *           NULL otherwise.
 */
mFILE *open_path_mfile(char *file, char *path, char *relative_to) {
    char *newsearch;
    char *ele;
    mFILE *fp;

    /* Use path first */
    if (!path)
	path = getenv("RAWDATA");
    if (NULL == (newsearch = tokenise_search_path(path)))
	return NULL;
    
    /*
     * Step through the search path testing out each component.
     * We now look through each path element treating some prefixes as
     * special, otherwise we treat the element as a directory.
     */
    for (ele = newsearch; *ele; ele += strlen(ele)+1) {
	char *ele2;

	/*
	 * '|' prefixing a path component indicates that we do not
	 * wish to perform the compression extension searching in that
	 * location.
	 *
	 * NB: this has been removed from the htslib implementation.
	 */
	if (*ele == '|') {
	    ele2 = ele+1;
	} else {
	    ele2 = ele;
	}

	if (0 == strncmp(ele2, "URL=", 4)) {
	    if ((fp = find_file_url(file, ele2+4))) {
		free(newsearch);
		return fp;
	    }
	} else if (!strncmp(ele2, "http:", 5) ||
		   !strncmp(ele2, "ftp:", 4)) {
	    if ((fp = find_file_url(file, ele2))) {
		free(newsearch);
		return fp;
	    }
	} else if ((fp = find_file_dir(file, ele2))) {
	    free(newsearch);
	    return fp;
	} 
   }

    free(newsearch);

    /* Look in the same location as the incoming 'relative_to' filename */
    if (relative_to) {
	char *cp;
	char relative_path[PATH_MAX+1];
	strcpy(relative_path, relative_to);
	if ((cp = strrchr(relative_path, '/')))
	    *cp = 0;
	if ((fp = find_file_dir(file, relative_path)))
	    return fp;
    }

    return NULL;
}
