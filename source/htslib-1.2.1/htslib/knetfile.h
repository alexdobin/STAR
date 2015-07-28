/* The MIT License

   Copyright (c) 2008 by Genome Research Ltd (GRL).
                 2010 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef KNETFILE_H
#define KNETFILE_H

#include <stdint.h>
#include <fcntl.h>

#ifndef _WIN32
#define netread(fd, ptr, len) read(fd, ptr, len)
#define netwrite(fd, ptr, len) write(fd, ptr, len)
#define netclose(fd) close(fd)
#else
#include <winsock2.h>
#define netread(fd, ptr, len) recv(fd, ptr, len, 0)
#define netwrite(fd, ptr, len) send(fd, ptr, len, 0)
#define netclose(fd) closesocket(fd)
#endif

// FIXME: currently I/O is unbuffered

#define KNF_TYPE_LOCAL 1
#define KNF_TYPE_FTP   2
#define KNF_TYPE_HTTP  3

typedef struct knetFile_s {
	int type, fd;
	int64_t offset;
	char *host, *port;

	// the following are for FTP only
	int ctrl_fd, pasv_ip[4], pasv_port, max_response, no_reconnect, is_ready;
	char *response, *retr, *size_cmd;
	int64_t seek_offset; // for lazy seek
    int64_t file_size;

	// the following are for HTTP only
	char *path, *http_host;
} knetFile;

#define knet_tell(fp) ((fp)->offset)
#define knet_fileno(fp) ((fp)->fd)

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
	int knet_win32_init();
	void knet_win32_destroy();
#endif

	knetFile *knet_open(const char *fn, const char *mode);

	/* 
	   This only works with local files.
	 */
	knetFile *knet_dopen(int fd, const char *mode);

	/*
	  If ->is_ready==0, this routine updates ->fd; otherwise, it simply
	  reads from ->fd.
	 */
	ssize_t knet_read(knetFile *fp, void *buf, size_t len);

	/*
	  This routine only sets ->offset and ->is_ready=0. It does not
	  communicate with the FTP server.
	 */
	off_t knet_seek(knetFile *fp, off_t off, int whence);
	int knet_close(knetFile *fp);

#ifdef __cplusplus
}
#endif

#endif
