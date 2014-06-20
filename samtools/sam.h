#ifndef BAM_SAM_H
#define BAM_SAM_H

#include "bam.h"

/*!
  @header

  This file provides higher level of I/O routines and unifies the APIs
  for SAM and BAM formats. These APIs are more convenient and
  recommended.

  @copyright Genome Research Ltd.
 */

/*! @typedef
  @abstract SAM/BAM file handler
  @field  type    type of the handler; bit 1 for BAM, 2 for reading and bit 3-4 for flag format
  @field  bam   BAM file handler; valid if (type&1) == 1
  @field  tamr  SAM file handler for reading; valid if type == 2
  @field  tamw  SAM file handler for writing; valid if type == 0
  @field  header  header struct
 */
typedef struct {
	int type;
	union {
		tamFile tamr;
		bamFile bam;
		FILE *tamw;
	} x;
	bam_header_t *header;
} samfile_t;

#ifdef __cplusplus
extern "C" {
#endif

	/*!
	  @abstract     Open a SAM/BAM file

	  @param fn SAM/BAM file name; "-" is recognized as stdin (for
	  reading) or stdout (for writing).

	  @param mode open mode /[rw](b?)(u?)(h?)([xX]?)/: 'r' for reading,
	  'w' for writing, 'b' for BAM I/O, 'u' for uncompressed BAM output,
	  'h' for outputing header in SAM, 'x' for HEX flag and 'X' for
	  string flag. If 'b' present, it must immediately follow 'r' or
	  'w'. Valid modes are "r", "w", "wh", "wx", "whx", "wX", "whX",
	  "rb", "wb" and "wbu" exclusively.

	  @param aux auxiliary data; if mode[0]=='w', aux points to
	  bam_header_t; if strcmp(mode, "rb")!=0 and @SQ header lines in SAM
	  are absent, aux points the file name of the list of the reference;
	  aux is not used otherwise. If @SQ header lines are present in SAM,
	  aux is not used, either.

	  @return       SAM/BAM file handler
	 */
	samfile_t *samopen(const char *fn, const char *mode, const void *aux);

	/*!
	  @abstract     Close a SAM/BAM handler
	  @param  fp    file handler to be closed
	 */
	void samclose(samfile_t *fp);

	/*!
	  @abstract     Read one alignment
	  @param  fp    file handler
	  @param  b     alignment
	  @return       bytes read
	 */
	int samread(samfile_t *fp, bam1_t *b);

	/*!
	  @abstract     Write one alignment
	  @param  fp    file handler
	  @param  b     alignment
	  @return       bytes written
	 */
	int samwrite(samfile_t *fp, const bam1_t *b);

	/*!
	  @abstract     Get the pileup for a whole alignment file
	  @param  fp    file handler
	  @param  mask  mask transferred to bam_plbuf_set_mask()
	  @param  func  user defined function called in the pileup process
	  #param  data  user provided data for func()
	 */
	int sampileup(samfile_t *fp, int mask, bam_pileup_f func, void *data);

	char *samfaipath(const char *fn_ref);
	int samthreads(samfile_t *fp, int n_threads, int n_sub_blks);

#ifdef __cplusplus
}
#endif

#endif
