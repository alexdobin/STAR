# Makefile rules useful for third-party code using htslib's public API.
#
#    Copyright (C) 2013-2014 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>

# The makefile fragment included below provides variables that can be used
# to express dependencies on headers supplied by an in-development htslib.
# If your source file foo.c #includes <htslib/hts.h> and <htslib/kstring.h>,
# you can write the correct prerequisites for foo.o as:
#
#	HTSDIR = <path to htslib top-level directory>
#	include $(HTSDIR)/htslib.mk
#
#	foo.o: foo.c $(htslib_hts_h) $(HTSDIR)/htslib/kstring.h
#
# Variables are not provided for k*.h, as those never include other headers.

HTSPREFIX = $(HTSDIR)/
include $(HTSDIR)/htslib_vars.mk

# Rules for rebuilding an in-development htslib's static and shared libraries.
# If your program foo links with libhts, adding the appropriate prerequisite
# will cause the library to be rebuilt as necessary:
#
#	foo: foo.o $(HTSDIR)/libhts.a
#
# or similarly if your target requires any of the tools supplied:
#
#	bar.bed.bgz.tbi: bar.bed.bgz $(HTSDIR)/tabix
#		$(HTSDIR)/tabix -p bed bar.bed.bgz

HTSLIB_PUBLIC_HEADERS = \
	$(HTSDIR)/htslib/bgzf.h \
	$(HTSDIR)/htslib/faidx.h \
	$(HTSDIR)/htslib/hfile.h \
	$(HTSDIR)/htslib/hts.h \
	$(HTSDIR)/htslib/hts_defs.h \
	$(HTSDIR)/htslib/khash.h \
	$(HTSDIR)/htslib/klist.h \
	$(HTSDIR)/htslib/knetfile.h \
	$(HTSDIR)/htslib/kseq.h \
	$(HTSDIR)/htslib/ksort.h \
	$(HTSDIR)/htslib/kstdint.h \
	$(HTSDIR)/htslib/kstring.h \
	$(HTSDIR)/htslib/sam.h \
	$(HTSDIR)/htslib/synced_bcf_reader.h \
	$(HTSDIR)/htslib/tbx.h \
	$(HTSDIR)/htslib/vcf.h \
	$(HTSDIR)/htslib/vcf_sweep.h \
	$(HTSDIR)/htslib/vcfutils.h

HTSLIB_ALL = \
	$(HTSLIB_PUBLIC_HEADERS) \
	$(HTSDIR)/bgzf.c \
	$(HTSDIR)/faidx.c \
	$(HTSDIR)/hfile_internal.h \
	$(HTSDIR)/hfile.c \
	$(HTSDIR)/hfile_net.c \
	$(HTSDIR)/hts.c \
	$(HTSDIR)/knetfile.c \
	$(HTSDIR)/kstring.c \
	$(HTSDIR)/sam.c \
	$(HTSDIR)/synced_bcf_reader.c \
	$(HTSDIR)/tbx.c \
	$(HTSDIR)/vcf.c \
	$(HTSDIR)/vcf_sweep.c \
	$(HTSDIR)/vcfutils.c \
	$(HTSDIR)/cram/cram.h \
	$(HTSDIR)/cram/cram_codecs.c \
	$(HTSDIR)/cram/cram_codecs.h \
	$(HTSDIR)/cram/cram_decode.c \
	$(HTSDIR)/cram/cram_decode.h \
	$(HTSDIR)/cram/cram_encode.c \
	$(HTSDIR)/cram/cram_encode.h \
	$(HTSDIR)/cram/cram_index.c \
	$(HTSDIR)/cram/cram_index.h \
	$(HTSDIR)/cram/cram_io.c \
	$(HTSDIR)/cram/cram_io.h \
	$(HTSDIR)/cram/cram_samtools.c \
	$(HTSDIR)/cram/cram_samtools.h \
	$(HTSDIR)/cram/cram_stats.c \
	$(HTSDIR)/cram/cram_stats.h \
	$(HTSDIR)/cram/cram_structs.h \
	$(HTSDIR)/cram/files.c \
	$(HTSDIR)/cram/mFILE.c \
	$(HTSDIR)/cram/mFILE.h \
	$(HTSDIR)/cram/md5.c \
	$(HTSDIR)/cram/md5.h \
	$(HTSDIR)/cram/misc.h \
	$(HTSDIR)/cram/open_trace_file.c \
	$(HTSDIR)/cram/open_trace_file.h \
	$(HTSDIR)/cram/os.h \
	$(HTSDIR)/cram/pooled_alloc.c \
	$(HTSDIR)/cram/pooled_alloc.h \
	$(HTSDIR)/cram/sam_header.c \
	$(HTSDIR)/cram/sam_header.h \
	$(HTSDIR)/cram/string_alloc.c \
	$(HTSDIR)/cram/string_alloc.h \
	$(HTSDIR)/cram/thread_pool.c \
	$(HTSDIR)/cram/thread_pool.h \
	$(HTSDIR)/cram/vlen.c \
	$(HTSDIR)/cram/vlen.h \
	$(HTSDIR)/cram/zfio.c \
	$(HTSDIR)/cram/zfio.h

$(HTSDIR)/libhts.a: $(HTSLIB_ALL)
	+cd $(HTSDIR) && $(MAKE) lib-static

$(HTSDIR)/libhts.so $(HTSDIR)/libhts.dylib: $(HTSLIB_ALL)
	+cd $(HTSDIR) && $(MAKE) lib-shared

$(HTSDIR)/bgzip: $(HTSDIR)/bgzip.c $(HTSLIB_PUBLIC_HEADERS)
	+cd $(HTSDIR) && $(MAKE) bgzip

$(HTSDIR)/tabix: $(HTSDIR)/tabix.c $(HTSLIB_PUBLIC_HEADERS)
	+cd $(HTSDIR) && $(MAKE) tabix

# Rules for phony targets.  You may wish to have your corresponding phony
# targets invoke these in addition to their own recipes:
#
#	clean: clean-htslib

clean-htslib install-htslib:
	+cd $(HTSDIR) && $(MAKE) $(@:-htslib=)

.PHONY: clean-htslib install-htslib
