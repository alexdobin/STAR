# Makefile variables useful for third-party code using htslib's public API.
#
#    Copyright (C) 2013-2014 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>

# These variables can be used to express dependencies on htslib headers.
# See htslib.mk for details.

htslib_bgzf_h = $(HTSPREFIX)htslib/bgzf.h
htslib_faidx_h = $(HTSPREFIX)htslib/faidx.h
htslib_hfile_h = $(HTSPREFIX)htslib/hfile.h $(htslib_hts_defs_h)
htslib_hts_h = $(HTSPREFIX)htslib/hts.h
htslib_hts_defs_h = $(HTSPREFIX)htslib/hts_defs.h
htslib_sam_h = $(HTSPREFIX)htslib/sam.h $(htslib_hts_h)
htslib_synced_bcf_reader_h = $(HTSPREFIX)htslib/synced_bcf_reader.h $(htslib_hts_h) $(htslib_vcf_h) $(htslib_tbx_h)
htslib_tbx_h = $(HTSPREFIX)htslib/tbx.h $(htslib_hts_h)
htslib_vcf_h = $(HTSPREFIX)htslib/vcf.h $(htslib_hts_h) $(HTSPREFIX)htslib/kstring.h
htslib_vcf_sweep_h = $(HTSPREFIX)htslib/vcf_sweep.h $(htslib_hts_h) $(htslib_vcf_h)
htslib_vcfutils_h = $(HTSPREFIX)htslib/vcfutils.h $(htslib_vcf_h)
