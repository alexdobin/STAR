# Makefile variables useful for third-party code using htslib's public API.
#
#    Copyright (C) 2013-2014 Genome Research Ltd.
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

# These variables can be used to express dependencies on htslib headers.
# See htslib.mk for details.

htslib_bgzf_h = $(HTSPREFIX)htslib/bgzf.h
htslib_faidx_h = $(HTSPREFIX)htslib/faidx.h
htslib_hfile_h = $(HTSPREFIX)htslib/hfile.h $(htslib_hts_defs_h)
htslib_hts_h = $(HTSPREFIX)htslib/hts.h
htslib_hts_defs_h = $(HTSPREFIX)htslib/hts_defs.h
htslib_regidx_h = $(HTSPREFIX)htslib/regidx.h
htslib_sam_h = $(HTSPREFIX)htslib/sam.h $(htslib_hts_h)
htslib_synced_bcf_reader_h = $(HTSPREFIX)htslib/synced_bcf_reader.h $(htslib_hts_h) $(htslib_vcf_h) $(htslib_tbx_h)
htslib_tbx_h = $(HTSPREFIX)htslib/tbx.h $(htslib_hts_h)
htslib_vcf_h = $(HTSPREFIX)htslib/vcf.h $(htslib_hts_h) $(HTSPREFIX)htslib/kstring.h
htslib_vcf_sweep_h = $(HTSPREFIX)htslib/vcf_sweep.h $(htslib_hts_h) $(htslib_vcf_h)
htslib_vcfutils_h = $(HTSPREFIX)htslib/vcfutils.h $(htslib_vcf_h)
