STAR 2.4.5b 2015/11/30
======================

Bug-fix release:

Fixed a problem with non-primary alignment flags with --outSAMmultNmax option.
Added counting of chimeric reads into Log.final.out .
Fixed a bug in --outSAMfilter KeepOnlyAddedReferences.
Fixed a minor bug that caused rare seg-faults.
Fixed a minor bug in STARlong extension at the ends of the read.
Fixed a seg-fault that occurred when non-default value of --genomeChrBinNbits was used.
Fixed a seg-fault that occurred when junctions where inserted after inserting reference sequences.


STAR 2.4.5a 2015/11/06
======================

**STAR now uses essential c++11 features and requires gcc 4.7.0 or later.**

Major new features:
-------------------
* Implemented on the fly insertion of the extra sequences into the genome indexes.
* Implemented --outSAMmultNmax parameter to limit the number of output alignments for multimappers.
* Implemented --outMultimapperOrder Random option to output multiple alignments in random order. 
    This also randomizes the choice of the primary alignment. Parameter --runRNGseed can be used to set the random generator seed.
    With this option, the ordering of multi-mapping alignments of each read, and the choice of the primary alignment will vary from run to run, unless only one thread is used and the seed is kept constant.

Minor new features:
-------------------
* Implemented --outSAMattrIHstart parameter. Setting it to 0 may be required for compatibility with downstream software such as Cufflinks or StringTie.
* Implemented --outSAMfilter KeepOnlyAddedReferences option.
* Implemented --help option - thanks to @yhoogstrate for the code.
* Implemented --alignEndsType Extend3pOfRead1 option for full extension of the 3' end of read 1.
* Implemented --alignSJstitchMismatchNmax option to allow for mismatches around non-canonical junctions.
* Implemented --chimSegmentReadGapMax parameter which defines the maximum gap in the read sequence between chimeric segments. By default it is set to 0 to replicate the behavior of the previous STAR versions.
* Implemented --chimFilter banGenomicN | None options to prohibit or allow the N characters in the vicinity of the chimeric junctions. By default, they are prohibited - the same behavior as in the previous versions.

Bug fixes:
----------
* For STARlong, increased compilation-time max read length to 500000 and max number of exons to 1000
* Fixed a bug which caused problems in some cases of genome generation without annotations.
* Fixed a bug in the --alignEndsType Extend5pOfRead1 option.

Code improvements:
------------------
* Improved compilation flags handling in Makefile - thanks to Christian Krause for the code.
* Improved treatment of the streams and files - thanks to Alex Finkel for the code.
* Merged pull request from Nathan S. Watson-Haigh: Makefile for manual;Travis-CI automated build; Update STAR-Fusion submodule to v0.3.1
* Merged pull request from Alex Finkel to allow 'parameter=value' option formatting, e.g. --runThreadN=8.



STAR 2.4.2a 2015/06/19
======================

* Implemented --quantMode GeneCounts option for counting number of reads per gene, similar to htseq-count.
* STARlong: fixed --outFilterIntronMotifs and --outSAMstrandField options.
* Yet another fix for --sjdbOverhang logic.
* Error message when shared memory and on the fly junction insertion are used together.
* Fixed a bug causing unnecessary 1 base soft-clipping in rare cases with sparse suffix array.
* Fixed a bug that caused problems with junction motifs in rare cases. Very few alignments affected, <1 per million.

STAR 2.4.1d 2015/05/19
======================

* Fixed problems with --sjdbOverhang default and user-defined values.
* Fixed problems with occasional non-adjacent output of multiple alignments into the unsorted BAM file and transcriptome BAM file.
* Fixed a bug causing seg-faults when shared memory options in --genomeLoad are used with --outStd SAM.
* Fixed a bug causing seg-faults for small values of --limitIObufferSize.
* Added STAR long pre-compiled executables.
* Fixed very minor issues with filtering into SJ.out.tab .
* Fixed some bugs in STARlong mapping algorithm.
* Fixed --outFilter BySJout filtering for STARlong.
* Fixed XS attrbutes in STARlong.
* Added --runDirPerm option for permissions of run-time directories.


STAR 2.4.1c 2015/04/24
======================

* Added latest version of STAR-Fusion as a separate directory.
* Fixed some compilation problems introduced in 2.4.1b.
* Added Mac executable.


STAR 2.4.1b 2015/04/23
======================

* Fixed a bug introduced in 2.4.1a causing serious problems for genomes generated without annotations.
      If you generated a genome without annotations with 2.4.1a please re-generate it.
* Fixed a bug causing seg-faults when generating genomes with a large (>500k) number of junctions.
* Fixed a bug causing seg-faults with --chimOutType WithinBAM for single-end reads.
* Fixed a bug with required --sjdbOverhang at the mapping step.

STAR 2.4.1a 2015/04/17
======================

* The annotations can now be included on the fly at the mapping step, without including them at the genome generation step.
* New option to activate on the fly "per sample" 2-pass method: "--twopassMode Basic".
* 2-pass mode can now be used with annotations, which can be included either at the run-time, or at the genome generation step.
* Included link (submodule) to Brian Haas' STAR-Fusion code for detecting fusion transcript from STAR chimeric output:  https://github.com/STAR-Fusion/STAR-Fusion
* Included Gery Vessere's shared memory implementation for POSIX and SysV. To compile STAR with POSIX shared memory, use `make POSIXSHARED`
* New option "--chimOutType WithinBAM" to include chimeric alignments together with normal alignments in the main (sorted or unsorted) BAM file(s).
* New option "--quantTranscriptomeBan Singleend" allows insertions, deletions ans soft-clips in the transcriptomic alignments, which are allowed by some expression quantification software (e.g. eXpress). 
* New option "--alignEndsTypeExtension Extend5pOfRead1" to enforce full extension of the 5p of the read1, while all other ends undergo local alignment and may be soft-clipped.

2.4.0k 03/30/2015
=================

* Implemented new BAM sorting algorithm that reduces significantly the required RAM.

2.4.0j 02/04/2015
=================

* Fixed a problem with scoring alignments for STARlong. STARlong alignments are slightly modified.
* Fixed a bug introduced in 2.4.0i that dropped a large number of aligmnents for --quantMode TranscriptomeSAM.
	Transcriptome alignments are now the same as in version 2.4.0h.
* Fixed a problem with lower case read sequences for --outSAMtype BAM options.
* Fixed a bug preventing parameter value to be "-".
* Fixed --genomeLoad LoadAndRemove option.

2.4.0i 01/14/2015
=================

* Fixed a bug with the _STARtmp temporary directory name for the 2-pass runs.
* Fixed a bug causing seg-faults for genome generation.
* Fixed a bug causing seg-faults for --quantMode TranscriptomeSAM

2.4.0h 12/09/2014
=================

* Fixed the problem causing Ubuntu error: "sh: 1: Syntax error: Bad fd number".
* Added --quantTranscriptomeBAMcompression option.
* Add newline at the end of STAR_VERSION string (contributed by Daniel Nicorici).
* Fixed a bug with parsing the last line of paired FASTA files (contributed by Alex Rolfe).

2.4.0g 11/26/2014
=================

* Fixed a bug with output score (AS attribute) of some chimeric alignments.
* Added --alignSoftClipAtReferenceEnds No option which prevents soft clipping of alignments at the reference (chromosome) ends, for compatibility with Cufflinks/Cuffmerge.
* Fixed wrong output of chimeric junctions boundaries in the Chimeric.out.junction file.
* Added --outSAMflagOR, --outSAMflagAND options to set specific bits of the SAM FLAG.
* --sjdbFileChrStartEnd can now accept multiple files which will be concatenated.
* Fixed the header of the Log.progress.out .
* Fixed a bug that prevented output of transcriptomic alignments (--quantMode TranscriptomeSAM)  with 1 base junction overhangs.
* Added --sysShell option to specify path to bash, in cases where bash is not the default shell.
* --outBAMcompression default changed to 1, which apparently does not change deflation ratio, but is much faster.
* Added --outBAMsortingThreadN option to specify number of threads for BAM sorting. By default (0) it's equal to min(6,runThreadN).

2.4.0f1 10/30/2014
==================

* Added read group (RG) BAM attributes to transcriptome BAM (contributed by https://github.com/godotgildor).
* Fixed a bug with double ID field in the RG header line (contributed by https://github.com/godotgildor).
* Fixed a bug in the 2-pass method (--twopass1readsN).
* Fixed a problem with RAM allocation for BAM sorting.

2.4.0e 10/24/2014
=================

* Added manual in PDF.
* New sub-directories: source, bin, doc.
* Output more information about read files into Log.out.
* Fixed some issues that may have caused dropping of multiple reads files.
* Added more thorough error checking for genome generation.
* Fixed a bug causing seg-faults with single-mate alignments for BAM sorting.
* Fixed some compilation issues on Mac OS X. Note that the default Clang lacks openMP support which is required for STAR compilation.
* Added Mac OS X executable.

2.4.0d 09/25/2014
=================

* Added .gitignore.
* Fixed the problem with 2nd field in the read name shorter than 3 bases (non-Illumina fastq).
* Added --outBAMcompression option.
* Added --bamRemoveDuplicatesType and --bamRemoveDuplicatesMate2basesN options.
* Added --outWigType wiggle read1_5p read2 options.
* Added --outWigNorm option.

2.4.0c 09/07/2014
=================

* Automated git version.
* Fixed a problem with overflowing SJ buffer.
* Implemented options --twopass1readsN, --twopassSJlimit, --readMapNumber.

2.4.0b 08/29/2014
=================

* Fixed problems with spaces in --outFilePrefixName.
* Fixed version information.

2.4.0a 08/11/2014
=================

* Implemented --outFilterMismatchNoverReadLmax option for a more consistent control of mismatches.

2.3.1z16 08/05/2014

Implemented --outWigReferencesPrefix option to filter references in the signal output.
Implemented --runMode inputAlignmentsFromBAM --inputBAMfile

2.3.1z15  
Implemented --outWigType bedGraph read1_5p option.
Fixed a problem with chimeric alignments with overlapping segments.
Fixed a problem with processing of fasta read input.

2.3.1z14 07/24/2014
Implemented 0x200 SAM flag for reads that did not pass Illumina filtering (i.e. contain “Y” as the 3rd character in the second field of the read name)
Implemented comma-separated lists in the --outSAMattrRGline read groups that will assign different read groups to multiple comma-separated read files in --readFilesIn

2.3.1z13 07/04/2014
Fixed problems with STARlong.

2.3.1z12 06/30/2014
Fixed problems with SAM/BAM output to stdout.

2.3.1z11 06/27/2014
Switched to htslib samtools library.
Fixed problem with indel near known splice junctions.
Fixed problem with FASTA reads input.

2.3.1z10 06/20/2014
Fixed problem with compilation, samtools/ZLIB related.

2.3.1z9 06/19/2014
2.3.1z8 
2.3.1z7 
Fixed problems with transcriptomic output.
Changed --sjdbFileChrStartEnd importing to allow direct import from SJ.out.tab

2.3.1z6 05/30/2014
2.3.1z5 05/30/2014
Fixed a bug causing problems with multiple zipped input files.
Preliminary release of BAM sorting and wiggle output

2.3.1z4 05/06/2014
Preliminary release with transcriptome output.

2.3.1z2 04/29/2014
Fixed a bug causing problems in some chimeric alignments.
Fixed a bug causing overflowing of SAM ISIZE.
Fixed chimeric output problems with --outFilterType BySJout option
Added extra Log.out messages for multi-threaded jobs.

2.3.1z1 03/13/2014
SAM header changes: 
    removed "cl:" attribute from the @PG line, output it as a separate comment line
    added --outSAMheaderHD, --outSAMheaderPG, --outSAMheaderCommentFile options
Added --outTmpDir, which set the path forSTAr temporary directory independent of --outFileNamePrefix

2.3.1z 02/05/2014
Fixed the incorrect behavior of --genomeLoad LoadAndRemove option.

2.3.1y 01/24/2014
Added read group sam attribute via --outSAMattrRGline parameter.
Fixed gcc 4.7.0. compilation problem.
Correct reverse complementarity of all IUPAC nucleotide codes in the SAM output.

2.3.1x 01/08/2014
Fixed a bug with --alignEndsType EndToEnd.

2.3.1v 12/21/2013
Added --alignEndsType EndToEnd option to align reads end-to-end, i.e. prohibit soft-clipping.
--outSAMattributes now allows to specify the SAM attributes in the desired combination and order.
Implemented standard (samtools-like) NM and MD tags.
Added --outSAMmapqUnique parameter (=255 by default), MAPQ value for unique mappers.

2.3.1u 11/23/2013
Added --outSAMreadID={Standard,Number} parameter to output read numbered read IDs.
Aded --outSAMmode NoQS option to suppress output of quality scores.

2.3.1t 11/20/2013
Fixed a bug that prevented alignment to the very beginning of the first reference.

2.3.1s 11/06/2013
Fixed a bug that produced incorrect placement of short deletions.

2.3.1r 10/01/2013
Compilation option to output "local alignment chanins".
Compilation option to output suffix array as a text file.

2.3.1q 08/15/2013
Fixed a problem with junction overhang in SJ.out.tab file for overlapping mates.

2.3.1p 04/13/2013
Fixed GCC 4.7 compatibility problems.
Changed min memory requirement for genome generation.

2.3.1o 04/13/2013
Fixed a bug with comma separated lists of input files.

2.3.1n 04/30/2013
Replaced incorretly released 2.3.1m.

2.3.1m 04/24/2013
Fixed a bug which in some cases caused problems with long reads.

2.3.1l 04/15/2013
Fixed a problem with --readFilesCommand.

2.3.1k 04/15/2013
Fixed chimeric output for single-end reads.
Fixed a problem with --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 settings.

2.3.1j 04/11/2013
Allowed spaces in paths: paths that contain spaces should be quoted with " ". Thanks for Tyler Hyndman for suggesting this patch.

2.3.1i 04/10/2013
Fixed problems with overflowing SJ buffer, new input parameters: --limitOutSJcollapsed, --limitOutSJoneRead.

2.3.1h 04/02/2013
Prevent output of strangely overlapping mates as chimeras.
Report error if --sjdbOverhang=0 with set --sjdbFileChrStartEnd or --sjdbGTFfile.

2.3.1g 03/29/2013
Implemented detection of proximal same-strand chimeras. Now it is possible to detect circular RNA (conisdered "chimeric").
More accurate treatment of overlapping mates.

2.3.1f 03/21/2013
New option --outSAMorder PairedKeepInputOrder to output alignments in the same order as they appear in the input FASTQ/A files.

2.3.1e 03/18/2013
Fixed possible problems with multi-threaded runs for small files which could have caused empty Chimeric.* and Unmapped.* output on some systems.

2.3.1d 03/17/2013
New option --outSAMprimaryFlag AllBestScore for marking all alignments with the best score as primary.
New parameter --limitOutSAMoneReadBytes 100000, limits the size of one SAM entry - important when a large number of multimappers is recorded.
Fixed a possible problem with Unmapped.* and Chimeric.* output which could generate empty or truncated output on some systems.
Coded a safer removal of the temporary directory _tmp which could have failed on some systems.
Fixed a bug which resulted in unexpected behavior for alignIntronMax < 7.

2.3.1c 03/01/2013
Fixed a bug which duplicated output in Chimeric.* and Unmapped.* when --outFilterType BySJout option is used.

2.3.1b 02/28/2013
Fixed possible issue which in some cases could have resulted in empty Chimeric.out.*

2.3.1a 02/25/2013
Fixed incorrect processing of --sjdbGTFchrPrefix.

2.3.0e

