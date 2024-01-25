STAR 2.7.11b --- 2024/01/24 ::: Minor in one parameter.
===========================================
* Replaced --quantTranscriptomeBan parameter with --quantTranscriptomeSAMoutput with more explicit naming of options. The default behavior is not affected.
* New option: --quantTranscriptomeSAMoutput BanSingleEnd_ExtendSoftclip : prohibit single-end alignments, extend softclips, allow indels.

STAR 2.7.11a --- 2023/08/15 ::: STARdiploid
===========================================
* Implemented STARdiploid option --genomeTransformType Diploid that generates personal diploid genome. At the mapping step, --genomeTransformOutput options will transform the alignments into reference genome coordinates.
* Implemented --soloCBtype String option for arbitrary cell barcode strings without passlist.
* Implemented STARsolo BAM tag sF, which outputs the feature type and number of genes for each read.
* Fixed a GstrandBit issue for the on-the-fly insertion of extra genomic sequences.
* Issue #1909: Fixed a bug causing wrong sequence length in the UB SAM tag for --soloType CB_UMI_Complex.
* Issue #1871: Fixed a bug which occurs when --soloCellReadStats Standard --twopassMode Basic are used together.
* Issue #1763: Fixed a bug causing segnmentation fault in rare cases for 2-pass mapping.
* Issue #1733: Fixed the issue with counting of intronicAS reads in STARsolo CellReads.stats output with --soloFeatures GeneFull_Ex50pAS option.
* Behavior change: for --wasp* ouput, the homozygous SNVs are filtered out from the VCF file.

STAR 2.7.10b --- 2022/11/01 ::: Bug-fix release.
===========================================================================
* PR #1638: Increased entropy of shmKey to avoid collisions between genomes. Many thanks to Jeff Hussmann (@jeffhussmann).
* Issue #1577: Reduced RAM usage for large STARsolo runs.
* Issue #1612: fixed the problem with Solo.out/SJ/raw/features.tsv sym-link.
* Issues #1469, #1602, #1608: fixed seg-faults introduced in 2.7.10a.
* Issue #1619: fixed memory leak in SoloFeature_cellFiltering.cpp
* Issue #1543: fixed a segfault occurring for STARsolo multimappers and large number of reads per cell.
* Issue #1558: fixed a bug with output of GX/GN BAM tags without CB/UB.
* Issue #1513: if --soloMultiMappers options are not requested, output "NoMulti" in the "Reads Mapped to Gene: Unique+Multiple Gene" line of the Summary.csv file.
* Issue #719: implemented auto re-allocation of the SJ output buffer.
* Fixed a bug with --soloMultiMappers for small number of cells cases.
* Fixed a problem with STARsolo CellReads.stats output for no-passlist runs.

STAR 2.7.10a --- 2022/01/14 ::: New features, behavior changes and bug fixes
===========================================================================
**New options and features:**
* Implemented --soloCellReadStats Standard option to output read statistics for each cell barcode.
* Allow to define --clip5pAdapterSeq with --clipAdapterType CellRanger4 option.
* Implemented --soloCBmatchWLtype ED2 to allow mismatches and one insertion+deletion (edit distance <=2) for --soloType CB_UMI_Complex.
* Implemented Solo BAM tags gx gn: output ';'-separated gene IDs and names for both unique- and multi-gene reads. Note that GX/GN tags are used to output gene ID/name for unique-gene reads.
* Implemented --soloFeatures GeneFull_ExonOverIntron GeneFull_Ex50pAS options which prioritize exonic over intronic overlaps for pre-mRNA counting.
* Added script extras/scripts/soloCountMatrixFromBAM.awk to re-create Solo count matrix from the BAM output.

**Changes in behavior:**
* Changed --soloType CB_samTagOut behavior: if barcode cennot be matched to the passlist, CB:Z:- will be recorded (previously CB tag was absent for such reads).
* Changed Solo summary statistics outputs in Barcodes.stats and Features.stats files.
* Changed Solo BAM tags GX GN behavior: for missing values, "-" is output instead of omitting the tag.
* Changed Solo BAM tags output for multiple --soloFeatures: now the first feature on the list is used for GX,GN,XB,UB tags.
* Changed Solo SJ behavior: it no longer depends on the whether the alignment is concordant to a Gene.
* Fixed a bug that resulted in slightly different solo counts if --soloFeatures Gene and GeneFull were used together with --soloCBmatchWLtype 1MM_multi_pseudocounts option.

**Bug fixes**

* PR #1425: Assign supplementary alignment to correct mate when mates fully overlap. Many thanks to Sebastian @suhrig for resolving this problem in the chimeric detection.
* Fixed a bug introduced in 2.7.9a for --quantMode TranscriptomeSAM output that resulted in both mapped and unmapped output for some reads. Many thanks to Diane Trout (@Caltech) for helping to track this bug.
* Issue #1223: fixed the N_unmapped value reported in ReadsPerGene.out.tab. The single-end (i.e. partially mapped alignment are not excluded from N_unmapped.
* Issues #535, #1350: fixed a long-standing problem that resulted in a seg-fault whem mapping to the rabbit genome.
* Issue #1316: fixed the seg-fault which occurred if --soloType CB_samTagOut and --soloCBwhitelist None are used together.
* Issue #1177: throw an error in case the BAM file does not contain NH and AS tags for duplication removal jobs (--runMode inputAlignmentsFromBAM --bamRemoveDuplicatesType UniqueIdenticalNotMulti).
* Issue #1262: fixed the bug that prevented EM matrix output when only EM option is specified in --soloMultiMappers.
* Issue #1230: fixed the bug that caused seg-faults for --runMode soloCellFiltering runs.

STAR 2.7.9a --- 2021/05/05 ::: STARsolo updates
=====================================================
**Major updates:**
* STARsolo can perform counting of multi-gene (multi-mapping) reads with --soloMultiMappers EM [Uniform Rescue PropUnqiue] options.
* PR #1163: [SIMDe](https://github.com/simd-everywhere/simde) takes care of correct SIMD extensions based on -m g++ flag: compilation option CXXFLAGS_SIMD is preset to -mavx2, but can be to the desired target architecture. Many thanks to Michael R. Crusoe @mr-c, Evan Nemerson @nemequ and Steffen Möller @smoe!

**New options and features:**
* New option: --soloUMIfiltering MultiGeneUMI_All to filter out all UMIs mapping to multiple genes (for uniquely mapping reads)
* New script extras/scripts/calcUMIperCell.awk to calculate total number of UMIs per cell and filtering status from STARsolo matrix.mtx
* New option: --outSJtype None   to omit outputting splice junctions to SJ.out.tab
* Simple script to convert BED spliced junctions (SJ.out.tab) to BED12 for UCSC display: extras/scripts/sjBED12.awk
* PR #1164: SOURCE_DATE_EPOCH to make the build more reproducible
* PR #1157: print STAR command line and version information to stdout

**Changes in behavior:**
* Minor changes to statistics output (Features.csv and Summary.csv) to accomodate multimappers.
* Modified option: ---limitIObufferSize now requires two numbers - separate sizes for input and output buffers

**Bug fixes**
* PR #1156: clean opal/opal.o
* Issue #1166: seg-fault for STARsolo --soloCBwhitelist None (no whitelist) with barcodes longer than 16b
* Issue #1167: STARsolo CR/UR SAM tags are scrambled in TranscriptomeSAM file Aligned.toTranscriptome.out.bam. This bug appeared in 2.7.7a.
* Issue #1177: Added file checks for the --inputBAMfile .
* Issue #1180: Output the actual number of alignments in NH attributes even if --outSAMmultNmax is set to a smaller value.
* Issue #1190: Allow GX/GN output for non-STARsolo runs.
* Issue #1220: corrupt SAM/BAM files for --outFilterType BySJout. The bug was introduced with the chnages in 2.7.7a.
* Issue #1211: scrambled CB tags in BAM output for --soloCBwhitelist None --soloFeatures Gene GeneFull.
* Fixed a bug causing seg-faults with --clipAdapterType CellRanger4 option.

STAR 2.7.8a --- 2021/02/20 ::: Major STARsolo updates
=====================================================

**Major new features:**
* ```--runMode soloCellFiltering``` option for cell filtering (calling) of the raw count matrix, without re-mapping
* Input from SAM/BAM for STARsolo, with options ```--soloInputSAMattrBarcodeSeq``` and ```--soloInputSAMattrBarcodeQual``` to specify SAM tags for the barcode read sequence and qualities
* ```--clipAdapterType CellRanger4``` option for 5' TSO adapter and 3' polyA-tail clipping of the reads to better match CellRanger >= 4.0.0 mapping results
* ```--soloBarcodeMate``` to support scRNA-seq protocols in which one of the paired-end mates contains both barcode sequence and cDNA (e.g. 10X 5' protocol)

**New options:**
* ```--soloCellFilter EmptyDrops_CR``` option for cell filtering (calling) nearly identical to that of CellRanger 3 and 4
* ```--readFilesSAMattrKeep``` to specify which SAM attributes from the input SAM to keep in the output
* ```--soloUMIdedup 1M_Directional_UMItools``` option matching the "directional" method in UMI-tools Smith, Heger and Sudbery (Genome Research 2017)
* ```--soloUMIdedup NoDedup``` option for counting reads per gene, i.e. no UMI deduplication
* ```--soloUMIdedup 1MM_CR``` option for 1 mismatch UMI deduplication similar to CellRanger >= 3.0
* ```--soloUMIfiltering MultiGeneUMI_CR``` option filters lower-count UMIs that map to more than one gene matching CellRanger >= 3.0
* ```--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts``` options which allows 1MM multimatching to WL for barcodes with N-bases (to better match CellRanger >= 3.0)

**Changes in behavior:**
* The UMI deduplication/correction specified in ```--soloUMIdedup``` is used for statistics output, filtering and UB tag in BAM output.
* If UMI or CB are not defined, the UB and CB tags in BAM output will contain "-" (instead of missing these tags).
* For ```--soloUMIfiltering MultiGeneUMI``` option, the reads with multi-gene UMIs will have UB tag "-" in BAM output.
* Different ```--soloUMIdedup``` counts, if requested, are recorded in separate .mtx files.
* Cell-filtered Velocyto matrices are generated using Gene cell filtering.
* Velocyto spliced/unspliced/ambiguous counts are reported in separate .mtx files.
* Read clipping options ```--clip*``` now require specifying the values for all read mates, even if they are identical.

**Bug fixes:**
* Issue #1107: fixed a bug causing seg-fault for ```--soloType SmartSeq``` with only one (pair of) fastq file(s)
* Issue #1129: fixed an issue with short barcode sequences and ```--soloBarcodeReadLength 0```
* Issue  #796: Fixed a problem with GX/GN tag output for ```--soloFeatures GeneFull``` option
* PR: #1012: fix the bug with ```--soloCellFilter TopCells``` option
* Fixed an issue that was causing slightly underestimated value of Q30 'Bases in RNA read' in ```Solo.out/Gene/Summary.csv```

STAR 2.7.7a --- 2020/12/28
==========================
**Major new feature: STARconsensus: mapping RNA-seq reads to consensus genome.**
* Insert (consensus) variants from a VCF file into the reference genome at the genome generation step with ```--genomeTransformVCF Variants.vcf --genomeTransformType Haploid```
* Map to the transformed genome. Alignments (SAM/BAM) and spliced junctions (SJ.out.tab) can be transformed back to the original (reference) coordinates with ```--genomeTransformOutput SAM and/or SJ```

**Minor bug fixes:**
* Deprecated ```--genomeConsensusFile``` option. Please use ```--genomeTransformVCF``` and ```--genomeTransformType``` options instead.
* Issue #1040: fixed a bug causing rare seg-faults for paired-end --soloType SmartSeq runs.
* Issue #1071: fixed a bug that can cause a crash for STARsolo runs with a small number of cells.

STAR 2.7.6a --- 2020/09/19
==========================
**Major new feature:**
Output multimapping chimeric alignments in BAM format using
```--chimMultimapNmax N>1 --chimOutType WithinBAM --outSAMtype BAM Unsorted [and/or] SortedByCoordinate```
Many thanks to Sebastian @suhrig who implemented this feature!
More detailed description from Sebastian in PR #802.

**Minor features and bug fixes:**
* Issue #1008: fixed the problem with Unmapped.out.mate? output for --soloType CB_samTagOut output.
* PR # 1012: fixed the bug with --soloCellFiltering TopCells option.
* Issue #786: fixed the bug causing the *Different SJ motifs problem* for overlapping mates.
* Issue #945: GX/GN can be output for all --soloType, as well as for non-solo runs.

STAR 2.7.5c --- 2020/08/16
==========================
Bug-fix release.
----------------

* Issue #988: proceed reading from GTF after a warning that exon end is past chromosome end.
* Issue #978: fixed corrupted transcriptInfo.tab in genome generation for cases where GTF file contains extra chromosomes not present in FASTA files.
* Issue #945: otuput GX/GN for --soloFeatures GeneFull .
* Implemented removal of control characters from the ends of input read lines, for compatibility with files pre-processed on Windows.

STAR 2.7.5b --- 2020/08/01
==========================
Bug-fix release.
----------------

* Issue #558: Fixed a bug that can cause a seg-fault in STARsolo run with paired-end reads that have protruding ends.
* Issue #952: Increased the maximum allowed length of the SAM tags in the input SAM files.
* Issue #955: fixed seg-fault-causing bug for --soloFeatures SJ option.
* Issue #963: When reading GTF file, skip any exons that extend past the end of the chromosome, and give a warning.
* Issue #965: output genome sizes with and without padding into Log.out.
* Docker build: switched to debian:stable-slim in the Dockerfile.
* --soloType CB_samTagOut now allows output of (uncorrected) UMI sequences and quality scores with SAM tags UR and UY.
* Throw an error if FIFO file cannot be created on non-Linux partitions.

STAR 2.7.5a 2020/06/16
======================
**Major new features:**
* Implemented STARsolo quantification for Smart-seq with --soloType SmartSeq option.
* Implemented --readFilesManifest option to input a list of input read files.

**Minor features and bug fixes:**
* Change in STARsolo SJ output behavior: junctions are output even if reads do not match genes.
* Fixed a bug with solo SJ output for large genomes.
* N-characters in --soloAdapterSequence are not counted as mismatches, allowing for multiple adapters (e.g. ddSeq).
* SJ.out.tab is sym-linked as features.tsv for Solo SJ output.
* Issue #882: 3rd field is now optional in Solo Gene features.tsv with --soloOutFormatFeaturesGeneField3.
* Issue #883: Patch for FreeBSD in SharedMemory and Makefile improvements.
* Issue #902: Fixed seg-fault for STARsolo CB/UB SAM attributes output with --soloFeatures GeneFull --outSAMunmapped Within options.
* Issue #934: Fixed a problem with annotated junctions that was causing very rare seg-faults.
* Issue #936: Throw an error if an empty whitelist is provided to STARsolo.

STAR 2.7.4a 2020/06/01
======================
Fixing multiple bugs and issues.
--------------------------------

**This version requires re-generation of the genome indexes**

* Fixed the long-standing seg-fault problem for small genomes.
* Issue #784: Fixed a seg-fault in STARsolo for cases where no cell barcodes matched whitelist.  
* Issue #798: Fixed the problem in Solo Q30 Bases in Summary.csv average (#798).
* Issue #843, #880: Throw an error if read file in --readFilesIn does not exist when using --readFilesCommand .
* Issue #864: Fixed seg-fault for STARsolo runs with very small number of reads or cells.
* Issue #881: Check if --genomeDir exists, create if necessary.
* Issue #882: Added 3rd column "Gene Expression" to solo features.tsv file for better compatibility with downstream tools.
* Issue #902: Fixed seg-fault for STARsolo CB/UB SAM attributes output with --soloFeatures GeneFull only option.
* Issue #907: Fixed the bug that prevented output of STARsolo GX/GN tags into the Aligned.out.bam if --quantMode TranscriptomeSAM is used.
* Issue #910: The output directory in --outFileNamePrefix is checked and created if it does not exist.
* If solo barcode read length is not checked (--soloBarcodeReadLength 0) and it is shorter than CB+UMI length, the barcode is padded with Ns and not counted.
* For genome generation runs, the Log.out file is moved into the --genomeDir directory.
* Fixed a bug with solo SJ output for large genomes.
* Implemented --seedMapMin option (previously hard-coded) to define minimum seed length.

STAR 2.7.3a 2019/10/08
======================
Major new features in STARsolo
------------------------------
* **Output enhancements:**
    * Summary.csv statistics output for raw and filtered cells useful for quick run quality assessment.
    * --soloCellFilter option for basic filtering of the cells, similar to the methods used by CellRanger 2.2.x.
* [**Better compatibility with CellRanger 3.x.x:**](docs/STARsolo.md#matching-cellranger-3xx-results)
    * --soloUMIfiltering MultiGeneUMI option introduced in CellRanger 3.x.x for filtering UMI collisions between different genes.
    * --soloCBmatchWLtype 1MM_multi_pseudocounts option, introduced in CellRanger 3.x.x, which slightly changes the posterior probability calculation for CB with 1 mismatch.
* [**Velocyto spliced/unspliced/ambiguous quantification:**](docs/STARsolo.md#velocyto-splicedunsplicedambiguous-quantification)
    * --soloFeatures Velocyto option to produce Spliced, Unspliced, and Ambiguous counts similar to the [velocyto.py](http://velocyto.org/) tool developed by [LaManno et al](https://doi.org/10.1038/s41586-018-0414-6). This option is under active development and the results may change in the future versions.
* [**Support for complex barcodes, e.g. inDrop:**](docs/STARsolo.md#barcode-geometry)
    * Complex barcodes in STARsolo with --soloType CB_UMI_Complex, --soloCBmatchWLtype --soloAdapterSequence, --soloAdapterMismatchesNmax, --soloCBposition,--soloUMIposition
* [**BAM tags:**](#bam-tags)
    * CB/UB for corrected CellBarcode/UMI
    * GX/GN for gene ID/name
* STARsolo most up-to-date [documentation](docs/STARsolo.md).


STAR 2.7.2d 2019/10/04
======================
* Fixed the problem with no header in Chimeric.out.sam

STAR 2.7.2c 2019/10/02
======================
* Fixed the problem with no output to Chimeric.out.sam

STAR 2.7.2b 2019/08/29
======================
Bug fixes in chimeric detection, contributed by Meng Xiao He (@mengxiao)
* Fix memory leak in handling chimeric multimappers: #721
* Ensure chimeric alignment score requirements are consistently checked: #722,#723.

STAR 2.7.2a 2019/08/13
======================
* Chimeric read reporting now requires that the chimeric read alignment score higher than the alternative non-chimeric alignment to the reference genome.  The Chimeric.out.junction file now includes the scores of the chimeric alignments and non-chimeric alternative alignments, in addition to the PEmerged bool attribute. (bhaas, Aug 2019)
* Fixed the problem with ALT=* in STAR-WASP.
* Implemented extras/scripts/soloBasicCellFilter.awk script to perform basic filtering of the STARsolo count matrices.
* Fixed a bug causing rare seg-faults with for --peOverlap* options and chimeric detection.
* Fixed a problem in STARsolo with unmapped reads counts in Solo.out/*.stats files.
* Fixed a bug in STARsolo with counting reads for splice junctions. Solo.out/matrixSJ.mtx output is slighlty changed.
* Fixed the problem with ALT=* in VCF files for STAR-WASP.

STAR 2.7.1a 2019/05/15
======================
**This version requires re-generation of the genome indexes**

* Implemented --soloFeatures GeneFull which counts reads overlapping full genes, i.e. includes reads that overlap introns. This can be combined with other features, e.g. --soloFeatures Gene SJ GeneFull .
* Implemented --soloCBwhitelist None option for solo* demultiplexing without CB whitelist. In this case error correction for CBs is not performed.
* Implemented Cell Barcodes longer than 16 bases (but shorter than 31 bases). Many thanks to Gert Hulselmans for implementing this feature (#588).
* Implemented collapsing of duplicate cell barcodes in the whitelist.
* Implemented --sjdbGTFtagExonParentGeneName and --sjdbGTFtagExonParentGeneType options to load gene name and biotype attributes from the GTF file.
* Fixed problems created by missing gene/transcript ID, name and biotype attributes in GTF files (issues #613, #628).
* Added warning for incorrectly scaled --genomeSAindexNbases parameter (issue #614).
* Added numbers of unmapped reads to the Log.final.out file (pull #622).
* Fixed a problem which may cause seg-faults for reads with many blocks (issue #342).

STAR 2.7.0f 2019/03/28
======================
* Fixed a problem in STARsolo with empty Unmapped.out.mate2 file. Issue #593.
* Fixed a problem with CR CY UR UQ SAM tags in solo output. Issue #593.
* Fixed problems with STARsolo and 2-pass.

STAR 2.7.0e 2019/02/25
======================
* Fixed problems with --quantMode GeneCounts and --parametersFiles options

STAR 2.7.0d 2019/02/19
======================
* Implemented --soloBarcodeReadLength option for barcode read length not equal to the UMI+CB length
* Enforced genome version rules for 2.7.0

STAR 2.7.0c 2019/02/08
======================
* This release is compiled with gcc-4.8.5, and requires at least gcc-4.8.5
* Fixed another problem in STARsolo genes.tsv output.
* Replaced tabs with spaces in STARsolo matrix.mtx output
* #559, #562 Fixed compilation problems.
* #550 (again, previous merge failed): Added correct header for the STARsolo matrix.mtx file, needed for python scipy mmread compatibility.

STAR 2.7.0b 2019/02/05
======================
* #550: Added correct header for the STARsolo matrix.mtx file, needed for python scipy mmread compatibility.
* #556: Fixed a problem with STARsolo genes.tsv file, which may also cause troubles with GTF files processing.
* Important: 2.7.0x releases require re-generation of the genome index.

STAR 2.7.0a 2019/01/23
======================
* This release introduces STARsolo for: mapping, demultiplexing and gene quantification for single cell RNA-seq.
* Multiple solo\* options control STARsolo algorithm. See the RELEASEnotes and the manual for more information.
* This release is compiled with gcc-5.3.0, and requires at least gcc-4.9.4


STAR 2.6.1d 2018/11/16
======================

* Fixed the problem causing BAM sorting error with large number of threads and small ulimit -n (github.com/alexdobin/STAR/issues/512).
* Fixed the bug causing inconsistent output for mate1/2 in the Unmapped files (github.com/alexdobin/STAR/issues/222).
* Fixed the non-thread safe error/exit (github.com/alexdobin/STAR/issues/514), and non-safe file size check (github.com/alexdobin/STAR/issues/516)
* Many thanks to Paul Menzel for helping to track and fix these problems.


STAR 2.6.1c 2018/10/17
======================

* Enforced the consistent choice of supplementary chimeric alignments for overlapping mates.


STAR 2.6.1b 2018/09/06
======================

* Fixed a problem with --outSAMfilter KeepOnlyAddedReferences option.
* Fixed a problem with output of an empty sorted BAM.


STAR 2.6.1a 2018/08/14
======================

* Process substitution can now be used with zipped VCF files, e.g. --varVCFfile <(zcat vcf.gz)
* Implemented fatal error exception if no SNPs are found in VCF files.
* Implemented --chimOutJunctionFormat 1 option to output some metadata (command lines and basic mapping statistics) at the end of Chimeric.out.junction file.
* The default value of --peOverlapMMp is reduced to 0.01 for less aggressive mate merging.
* Fixed the problem with control characters (ASCII<32) in genome and input read sequences. They used to be converted to N, now they are removed.
* Fixed a bug that caused serious problems with --sjdbInsertSave All option.
* Fixed a bug in merging mates (--peOverlap*) algorithm that was causing rare seg-faults.
* Fixed the GtstrandBit problem.
* Fixed a bug with multiple RG lines when inputting reads in SAM format.
* Fixed a bug causing seg-faults with shared memory and --outStd options.
* Fixed a bug with --outTmpDir and fifo files.

STAR 2.6.0c 2018/05/10
======================

* Fixed bugs in merging mates (--peOverlap*) and WASP filtering algorithms. Please see CHANGES and RELEASEnotes from 2.6.0a.


STAR 2.6.0b 2018/05/02
======================

* Fixed bugs introduced in 2.6.0a. Please see CHANGES and RELEASEnotes from 2.6.0a.


STAR 2.6.0a 2018/04/23
======================

Major new features:
-------------------
* Merging and mapping of overlapping paired-end reads with new options --peOverlapNbasesMin and --peOverlapMMp. The developmment of this algorithm was supported by Illumina, Inc. Many thanks to June Snedecor, Xiao Chen, and Felix Schlesinger for their extensive help in developing this feature.
* --varVCFfile option to input variant VCF file.
* New SAM attributes in the --outSAMattributes, vG, vA, and vW to report variants overlapping alignments.
* --waspOutputMode option for filtering allele specific alignments. This is re-implementation of the original WASP algorithm by Bryce van de Geijn, Graham McVicker, Yoav Gilad & Jonathan K Pritchard. Please cite the original WASP paper: Nature Methods 12, 1061–1063 (2015), https://www.nature.com/articles/nmeth.3582 . Many thanks to Bryce van de Geijn for fruitful discussions.
* Detection of multimapping chimeras, with new options --chimMultimapNmax, --chimMultimapScoreRange and --chimNonchimScoreDropMin . Many thanks to Brian Haas for testing and feedback.


Minor new features:
-------------------
* --alignInsertionFlush option which defines how to flush ambiguous insertion positions: None: old method, insertions are not flushed; Right: insertions are flushed to the right.
* --outSAMtlen option to select the calculation method for the TLEN field in the SAM/BAM files.
* --outBAMsortingBinsN option to control the number of sorting bins. Increasing this number reduces the amount of RAM required for sorting.

STAR 2.5.4b 2018/02/09
======================

* Recompiled Linux executables with the correct version tag.
* Updated manual.

STAR 2.5.4a 2018/01/23
======================

### New features:
* Implemented read group ID output as the last column of the Chimeric.out.junction file.
* Implemented --readFilesPrefix option for specifying prefix (e.g. directory path) for the file names in --readFilesIn .
* Implemented standard SAM attribute "MC" to output the mate's CIGAR. Add MC to the list of attributes in the --outSAMattribute option.
* Implemented the ability to input the reads from unmapped SAM/BAM file: --readFilesType SAM SE[PE] for single-end [paired-end] reads to read from the SAM file specified, as usual in --readFilesIn. For BAM files, in addition, specify --readFilesCommand samtools view -h .
* Implemented --seedSplitMin option which was previously hardcoded at 12. his will allow mapping of mates shorter than 12nt.
* Implemented --outFilterIntronStrands None option to switch off filtering by strand consistency of junctions.
* Added new scripts extras/scripts/mergeLogFinal.awk,mergeSuperContig.awk,sjMotif.m

### Bug fixes:
* Fixed a bug in chimeric detection code which sometimes led to uninitialized memory access. The chimeric output may change for a very small number of reads.
* Fixed a problem with --alignEndsProtrude implementation which prevented the output of alignments with protruded ends.
* Fixed a bug which set non-primary bit 0x100 in the SAM FLAG for unmapped mates.
* Fixed a bug in liftOver command that output an extra field in the GTF file.
* Fixed a problem that can arise for very small genomes while using --alignIntronMax 1.

STAR 2.5.3a 2017/03/17
======================

* Fixed occasional seg-faults after the completion of the mapping runs with shared memory.
* Implemented --genomeFileSizes option to supply sizes of the genome index files. This allows for streaming of index files.
* Implemented extra references input in the SAM/AM header from user-created "extraReferences.txt" file in the genome directory.
* Implemented --chimOutType HardClip OR SoftClip options to output hard (default) / soft clipping in the BAM CIGAR for supplementary chimeric alignments.
* Implemented --chimMainSegmentMultNmax parameters, which may be used to prohibit chimeric alignments with multimapping main segments to reduce false positive chimeras.
* Implemented new SAM attribute 'ch' to mark chimeric aligmments in the BAM file for --chimOutType WithinBAM option.
* Fixed a problem with RNEXT field in the Chimeric.out.sam file: RNEXT now always points to the other mate start.
* Implemented --bamRemoveDuplicatesType UniqueIdenticalNotMulti option, which (unlike the UniqueIdentical optipon) will NOT mark multi-mappers as duplicates.
* For --bamRemoveDuplicatesType UniqueIdentical, the unmmapped reads are no longer marked as duplicates.

STAR 2.5.2b 2016/08/19
======================

* Fixed a problem with --outSAMmultNmax 1 not working for transcriptomic output.
* Fixed a bug with chimeric BAM output for --chimOutType WithinBAM option.
* Fixed a bug that could cause non-stable BAM sorting if the gcc qsort is unstable.
* Fixed a bug with causing seg-faults when combining  --twopassMode Basic --outSAMorder PairedKeepInputOrder .
* Fixed a problem with SAM header in cases where reference sequences are added at the mapping stage.

STAR 2.5.2a 2016/05/10
======================

* Fixed the "GstrandBit" problem.
* Fixed a bug introduced in 2.5.1a that caused problems with single-end alignments output in some cases.
* Fixed a bug that can cause STARlong seg-faults in rare cases.
* Fixed a bug that caused output of unmapped mates for single end alignments even with --outSAMunmapped None .
* Implemented --winReadCoverageRelativeMin and --winReadCoverageBasesMin to control coverage of the alignment windows for STARlong.
* Implemented --outSAMfilter KeepAllAddedReferences option which will keep all alignments to the added references.
* Implemented --alignEndsProtrude option to control output of alignments with protruding ends.
* Implemented --outTmpKeep All option to keep the temporary files.
* Implemented --alignEndsType Extend5pOfReads12 option for full extension of 5' ends of both mates.


STAR 2.5.1b 2016/01/22
======================

* Fixed a bug in signal generation with --outWigType introduced in 2.5.1a


STAR 2.5.1a 2016/01/19
======================

* Fixed a bug in --quantMode TranscriptomeSAM that prevented output to Aligned.toTranscriptome.out.bam of the reads mapped to the very last annotated transcript.
* Cleaned up the code to remove compilation warnings (thanks to github.com/yhoogstrate).
* Implemented --outSAMunmapped Within KeepPairs option to record unmapped mate adjacent to the mapped one, in case single-end alignments are allowed.
  For multi-mappers, the unmapped mate will be recored mulitple times adjacent to the mappet mate of each alignment.


STAR 2.5.0c 2015/12/23
======================

* Implemented --genomeSuffixLengthMax option to control max suffix length at the genome generation step.
* Fixed a bug that caused genome generation stalling in some cases.
* In Aligned.toTranscriptome.out.bam (--quantMode TranscriptomeSAM), non-primary SAM flag is assigned to all but one randomly selected alignment in Aligned.toTranscriptome.out.bam .
* Fixed a bug that filtered out some chimeric junctions.
* Fixed a bug that prevented chimeric output for some of the "circular" configurations.

STAR 2.5.0b 2015/11/30
======================

**Bug-fix release:**

* Fixed a problem with non-primary alignment flags with --outSAMmultNmax option.
* Added counting of chimeric reads into Log.final.out .
* Fixed a bug in --outSAMfilter KeepOnlyAddedReferences.
* Fixed a minor bug that caused rare seg-faults.
* Fixed a minor bug in STARlong extension at the ends of the read.
* Fixed a seg-fault that occurred when non-default value of --genomeChrBinNbits was used.
* Fixed a seg-fault that occurred when junctions where inserted after inserting reference sequences.


STAR 2.5.0a 2015/11/06
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
