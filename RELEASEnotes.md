STAR 2.6.0a 2018/04/23
======================

Major new features:
-------------------

**1. Merging and mapping of overlapping paired-end reads.**  
This feature improves mapping accuracy  for paired-end libraries with short insert sizes, where many reads have overlapping mates. Importantly, it allows detection of chimeric junction in the overlap region.
STAR will search for an overlap between mates larger or equal to --peOverlapNbasesMin bases with proportion of mismatches in the overlap area not exceeding --peOverlapMMp .
If the overlap is found, STAR will map merge the mates and attempt to map the resulting (single-end) sequence.
If requested, the chimeric detection will be performed on the merged-mate sequence, thus allowing chimeric detection in the overlap region.
If the score of this alignment higher than the original one, or if a chimeric alignment is found, STAR will report the merged-mate aligment instead of the original one.
In the output, the merged-mate aligment will be converted back to paired-end format.  
The developmment of this algorithm was supported by Illumina, Inc. 
Many thanks to June Snedecor, Xiao Chen, and Felix Schlesinger for their extensive help in developing this feature.


**2. Detection of personal variants overlapping alignments.**  
Option --varVCFfile /path/to/vcf/file is used to input VCF file with personal variants. Only single nucleotide variants (SNVs) are supported at the moment. 
Each variant is expected to have a genotype with two alleles.
To output variants that overlap alignments, vG and vA have to be added to --outSAMattributes list. 
SAM attribute vG outputs the genomic coordinate of the variant, allowing for identification of the variant.
SAM attribute vA outputs which allele is detected in the read: 1 or 2 match one of the genotype alleles, 3 - no match to genotype.

**3. WASP filtering of allele specific alignments.**  
This is re-implementation of the original WASP algorithm by Bryce van de Geijn, Graham McVicker, Yoav Gilad & Jonathan K Pritchard. Please cite the original [WASP paper: Nature Methods 12, 1061–1063 (2015)   ](https://www.nature.com/articles/nmeth.3582).
WASP filtering is activated with --waspOutputMode SAMtag, which will add vW tag to the SAM output: 
vW:i:1 means alignment passed WASP filtering, while all other values mean it did not pass.  
Many thanks to Bryce van de Geijn for fruitful discussions.

**4. Detection of multimapping chimeras.**  
Previous STAR chimeric detection algorithm only detected uniquely mapping chimeras, which reduced its sensitivity in some cases.
The new algorithm can detect and output multimapping chimeras. Presently, the only output into Chimeric.out.junction is supported.
This algorithm is activated with >0 value in --chimMultimapNmax, which defines the maximum number of chimeric multi-alignments.
The --chimMultimapScoreRange (=1 by default) parameter defines the score range for multi-mapping chimeras below the best chimeric score, similar to the --outFilterMultimapScoreRange parameter for normal alignments.
The --chimNonchimScoreDropMin (=20 by default) defines the threshold triggering chimeric detection: the drop in the best non-chimeric alignment score with respect to the read length has to be smaller than this value.  
Many thanks to Brian Haas for testing and feedback.


Minor new features:
-------------------
* --outSAMtlen 1/2 option to select the calculation method for the TLEN field in the SAM/BAM files:
              1 ... leftmost base of the (+)strand mate to rightmost base of the (-)mate. (+)sign for the (+)strand mate
              2 ... leftmost base of any mate to rightmost base of any mate. (+)sign for the mate with the leftmost base. This is different from 1 for overlapping mates with protruding ends
* --alignInsertionFlush option which defines how to flush ambiguous insertion positions: None: old method, insertions are not flushed; Right: insertions are flushed to the right.
* --outBAMsortingBinsN option to control the number of sorting bins. Increasing this number reduces the amount of RAM required for sorting.


STAR 2.5.0a 2015/11/06
======================

**STAR now uses essential c++11 features. Compiling from sources requires gcc 4.7.0 or later.**

Major new features:
-------------------
1. It is now possible to add extra sequences to the reference genome ont the fly (without re-generating the genome) by specifying 
_--genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2_ at the mapping stage. 

2. By default, the order of the multi-mapping alignments for each read is not truly random.
The _--outMultimapperOrder Random_ option outputs multiple alignments for each read in random order, 
and also also randomizes the choice of the primary alignment from the highest scoring alignments. 
Parameter _--runRNGseed_ can be used to set the random generator seed. 
With this option, the ordering of multi-mapping alignments of each read, 
and the choice of the primary alignment will vary from run to run, unless only one thread is used and the seed is kept constant.

3. The --outSAMmultNmax parameter limits the number of output alignments (SAM lines) for multimappers. 
For instance, _--outSAMmultNmax 1_ will output exactly one SAM line for each mapped read.


STAR 2.4.2a 2015/06/19
======================

New features:

Counting reads per gene while mapping with --quantMode GeneCounts option.
A read is counted if it overlaps (1nt or more) one and only one gene. Both ends of the paired-end read are checked for overlaps.
The counts coincide with those produced by htseq-count with default parameters.

Requires annotations (GTF or GFF with --sjdbGTFfile option) used at the genome generation step, or at the mapping step.

Outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which correspond to different strandedness options:
column 1: gene ID 
column 2: counts for unstranded RNA-seq
column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
Select the output according to the strandedness of your data.
Note, that if you have stranded data and choose one of the columns 3 or 4, the other column (4 or 3) will give you the count of antisense reads.
 
With --quantMode TranscriptomeSAM GeneCounts, and get both the Aligned.toTranscriptome.out.bam and ReadsPerGene.out.tab outputs.



STAR 2.4.1a 2015/04/17
======================

New features:

1. The annotations can now be included on the fly at the mapping step, without including them at the genome generation step.
   At the mapping step, specify --sjdbGTFfile /path/to/ann.gtf and/or --sjdbFileChrStartEnd /path/to/sj.tab, as well as --sjdbOverhang, and any other --sjdb* options.
   The genome indices can be generated with or  without another set of annotations/junctions. In the latter case the new junctions will added to the old ones.
   STAR will insert the junctions into genome indices on the fly before mapping, which takes 1~2 minutes.
   The on the fly genome indices can be saved (for reuse) with "--sjdbInsertSave All", into _STARgenome directory inside the current run directory.
   Default --sjdbOverhang is now set at 100, and does not have to be specified unless you need to change this value.

   The "all-sample" 2-pass method can be simplified using this on the fly junction insertion option: 
   (i) run the 1st pass for all samples as usual, with or without annotations
   (ii) run 2nd pass for all samples, listing SJ.out.tab files from all samples in --sjdbFileChrStartEnd /path/to/sj1.tab /path/to/sj2.tab ...

2. New option to activate on the fly "per sample" 2-pass method: "--twopassMode Basic".
   Default --twopass1readsN is now -1, i.e. using all reads in the 1st pass.
   2-pass mode can now be used with annotations, which can be included either at the run-time (see #1), or at the genome generation step.
   Annotated junctions will be included in both the 1st and 2nd passes.

3. Included link (submodule) to Brian Haas' STAR-Fusion code for detecting fusion transcript from STAR chimeric output:
   https://github.com/STAR-Fusion/STAR-Fusion

4. Included Gery Vessere's shared memory implementation for POSIX and SysV. 
   To compile STAR with POSIX shared memory, use `make POSIXSHARED`

5. New option "--chimOutType WithinBAM" to include chimeric alignments together with normal alignments in the main (sorted or unsorted) BAM file(s).
   Formatting of chimeric alignments follows the latest SAM/BAM specifications. Thanks to Felix Schlesinger for thorough testing of this option.

6. New option "--quantTranscriptomeBan Singleend" allows insertions, deletions ans soft-clips in the transcriptomic alignments, which can be used by some expression quantification software (e.g. eXpress). 
   
7. New option "--alignEndsTypeExtension Extend5pOfRead1" to enforce full extension of the 5p of the read1, while all other ends undergo local alignment and may be soft-clipped.
