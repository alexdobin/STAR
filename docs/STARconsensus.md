STARconsensus: mapping RNA-seq reads to consensus genome.
=========================================================

* Introduced in STAR 2.7.7a (2020/12/28)

* Provide the VCF file with consensus SNVs and InDels at the genome generation stage with ```--genomeTransformVCF Variants.vcf --genomeTransformType Haploid```.
The alternative alleles in this VCF will be inserted to the reference genome to create a "transformed" genome.
Both the genome sequence and transcript/gene annotations are transformed.

* At the mapping stage, the reads will be mapped to the tranformed (consensus) genome.
The quantification in the transformed annotations can be performed with standard ```--quantMode TranscriptomeSAM and/or GeneCounts``` options.
If desired, alignments (SAM/BAM) and spliced junctions (SJ.out.tab) can be transformed back to the original (reference) coordinates with ```--genomeTransformOutput SAM and/or SJ```.
This is useful if downstream processing relies on reference coordinates.
