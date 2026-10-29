# Test targets

The alignment path exercised by `--runMode alignReads`, the functions it calls,
and the options they read, grouped into four areas. A row is ticked (`[x]`)
when tests check the row's behaviour against simulation truth or a
hand-worked value, `[~]` when they cover part of it (some modes, or a check
that a broken implementation could pass), `[ ]` when nothing does. A test
that only checks that STAR exits 0 does not count. "Kept out" marks
behaviour where STAR 2.7.11b contradicts its documentation: no suite test
asserts it, and the discrepancy is reported separately.

Areas: **A** input, seeds, stitching, scoring, filters, multimappers;
**B** SAM/BAM output, attributes, ordering, sorting, unmapped output;
**C** junctions, sjdb/2-pass, GeneCounts, TranscriptomeSAM, Log.final.out;
**D** STARsolo, chimeric.

Out of scope: genomeGenerate internals beyond fixtures, WASP
(`ReadAlign_waspMap`), variant/diploid genomes (`*_transformGenome`,
`Transcript_variation*`, `Variation`), `--inputBAMfile` modes
(`signalFromBAM`, `bamRemoveDuplicates`), wiggle output, liftOver, STARlong
(`ReadAlign_stitchWindowSeeds`, `COMPILE_FOR_LONG_READS`), the splice graph
(`SpliceGraph*`, `*SpliceGraph*`), `--genomeFastaFiles` at the mapping step
(`Genome_insertSequences`, `insertSeqSA`, `--outSAMfilter`), shared memory.

Call path per read: `ReadAlignChunk::processChunks` (reads a chunk of input)
-> `ReadAlignChunk::mapChunk` -> `ReadAlign::oneRead` -> `readLoad` (+
`ClipMate::clip`) -> `mapOneRead` (seeds: `maxMappableLength2strands`,
`storeAligns`) -> `stitchPieces` (`createExtendWindowsWithAlign`,
`assignAlignToWindow`, `stitchWindowAligns` -> `stitchAlignToTranscript`,
`extendAlign`) -> `peOverlapMergeMap` -> `multMapSelect` -> `mappedFilter` ->
`chimericDetection` -> `outputAlignments` (`recordSJ`, `quantTranscriptome`,
`soloRead->record`, `writeSAM` -> `outputTranscriptSAM` / `alignBAM`,
`outReadsUnmapped`). After all reads: `outputSJ`, `bamSortByCoordinate`,
`Quantifications`/`Transcriptome::quantsOutput`, `Solo::processAndOutput`,
`Stats::reportFinal`. `twoPassRunPass1` and `sjdbInsertJunctions` run before
the main pass.

Tests are named `<file>::<test>` with these file abbreviations: `aln`
test_alignment.py, `sam` test_sam_output.py, `jq` test_junctions_quant.py,
`sc` test_solo_chimeric.py, `a_input` test_a_input.py, `a_ext`
test_a_extension.py, `a_filt` test_a_filters.py, `a_mm` test_a_multimap.py,
`a_pe` test_a_pe_overlap.py, `a_score` test_a_scoring.py, `b_attr`
test_b_attributes.py, `b_flags` test_b_flags_unmapped.py, `b_sort`
test_b_sorting.py, `c_sj` test_c_sj_out.py, `c_sjdb` test_c_sjdb_twopass.py,
`c_quant` test_c_quant.py, `c_log` test_c_log_final.py, `d_chim`
test_d_chimeric.py, `d_bc` test_d_solo_barcodes.py, `d_umi`
test_d_solo_umi.py, `d_mm` test_d_solo_multimap.py, `d_feat`
test_d_solo_features.py, `d_cf` test_d_solo_cellfilter.py, `d_out`
test_d_solo_outputs.py, `d_bam` test_d_solo_bam.py.

## Functions

| Done | Area | File | Function(s) | What to check | Covered by |
|---|---|---|---|---|---|
| [x] | A | Parameters_readFilesInit.cpp | `Parameters::readFilesInit` | comma-separated `--readFilesIn` (file index, read groups), `--readFilesPrefix`, `--readFilesManifest` | a_input::test_read_files_manifest_and_prefix, b_attr::test_read_groups, sam::test_read_groups_from_comma_separated_inputs |
| [x] | A | Parameters_openReadsFiles.cpp, Parameters_closeReadsFiles.cpp | `openReadsFiles`, `closeReadsFiles` | `--readFilesCommand` (gzip -cd, cat), one process per file | aln::test_read_input_formats[gz], a_input::test_paired_end_name_from_mate1, a_input::test_read_files_manifest_and_prefix (`cat` not run) |
| [x] | A | ReadAlignChunk_processChunks.cpp | `processChunks`, `fastqReadOneLine`, `removeStringEndControl` | chunk filling, FASTA vs FASTQ, multi-line FASTA, name trimming at space and `--readNameSeparator`, `--outSAMreadID Number`, `--readMapNumber`, CRLF line ends, `--outSAMorder PairedKeepInputOrder` | a_input::test_fastq_crlf_and_multiline_fasta, a_input::test_read_name_trimming, b_flags::test_sam_order_keep_input_order (CRLF FASTQ without a name comment kept out) |
| [~] | A | readLoad.cpp | `readLoad` | sequence/quality loading, length mismatch error, `--readQualityScoreBase`, `--outQSconversionAdd`, lower-case and N bases | a_input::test_quality_conversion_add, a_input::test_quality_length_mismatch_is_fatal, a_input::test_fastq_crlf_and_multiline_fasta (open: `--readQualityScoreBase`, read only by SoloFeature_statsOutput.cpp) |
| [~] | A | ParametersClip_initialize.cpp, ClipMate_initialize.cpp, ClipMate_clip.cpp, ClipMate_clipChunk.cpp, ClipCR4.cpp | `ParametersClip::initialize*`, `ClipMate::clip`, `clipChunk`, `ClipCR4::*` | `--clip5pNbases`, `--clip3pNbases` (per mate), `--clip3pAdapterSeq` (per mate, polyA), `--clip3pAdapterMMp`, `--clip3pAfterAdapterNbases`, `--clipAdapterType CellRanger4/None`, filters normalised to clipped length | a_input::test_clip_bases_per_mate, a_input::test_clip3p_adapter_per_mate, a_input::test_clip3p_adapter_mismatch_proportion (open: CellRanger4; None and single values for paired reads kept out) |
| [x] | A | SequenceFuns.cpp | `convertNucleotidesToNumbers`, `complementSeqNumbers`, `revComplementNucleotides`, `localSearch`, `qualitySplit` | numeric encoding, reverse complement of SEQ/QUAL, adapter search, splitting reads at N | aln::test_se_reads_align_to_their_origin, a_input::test_seed_split_min, a_input::test_clip3p_adapter_mismatch_proportion |
| [x] | A | ReadAlign_oneRead.cpp | `ReadAlign::oneRead` | mate concatenation, `Lread`, mismatch budget (`--outFilterMismatchNoverReadLmax`) | aln::test_paired_end_reads, a_filt::test_mismatch_budget_paired_end, a_filt::test_score_filter_paired_end_sum_of_mates |
| [~] | A | ReadAlign_mapOneRead.cpp | `ReadAlign::mapOneRead` | seed start points (`--seedSearchStartLmax`, `--seedSearchStartLmaxOverLread`), `--seedSplitMin`, `--seedMapMin`, reads with no seeds (unmapped other) | a_input::test_seed_split_min, a_ext::test_seed_map_min, a_mm::test_seed_search_options_keep_true_alignments (start options only shown not to lose true alignments; seedMapMin boundary kept out) |
| [~] | A | ReadAlign_maxMappableLength2strands.cpp, SuffixArrayFuns.cpp, binarySearch2.cpp | `maxMappableLength2strands`, `compareSeqToGenome`, `maxMappableLength`, `suffixArraySearch1`, `findMultRange`, `compareRefEnds`, `binarySearch2` | exact seed positions on both strands, SAindex lookup, multimapper ranges, `--seedSearchLmax` | aln::test_multimappers_report_all_loci, aln::test_spliced_reads_found_without_annotation, a_mm::test_seed_search_options_keep_true_alignments (`--seedSearchLmax` checked as weakly as the start options) |
| [~] | A | ReadAlign_storeAligns.cpp | `ReadAlign::storeAligns` | `--seedPerReadNmax`, `--seedMultimapNmax`, duplicate seeds | a_mm::test_seed_multimap_nmax (open: `--seedPerReadNmax`; a seed with as many loci as the limit kept out) |
| [~] | A | ReadAlign_createExtendWindowsWithAlign.cpp, ReadAlign_assignAlignToWindow.cpp | window creation and seed assignment | `--winAnchorMultimapNmax`, `--winBinNbits`, `--winAnchorDistNbins`, `--winFlankNbins`, `--seedPerWindowNmax`, `--seedNoneLociPerWindow`, `--alignWindowsPerReadNmax` | a_mm::test_seed_per_window_nmax, a_ext::test_window_bins_bound_intron_length, a_mm::test_primary_is_best_when_written_second (holds the winAnchorMultimapNmax check) (open: winFlankNbins, seedNoneLociPerWindow; alignWindowsPerReadNmax kept out) |
| [~] | A | ReadAlign_stitchPieces.cpp | `ReadAlign::stitchPieces` | windows to transcripts, `--alignTranscriptsPerWindowNmax`, `--alignTranscriptsPerReadNmax` | a_mm::test_transcripts_per_window_nmax, a_filt::test_unmapped_categories (alignTranscriptsPerReadNmax kept out) |
| [x] | A | stitchWindowAligns.cpp, blocksOverlap.cpp | `stitchWindowAligns`, `blocksOverlap` | recursive seed inclusion, mate consistency, genomic-length score, `--alignMatesGapMax`, `--alignEndsProtrude`, `--alignSplicedMateMapLmin(OverLmate)`, `--outFilterIntronMotifs`, `--outFilterIntronStrands`, `--outSAMstrandField` filtering, `--alignSoftClipAtReferenceEnds`, duplicate transcripts | a_ext::test_mates_gap_max, a_ext::test_spliced_mate_map_min, a_filt::test_intron_strand_filter |
| [~] | A | stitchAlignToTranscript.cpp | `stitchAlignToTranscript` | junction placement (`--scoreStitchSJshift`, repeats), motif penalties `--scoreGap*`, `--sjdbScore`, `--alignIntronMin` (deletion vs intron), `--alignIntronMax`, `--alignSJoverhangMin`, `--alignSJDBoverhangMin`, `--alignSJstitchMismatchNmax`, `--alignInsertionFlush`, `--alignEndsType` | a_score::test_stitch_sj_shift_and_mismatch_limit, a_ext::test_intron_min_separates_deletions_from_introns, c_sjdb::test_sjdb_score_breaks_tie_between_copies (open: junctions that can slide over a repeat) |
| [~] | A | stitchGapIndel.cpp | `stitchGapIndel` | indel placement, `--scoreDelOpen/Base`, `--scoreInsOpen/Base` | a_score::test_indel_penalties, aln::test_alignment_score_hand_worked, a_score::test_insertion_flush_right (open: placement of ambiguous indels) |
| [x] | A | extendAlign.cpp | `extendAlign` | soft-clip vs extension through mismatches, `--alignEndsType EndToEnd/Extend5pOfRead1/Extend5pOfReads12`, extension at chromosome ends | a_ext::test_end_mismatches_by_align_ends_type, a_ext::test_mate_5p_extension_by_align_ends_type, a_ext::test_soft_clip_at_reference_ends |
| [x] | A | Transcript.cpp, Transcript_alignScore.cpp | `Transcript::add/reset/extractSpliceJunctions`, `alignScore` | AS, nM | aln::test_alignment_score_hand_worked, b_attr::test_attribute_values |
| [x] | A | ReadAlign_multMapSelect.cpp | `ReadAlign::multMapSelect` | `--outFilterMultimapScoreRange`, `--outMultimapperOrder Random` with `--runRNGseed`, primary choice (`--outSAMprimaryFlag`) | aln::test_multimap_score_range, a_mm::test_random_multimapper_order, b_flags::test_primary_flag_best_score |
| [x] | A | ReadAlign_mappedFilter.cpp | `ReadAlign::mappedFilter` | `--outFilterScoreMin(OverLread)`, `--outFilterMatchNmin(OverLread)`, `--outFilterMismatchNmax`, `--outFilterMismatchNoverLmax`, `--outFilterMultimapNmax`; unmapType 0-3 | a_filt::test_score_and_match_filters, a_filt::test_mismatch_filters, a_filt::test_unmapped_categories (single-end OverLread rounding kept out) |
| [x] | A | ReadAlign_peOverlapMergeMap.cpp | `peOverlapMergeMap`, `peMergeMates`, `peOverlapSEtoPE` | `--peOverlapNbasesMin`, `--peOverlapMMp` | a_pe::test_overlap_merging_rescues_short_overhangs, a_pe::test_overlap_mismatch_proportion, d_chim::test_pe_overlap_merged_chimera |
| [x] | A/C | ReadAlign_outputAlignments.cpp | `outFilterBySJout` | `--outFilterType BySJout` two-stage filtering | c_sj::test_filter_type_bysjout |
| [x] | A/B | mapThreadsSpawn.cpp, ThreadControl.cpp, ReadAlignChunk.cpp | thread spawn, `chunkFstreamCat`, `chunkFilesCat` | `--runThreadN` 1 vs >1 give the same outputs, with several input chunks | sam::test_thread_count_does_not_change_output, a_mm::test_thread_count_area_a_options, c_sj::test_outputs_do_not_depend_on_input_chunks, d_mm::test_solo_outputs_do_not_depend_on_input_chunks |
| [x] | B | ReadAlign_outputAlignments.cpp | `outputAlignments`, `writeSAM` | `--outSAMmultNmax`, `--outSAMunmapped Within [KeepPairs]`, mate-unmapped records, order of multimapper records | b_flags::test_unmapped_reads_in_alignment_files, b_flags::test_mult_nmax_writes_top_score_first, a_mm::test_old_multimapper_order (KeepPairs in sorted BAM kept out) |
| [x] | B | ReadAlign_outputTranscriptSAM.cpp | `outputTranscriptSAM` | SAM text: FLAG, MAPQ (`--outSAMmapqUnique`), RNEXT/PNEXT/TLEN, SEQ/QUAL (`--outSAMmode NoQS`), `--outSAMflagOR/AND`, `--outSAMattrIHstart`, every attribute, uT for unmapped reads | b_attr::test_attribute_values[sam], b_attr::test_mapq_ladder_and_hi_start[sam], b_flags::test_flag_and_or[sam] (NoQS and flagOR/AND on unmapped SAM records kept out) |
| [x] | B | ReadAlign_alignBAM.cpp | `alignBAM`, `samAttrNM_MD` | the same fields in the BAM writer (separate code), `--outSAMtlen 2`, SEQ encoding, integer tag types | b_attr::test_bam_integer_tag_types, b_sort::test_bam_seq_odd_length_and_lower_case, b_sort::test_sam_and_bam_writers_agree, b_flags::test_tlen_method2_for_non_protruding_pairs |
| [x] | B | ReadAlign_calcCIGAR.cpp | `calcCIGAR` | CIGAR of each mate, soft clips | aln::test_paired_end_reads, a_input::test_clip_bases_per_mate, a_score::test_insertion_flush_right |
| [x] | B | BAMfunctions.cpp | `bamAttrArrayWrite*`, `outBAMwriteHeader`, `bam_cigarString`, `reg2bin` | BAM tag types, header, bin field | b_attr::test_bam_integer_tag_types, b_sort::test_sam_and_bam_writers_agree, b_sort::test_header_lines |
| [x] | B | BAMoutput.cpp | `unsortedOneAlign`, `unsortedFlush`, `coordOneAlign`, `coordBins`, `coordFlush`, `coordUnmappedPrepareBySJout` | `--outBAMcompression`, `--outBAMsortingBinsN`, `--limitBAMsortRAM`, `--limitIObufferSize` | b_sort::test_bam_compression_levels, b_sort::test_sorted_bam_invariant_to_sorting_options, c_sj::test_filter_type_bysjout[sorted] (compression level 10 kept out) |
| [x] | B | bamSortByCoordinate.cpp, BAMbinSortByCoordinate.cpp, BAMbinSortUnmapped.cpp, bam_cat.c | sorting and concatenation | order by (ref, pos), ties, unmapped reads last, `--outBAMsortingThreadN` | b_sort::test_sorted_bam_order_ties_and_unmapped, b_sort::test_sorted_bam_invariant_to_sorting_options |
| [x] | B | samHeaders.cpp, Parameters_readSAMheader.cpp | `samHeaders` | @HD/@SQ/@PG/@CO/@RG, `--outSAMheaderHD`, `--outSAMheaderPG`, `--outSAMheaderCommentFile`, `--outSAMattrRGline` | b_sort::test_header_lines, b_attr::test_read_groups, a_input::test_read_files_manifest_and_prefix |
| [~] | B | Parameters_samAttributes.cpp | `Parameters::samAttributes` | presets None/Standard/All, user order, RG added when an RG line is given, BAM-only tags rejected for SAM | b_attr::test_attribute_presets, b_attr::test_read_groups, sam::test_attribute_order_and_junction_tags (All and BAM-only tags with SAM output kept out) |
| [x] | B | ReadAlign_outputAlignments.cpp | `outReadsUnmapped` | `--outReadsUnmapped Fastx`, mate status 00/10/01 | sam::test_unmapped_fastx_output, b_flags::test_unmapped_fastx_pairs_and_names |
| [~] | B | STAR.cpp, Parameters.cpp | `--outStd SAM/BAM_Unsorted/BAM_SortedByCoordinate/BAM_Quant` | alignments on stdout (run.stdout_bytes) | b_sort::test_out_std (open: BAM_Quant) |
| [x] | C | ReadAlign_outputTranscriptSJ.cpp, ReadAlign_outputAlignments.cpp | `outputTranscriptSJ`, `recordSJ` | per-read junction records, junctions crossed by both mates counted once, `--outSJfilterReads Unique`, `--limitOutSJoneRead` | c_sj::test_sj_out_columns_every_motif, c_sj::test_sj_counts_read_pairs_once, c_sj::test_sj_filter_reads_unique, c_sj::test_sj_collapse_with_small_buffers |
| [x] | C | OutSJ.cpp | `OutSJ::collapseSJ`, `Junction::collapseOneSJ`, `compareSJ`, `outputStream` | unique/multi counts, max overhang, annotation flag, junctions sharing a donor or acceptor | c_sj::test_sj_out_columns_every_motif, c_sj::test_sj_filters, c_sj::test_sj_collapse_with_small_buffers |
| [x] | C | outputSJ.cpp | `outputSJ` | merge across threads, `--outSJfilterOverhangMin`, `--outSJfilterCountUniqueMin`, `--outSJfilterCountTotalMin`, `--outSJfilterDistToOtherSJmin`, `--outSJfilterIntronMaxVsReadN`, annotated junctions exempt, `--outSJtype None`, `--limitOutSJcollapsed` | c_sj::test_sj_filters, c_sj::test_sj_type_none, c_sj::test_outputs_do_not_depend_on_input_chunks |
| [~] | C | GTF.cpp, GTF_transcriptGeneSJ.cpp | `GTF::GTF`, `transcriptGeneSJ` | GTF loading at the mapping step: `--sjdbGTFchrPrefix`, `--sjdbGTFfeatureExon`, `--sjdbGTFtagExonParent*`; transcripts nested in another transcript | c_sjdb::test_sjdb_gtf_chr_prefix_feature_and_tags, c_quant::test_nested_annotation (open: `--sjdbGTFtagExonParentGeneName`, `--sjdbGTFtagExonParentGeneType`) |
| [x] | C | sjdbInsertJunctions.cpp, sjdbPrepare.cpp, sjdbBuildIndex.cpp, sjdbLoadFromFiles.cpp, sjdbLoadFromStream.cpp | on-the-fly insertion | `--sjdbGTFfile`, `--sjdbFileChrStartEnd`, `--sjdbOverhang`, `--sjdbInsertSave All`, `--limitSjdbInsertNsj` (the inserted suffix array is checked only through alignments) | c_sjdb::test_sjdb_file_chr_start_end_on_the_fly, c_sjdb::test_sjdb_overhang_and_insert_save_all, c_sjdb::test_limit_sjdb_insert_nsj |
| [x] | C | sjAlignSplit.cpp | `sjAlignSplit` | alignments to inserted junction sequences turned back into spliced genomic alignments | c_sjdb::* (4-base-overhang reads), jq::test_sjdb_inserted_on_the_fly_matches_prebuilt_index |
| [x] | C | twoPassRunPass1.cpp | `twoPassRunPass1` | `--twopassMode Basic`, `--twopass1readsN`, junctions annotated in pass 2, filtered pass-1 junctions not inserted, quantification in pass 2 only | c_sjdb::test_twopass1readsN, c_sjdb::test_twopass_inserts_filtered_junctions_and_counts_once |
| [x] | C | Transcriptome.cpp, Transcriptome_geneCountsAddAlign.cpp, Quantifications.cpp | `Transcriptome::Transcriptome`, `quantsAllocate`, `quantsOutput`, `geneCountsAddAlign`, `Quantifications::addQuants` | ReadsPerGene.out.tab: three strandedness columns, N_* rows, ambiguous and no-feature reads, exon edges, paired-end reads, genes inside another gene's exon | c_quant::test_gene_counts_overlap_boundaries, c_quant::test_gene_counts_and_transcriptome_paired_end, c_quant::test_nested_annotation |
| [x] | C | ReadAlign_quantTranscriptome.cpp, Transcriptome_quantAlign.cpp | `quantTranscriptome`, `quantAlign`, `alignToTranscript` | transcriptomic coordinates, NH over transcripts, nested transcripts, `--quantTranscriptomeSAMoutput` (BanSingleEnd_BanIndels_ExtendSoftclip, BanSingleEnd, BanSingleEnd_ExtendSoftclip), `--quantTranscriptomeBAMcompression` | c_quant::test_transcriptome_sam_output_modes, c_quant::test_nested_annotation, c_quant::test_transcriptome_bam_compression |
| [x] | C | Stats.cpp | `Stats::addStats`, `transcriptStats`, `reportFinal` | every Log.final.out field (splice types, mismatch/indel rates, read lengths, percentages), merged over input chunks | c_log::test_log_final_every_field, c_log::test_log_final_too_many_mismatches, c_sj::test_outputs_do_not_depend_on_input_chunks |
| [x] | D | ParametersSolo.cpp | `ParametersSolo::initialize` | `--soloType CB_UMI_Simple`, CB/UMI geometry, whitelist loading, `--soloFeatures`, `--soloStrand`, `--soloOutFileNames` | d_bc::test_barcode_geometry_umi_before_cb, d_feat::test_gene_features_per_strand, d_out::test_out_file_names_and_field3 |
| [x] | D | SoloReadBarcode.cpp, SoloReadBarcode_getCBandUMI.cpp, SoloBarcode.cpp | `getCBandUMI`, `matchCBtoWL`, `convertCheckUMI`, `addStats`, `statsOut` | `--soloCBmatchWLtype` (Exact, 1MM, 1MM_multi, pseudocounts, Nbase), N in CB/UMI, homopolymer UMIs, `--soloBarcodeReadLength`, `--soloBarcodeMate`, Barcodes.stats | d_bc::test_cb_whitelist_matching, d_bc::test_barcode_read_length_check, d_out::test_barcode_inside_cdna_mate |
| [x] | D | SoloRead.cpp, SoloRead_record.cpp, SoloReadFeature.cpp, SoloReadFeature_record.cpp, SoloReadFeature_inputRecords.cpp | per-read feature records | Features.stats categories | d_feat::test_gene_features_per_strand, d_bc::test_cb_whitelist_matching, d_mm::test_unique_matrix_ignores_multigene_reads |
| [x] | D | Transcriptome_classifyAlign.cpp, Transcriptome_geneFullAlignOverlap.cpp, Transcriptome_geneFullAlignOverlap_ExonOverIntron.cpp, Transcriptome_alignExonOverlap.cpp | feature assignment | Gene, GeneFull, GeneFull_ExonOverIntron, GeneFull_Ex50pAS, SJ, Velocyto features; intronic, antisense, nested and overlapping genes; `--soloStrand Reverse/Unstranded` | d_feat::test_gene_features_per_strand, d_feat::test_velocyto_spliced_unspliced_ambiguous, d_feat::test_sj_feature_counts_junctions |
| [x] | D | SoloFeature_processRecords.cpp, SoloFeature_sumThreads.cpp, SoloFeature_countCBgeneUMI.cpp, SoloFeature_redistributeReadsByCB.cpp | record collection and per-cell counting | thread merge over several input chunks | d_mm::test_solo_outputs_do_not_depend_on_input_chunks |
| [x] | D | SoloFeature_collapseUMIall.cpp, SoloFeature_collapseUMI_Graph.cpp | UMI collapsing | `--soloUMIdedup` (1MM_All, 1MM_Directional, 1MM_Directional_UMItools, Exact, NoDedup, 1MM_CR), `--soloUMIfiltering`, `--soloMultiMappers` (MultiGeneUMI_All: documented behaviour not met, kept out) | d_umi::test_umi_dedup_methods_in_one_run, d_umi::test_umi_filtering_multigene, d_mm::test_unique_and_mult_matrices |
| [~] | D | SoloFeature_cellFiltering.cpp, SoloFeature_emptyDrops_CR.cpp | cell calling | `--soloCellFilter CellRanger2.2/TopCells/EmptyDrops_CR/None`, filtered/ matrices | d_cf::* (TopCells with fewer cells than detected: documented behaviour not met, kept out; EmptyDrops_CR extra cells pinned, not derived) |
| [x] | D | SoloFeature_outputResults.cpp, SoloFeature_statsOutput.cpp, Solo.cpp | `outputResults`, `statsOutput`, `processAndOutput` | matrix.mtx, barcodes.tsv, features.tsv (`--soloOutFormatFeaturesGeneField3`), Features.stats, Summary.csv, UMIperCellSorted.txt, `--soloCellReadStats` | d_cf::test_filtered_summary_statistics, d_out::test_out_file_names_and_field3, d_out::test_cell_read_stats |
| [x] | D | SoloFeature_addBAMtags.cpp | `addBAMtags` | CB/UB/GX/GN/CR/UR tags in sorted BAM | d_bam::test_solo_bam_tags, d_mm::test_solo_outputs_do_not_depend_on_input_chunks |
| [x] | D | ReadAlign_chimericDetection.cpp, ReadAlign_chimericDetectionOld.cpp, ReadAlign_chimericDetectionOldOutput.cpp | old chimeric scheme (`--chimMultimapNmax 0`) | Chimeric.out.junction columns, WithinBAM HardClip/SoftClip, SA tags, `--chimSegmentMin`, `--chimScore*`, `--chimJunctionOverhangMin`, `--chimSegmentReadGapMax`, `--chimFilter`, `--chimMainSegmentMultNmax`, PE encompassing junctions (type -1), spliced segments | d_chim::* , sc::test_chimeric_within_bam |
| [x] | D | ChimericDetection_chimericDetectionMult.cpp, ChimericAlign*.cpp, ChimericSegment.cpp, ReadAlign_chimericDetectionPEmerged.cpp | multimapping chimeras (`--chimMultimapNmax > 0`) | extra columns, `--chimMultimapScoreRange`, `--chimOutJunctionFormat 1`, segments in one window, merged mates | d_chim::test_new_scheme_columns, d_chim::test_multimap_score_range, d_chim::test_intrachromosomal_chimeras, d_chim::test_pe_overlap_merged_chimera |

## Options

Defaults are from source/parametersDefault. "Covered by" lists tests that set
the option to a non-default value or check its default effect by name.

| Done | Area | Option | Default | Covered by |
|---|---|---|---|---|
| [x] | A | runThreadN | 1 | sam::test_thread_count_does_not_change_output, a_mm::test_thread_count_area_a_options (every run sets it) |
| [x] | A | readFilesIn | Read1 Read2 | sam::test_read_groups_from_comma_separated_inputs, b_attr::test_read_groups, b_flags::test_read_id_number |
| [x] | A | readFilesCommand | - | aln::test_read_input_formats[gz], a_input::test_paired_end_name_from_mate1 |
| [x] | A | readFilesPrefix | - | a_input::test_read_files_manifest_and_prefix |
| [x] | A | readFilesManifest | - | a_input::test_read_files_manifest_and_prefix |
| [x] | A | readMapNumber | -1 | a_input::test_read_map_number |
| [ ] | A | readMatesLengthsIn | NotEqual | (parsed but read nowhere else in source/) |
| [x] | A | readNameSeparator | / | a_input::test_read_name_trimming, aln::test_read_input_formats[slash] |
| [ ] | A | readQualityScoreBase | 33 | (read only by SoloFeature_statsOutput.cpp) |
| [ ] | A | clipAdapterType | Hamming | (CellRanger4 untested; None kept out) |
| [x] | A | clip3pNbases | 0 | aln::test_clip_fixed_bases, a_input::test_clip_bases_per_mate |
| [x] | A | clip5pNbases | 0 | aln::test_clip_fixed_bases, a_input::test_clip_bases_per_mate |
| [x] | A | clip3pAdapterSeq | - | aln::test_clip_3p_adapter_and_polya, a_input::test_clip3p_adapter_per_mate |
| [x] | A | clip3pAdapterMMp | 0.1 | a_input::test_clip3p_adapter_mismatch_proportion |
| [x] | A | clip3pAfterAdapterNbases | 0 | a_input::test_clip3p_after_adapter_bases |
| [x] | A | outFilterType | Normal | c_sj::test_filter_type_bysjout |
| [x] | A | outFilterMultimapScoreRange | 1 | aln::test_multimap_score_range, a_mm::test_primary_is_best_scoring_alignment |
| [x] | A | outFilterMultimapNmax | 10 | aln::test_too_many_loci, a_filt::test_unmapped_categories |
| [x] | A | outFilterMismatchNmax | 10 | a_filt::test_mismatch_filters, a_filt::test_unmapped_categories |
| [x] | A | outFilterMismatchNoverLmax | 0.3 | a_filt::test_mismatch_filters |
| [x] | A | outFilterMismatchNoverReadLmax | 1.0 | a_filt::test_mismatch_filters, a_filt::test_mismatch_budget_paired_end |
| [x] | A | outFilterScoreMin | 0 | a_filt::test_score_and_match_filters |
| [x] | A | outFilterScoreMinOverLread | 0.66 | a_filt::test_score_and_match_filters, a_filt::test_score_filter_paired_end_sum_of_mates (single-end rounding kept out) |
| [x] | A | outFilterMatchNmin | 0 | a_filt::test_score_and_match_filters |
| [x] | A | outFilterMatchNminOverLread | 0.66 | a_filt::test_score_and_match_filters |
| [x] | A | outFilterIntronMotifs | None | a_filt::test_intron_motif_filter |
| [x] | A | outFilterIntronStrands | RemoveInconsistentStrands | a_filt::test_intron_strand_filter |
| [x] | A | scoreGap, scoreGapNoncan, scoreGapGCAG, scoreGapATAC | 0 -8 -4 -8 | a_score::test_junction_motif_penalties, a_score::test_stitch_mismatch_limit_atac |
| [x] | A | scoreGenomicLengthLog2scale | -0.25 | a_score::test_genomic_length_log2scale |
| [x] | A | scoreDelOpen, scoreDelBase, scoreInsOpen, scoreInsBase | -2 -2 -2 -2 | a_score::test_indel_penalties, a_ext::test_intron_min_separates_deletions_from_introns |
| [x] | A | scoreStitchSJshift | 1 | a_score::test_stitch_sj_shift_and_mismatch_limit |
| [~] | A | seedSearchStartLmax, seedSearchStartLmaxOverLread, seedSearchLmax | 50 1.0 0 | a_mm::test_seed_search_options_keep_true_alignments (only shows true alignments survive) |
| [~] | A | seedMultimapNmax, seedPerReadNmax, seedPerWindowNmax, seedNoneLociPerWindow | 10000 1000 50 10 | a_mm::test_seed_multimap_nmax, a_mm::test_seed_per_window_nmax (open: seedPerReadNmax, seedNoneLociPerWindow) |
| [x] | A | seedSplitMin, seedMapMin | 12 5 | a_input::test_seed_split_min, a_ext::test_seed_map_min (a piece of exactly seedMapMin bases kept out) |
| [x] | A | alignIntronMin, alignIntronMax | 21 0 | a_ext::test_intron_min_separates_deletions_from_introns, a_ext::test_intron_max |
| [x] | A | alignMatesGapMax | 0 | a_ext::test_mates_gap_max |
| [x] | A | alignSJoverhangMin | 5 | a_ext::test_sj_overhang_min_unannotated, a_ext::test_seed_map_min |
| [x] | A | alignSJDBoverhangMin | 3 | a_ext::test_sjdb_overhang_min, c_sjdb::test_sjdb_file_chr_start_end_on_the_fly |
| [x] | A | alignSJstitchMismatchNmax | 0 -1 0 0 | a_score::test_stitch_sj_shift_and_mismatch_limit, a_score::test_stitch_mismatch_limit_atac |
| [x] | A | alignSplicedMateMapLmin, alignSplicedMateMapLminOverLmate | 0 0.66 | a_ext::test_spliced_mate_map_min |
| [~] | A | alignWindowsPerReadNmax, alignTranscriptsPerWindowNmax, alignTranscriptsPerReadNmax | 10000 100 10000 | a_mm::test_transcripts_per_window_nmax (the other two kept out) |
| [x] | A | alignEndsType | Local | a_ext::test_end_mismatches_by_align_ends_type, a_ext::test_mate_5p_extension_by_align_ends_type, c_log::test_log_final_too_many_mismatches |
| [~] | A | alignEndsProtrude | 0 ConcordantPair | b_flags::test_concordant_protruding_pairs_are_proper, sam::test_tlen_methods_for_protruding_mates (DiscordantPair kept out) |
| [x] | A | alignSoftClipAtReferenceEnds | Yes | a_ext::test_soft_clip_at_reference_ends |
| [x] | A | alignInsertionFlush | None | a_score::test_insertion_flush_right |
| [x] | A | peOverlapNbasesMin, peOverlapMMp | 0 0.01 | a_pe::test_overlap_merging_rescues_short_overhangs, a_pe::test_overlap_mismatch_proportion, d_chim::test_pe_overlap_merged_chimera |
| [~] | A | winAnchorMultimapNmax, winBinNbits, winAnchorDistNbins, winFlankNbins | 50 16 9 4 | a_mm::test_primary_is_best_when_written_second (winAnchorMultimapNmax part), a_ext::test_window_bins_bound_intron_length (open: winFlankNbins) |
| [x] | A | outMultimapperOrder, runRNGseed | Old_2.4 777 | a_mm::test_random_multimapper_order, a_mm::test_old_multimapper_order |
| [x] | B | outSAMtype | SAM | sam::test_sam_unsorted_and_sorted_bam_agree, b_sort::test_sam_and_bam_writers_agree, b_flags::test_no_alignment_output |
| [x] | B | outSAMmode | Full | b_flags::test_sam_mode_noqs, b_flags::test_no_alignment_output (NoQS on unmapped SAM records kept out) |
| [x] | B | outSAMstrandField | None | sam::test_strand_field_intron_motif, b_attr::test_strand_field_intron_motif |
| [x] | B | outSAMattributes | Standard | b_attr::test_attribute_presets, b_attr::test_attribute_values, sam::test_attribute_order_and_junction_tags (All with SAM output kept out) |
| [x] | B | outSAMattrIHstart | 1 | b_attr::test_mapq_ladder_and_hi_start |
| [x] | B | outSAMunmapped | None | b_flags::test_unmapped_reads_in_alignment_files, sam::test_unmapped_within_reports_reason (KeepPairs in sorted BAM kept out) |
| [x] | B | outSAMorder | Paired | b_flags::test_sam_order_keep_input_order |
| [x] | B | outSAMprimaryFlag | OneBestScore | sam::test_primary_flag_all_best_score, b_flags::test_primary_flag_best_score, a_mm::test_primary_is_best_when_written_second |
| [x] | B | outSAMreadID | Standard | b_flags::test_read_id_number, a_input::test_read_name_trimming |
| [x] | B | outSAMmapqUnique | 255 | b_attr::test_mapq_ladder_and_hi_start |
| [x] | B | outSAMflagOR, outSAMflagAND | 0 65535 | b_flags::test_flag_and_or (unmapped SAM records kept out) |
| [x] | B | outSAMattrRGline | - | sam::test_read_groups_from_comma_separated_inputs, b_attr::test_read_groups, b_sort::test_header_lines |
| [x] | B | outSAMheaderHD, outSAMheaderPG, outSAMheaderCommentFile | - | b_sort::test_header_lines |
| [x] | B | outSAMmultNmax | -1 | sam::test_mult_nmax_limits_records_not_nh, b_flags::test_mult_nmax_writes_top_score_first |
| [~] | B | outSAMtlen | 1 | sam::test_tlen_methods_for_protruding_mates, b_flags::test_tlen_method2_for_non_protruding_pairs (2 checked for BAM only; the SAM writer ignoring it kept out) |
| [x] | B | outBAMcompression | 1 | b_sort::test_bam_compression_levels (level 10 kept out) |
| [x] | B | outBAMsortingThreadN, outBAMsortingBinsN, limitBAMsortRAM | 0 50 0 | b_sort::test_sorted_bam_invariant_to_sorting_options |
| [~] | B | outStd | Log | b_sort::test_out_std (open: BAM_Quant) |
| [x] | B | outReadsUnmapped | None | sam::test_unmapped_fastx_output, b_flags::test_unmapped_fastx_pairs_and_names |
| [x] | B | outQSconversionAdd | 0 | a_input::test_quality_conversion_add, b_flags::test_qs_conversion_add |
| [x] | B | outTmpDir, outTmpKeep | - None | b_sort::test_tmp_dir_and_keep |
| [x] | C | outSJtype | Standard | c_sj::test_sj_type_none |
| [x] | C | outSJfilterReads | All | c_sj::test_sj_filter_reads_unique |
| [x] | C | outSJfilterOverhangMin | 30 12 12 12 | c_sj::test_sj_filters[overhang, overhang_off] |
| [x] | C | outSJfilterCountUniqueMin, outSJfilterCountTotalMin | 3 1 1 1, 3 1 1 1 | c_sj::test_sj_filters[counts, counts_off] |
| [x] | C | outSJfilterDistToOtherSJmin | 10 0 5 10 | c_sj::test_sj_filters[dist] |
| [x] | C | outSJfilterIntronMaxVsReadN | 50000 100000 200000 | c_sj::test_sj_filters[intronmax, intronmax_edge] |
| [x] | C | sjdbGTFfile, sjdbOverhang (mapping step) | - 100 | c_sjdb::test_sjdb_overhang_and_insert_save_all |
| [x] | C | sjdbFileChrStartEnd | - | c_sjdb::test_sjdb_file_chr_start_end_on_the_fly |
| [~] | C | sjdbGTFchrPrefix, sjdbGTFfeatureExon, sjdbGTFtagExonParentTranscript, sjdbGTFtagExonParentGene, sjdbGTFtagExonParentGeneName, sjdbGTFtagExonParentGeneType | - exon transcript_id gene_id gene_name gene_type/gene_biotype | c_sjdb::test_sjdb_gtf_chr_prefix_feature_and_tags (GeneName, GeneType open) |
| [x] | C | sjdbScore | 2 | c_sjdb::test_sjdb_score_breaks_tie_between_copies |
| [x] | C | sjdbInsertSave, limitSjdbInsertNsj | Basic 1000000 | c_sjdb::test_sjdb_overhang_and_insert_save_all, c_sjdb::test_limit_sjdb_insert_nsj |
| [x] | C | twopassMode | None | c_sjdb::test_twopass_inserts_filtered_junctions_and_counts_once |
| [x] | C | twopass1readsN | -1 | c_sjdb::test_twopass1readsN |
| [x] | C | quantMode | - | c_quant::* |
| [x] | C | quantTranscriptomeSAMoutput | BanSingleEnd_BanIndels_ExtendSoftclip | c_quant::test_transcriptome_sam_output_modes |
| [~] | C | quantTranscriptomeBAMcompression | 1 | c_quant::test_transcriptome_bam_compression (-2, 0, 1, 9; level 10 does not behave as documented, kept out) |
| [x] | C | limitOutSJoneRead, limitOutSJcollapsed | 1000 1000000 | c_sj::test_sj_collapse_with_small_buffers |
| [x] | D | soloType | None | d_*::* (CB_UMI_Simple) |
| [x] | D | soloCBwhitelist, soloCBstart, soloCBlen, soloUMIstart, soloUMIlen | - 1 16 17 10 | d_bc::test_barcode_geometry_umi_before_cb, d_bc::test_no_whitelist_counts_every_barcode |
| [x] | D | soloBarcodeReadLength, soloBarcodeMate | 1 0 | d_bc::test_barcode_read_length_check, d_out::test_barcode_inside_cdna_mate |
| [x] | D | soloCBmatchWLtype | 1MM_multi | d_bc::test_cb_whitelist_matching |
| [x] | D | soloStrand | Forward | d_feat::test_gene_features_per_strand |
| [x] | D | soloFeatures | Gene | d_feat::* |
| [x] | D | soloMultiMappers | Unique | d_mm::test_unique_and_mult_matrices |
| [x] | D | soloUMIdedup | 1MM_All | d_umi::test_umi_dedup_methods_in_one_run |
| [~] | D | soloUMIfiltering | - | d_umi::test_umi_filtering_multigene (MultiGeneUMI_All, and MultiGeneUMI for a UMI with one read in each of two genes, do not behave as documented, kept out) |
| [x] | D | soloOutFileNames | Solo.out/ features.tsv barcodes.tsv matrix.mtx | d_out::test_out_file_names_and_field3 |
| [~] | D | soloCellFilter | CellRanger2.2 3000 0.99 10 | d_cf::* (TopCells N below the detected cells kept out; EmptyDrops_CR result pinned) |
| [x] | D | soloOutFormatFeaturesGeneField3 | "Gene Expression" | d_out::test_out_file_names_and_field3, d_out::test_field3_custom_string |
| [x] | D | soloCellReadStats | None | d_out::test_cell_read_stats |
| [~] | D | chimSegmentMin | 0 | d_chim::test_segment_length_and_overhang (values 49 and 50 for 50-base segments do not behave as documented, kept out) |
| [x] | D | chimOutType | Junctions | d_chim::test_separate_sam_old, d_chim::test_within_bam_soft_clip, sc::test_chimeric_within_bam |
| [~] | D | chimScoreMin, chimScoreDropMax, chimScoreSeparation, chimScoreJunctionNonGTAG | 0 20 10 -1 | d_chim::test_score_min_and_drop_max, d_chim::test_score_separation, d_chim::test_junction_non_gtag_penalty (chimScoreSeparation at the boundary kept out) |
| [x] | D | chimJunctionOverhangMin, chimSegmentReadGapMax, chimFilter, chimMainSegmentMultNmax | 20 0 banGenomicN 10 | d_chim::test_segment_length_and_overhang, d_chim::test_segment_read_gap_max, d_chim::test_chim_filter_genomic_n, d_chim::test_main_segment_multimapping |
| [~] | D | chimMultimapNmax, chimMultimapScoreRange, chimNonchimScoreDropMin, chimOutJunctionFormat | 0 1 20 0 | d_chim::test_multimapping_chimera_nmax, d_chim::test_multimap_score_range, d_chim::test_nonchim_score_drop_min, d_chim::test_junction_format_comment_lines (chimNonchimScoreDropMin at the boundary kept out) |

## Counts

Rows (a row may group related functions or options); rows marked A/B or A/C
are counted under A. Done / partly / open after the area reviews; the pilot
counts were A 6/14/2 and 10/7/31, B 2/8/1 and 7/3/11, C 4/6/0 and 4/2/11,
D 1/7/3 and 3/3/12.

| Area | Function rows | done / partly / open | Option rows | done / partly / open |
|---|---|---|---|---|
| A | 22 | 13 / 9 / 0 | 48 | 40 / 5 / 3 |
| B | 11 | 9 / 2 / 0 | 21 | 19 / 2 / 0 |
| C | 10 | 9 / 1 / 0 | 17 | 15 / 2 / 0 |
| D | 11 | 10 / 1 / 0 | 18 | 13 / 5 / 0 |
