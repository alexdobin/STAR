# STAR unit tests

A pytest suite for STAR's read-alignment path. Each test simulates reads from
a small synthetic genome whose sequence and annotation are known, runs the
STAR binary on them, and compares STAR's outputs with the simulated truth or
with values worked out by hand. The suite needs Python >= 3.8 and pytest;
samtools and pysam are optional and only used to cross-check the harness's own
BAM decoder.

## Build and run

```bash
make -C source STAR                      # builds source/STAR
python3 -m pytest -q tests               # uses source/STAR
STAR_BIN=/path/to/STAR python3 -m pytest -q tests
STAR_TEST_THREADS=4 python3 -m pytest -q tests   # default --runThreadN for every run
```

Without a binary at `$STAR_BIN` or `source/STAR`, every test is skipped with a
message saying how to build one. The suite has 361 tests; on a 256-core
Linux server under load it took 23 s with `STAR_TEST_THREADS=1` and
48 s with `STAR_TEST_THREADS=4`. Run one area with, for example, `python3 -m pytest -q
tests/test_c_*.py tests/test_junctions_quant.py`.

Keep pytest's base temporary directory on a local disk (the default,
`/tmp/pytest-of-$USER`). With run directories on a network file system a
coordinate-sorted BAM run took 1.1 s instead of 0.06 s.

## Known bugs are kept out of the suite

Every test asserts behaviour that STAR's documentation (source/parametersDefault,
doc/STARmanual.pdf, docs/STARsolo.md, CHANGES.md) states, or that follows from
the simulated truth. Where STAR 2.7.11b contradicts its documentation, the test
that would fail is not part of the suite; such discrepancies are reported as
separate issues. Docstrings near such cases say what the documentation
promises and which neighbouring values are tested. A few tests pin the current
output of behaviour the documentation leaves open; their docstrings say "pins
the current output".

## Threads and input chunks

`run_star` always passes `--runThreadN` explicitly: the `threads=` argument if
given, else `$STAR_TEST_THREADS`, else 1. `threads` must be between 1 and
`MAX_THREADS` (8), and `--runThreadN` in `extra_args` is rejected, so STAR
never sizes itself to the machine. BAM sorting uses `min(6, runThreadN)`
threads by default. Tests that compare thread counts pass `threads=1` and
`threads=4` explicitly; all other tests should leave `threads` unset so the
whole suite can be run at 1 and at 4 threads.

STAR reads its input in chunks of up to `limitIObufferSize[0] / mates - 2 *
651 - 2 * 50000` bytes per mate (Parameters.cpp, IncludeDefine.h), about 30
MB by default, so every test read set fits into one chunk and only one thread
maps reads. Tests of the merging of per-thread results therefore pass
`--limitIObufferSize 110000 50000000` (220000 for two input files), which
gives chunks of 8698 bytes, with `threads=4`, and compare with `threads=1`.

## Files

| File | Covers |
|---|---|
| `star_helpers.py` | the harness: synthetic genome, read simulator, `run_star`, parsers, pure-Python BAM decoder |
| `conftest.py` | session fixtures (binary, genomes, simulator, whitelist) |
| `pytest.ini` | disables the cache provider so runs leave nothing in the source tree |
| `test_harness.py` | the harness itself: genome layout, simulator truth, parsers, the BAM decoder on a hand-encoded BAM and against samtools/pysam, Solo matrix keys |
| `test_alignment.py` | area A: input formats, clipping, seeds and stitching, indels, scores, multimappers, filters, paired-end reads |
| `test_a_input.py` | area A: FASTQ/FASTA parsing, read names, `--readMapNumber`, `--outQSconversionAdd`, `--seedSplitMin`, 5'/3' clipping |
| `test_a_extension.py` | area A: `--alignEndsType`, soft-clipping at reference ends, overhang minimums, intron and mate-gap limits, spliced-mate minimums |
| `test_a_scoring.py` | area A: junction and indel penalties, genomic-length score, `--scoreStitchSJshift`, `--alignSJstitchMismatchNmax`, `--alignInsertionFlush` |
| `test_a_filters.py` | area A: unmapped categories, score/match/mismatch filters, `--outFilterIntronMotifs`, `--outFilterIntronStrands` |
| `test_a_multimap.py` | area A: multimapper order and primary choice, window and seed limits, seed search start options, thread equivalence |
| `test_a_pe_overlap.py` | area A: merging of overlapping mates (`--peOverlapNbasesMin`, `--peOverlapMMp`) |
| `helpers_a.py` | area A helpers: the STAR rules the area A expected values are worked from |
| `test_sam_output.py` | area B: SAM vs BAM vs sorted BAM, unmapped output, attributes, primary flag, `--outSAMmultNmax`, read groups, XS, TLEN, thread equivalence |
| `test_b_attributes.py` | area B: attribute values, presets, `--outSAMattrIHstart`, MAPQ, strand field and read groups in all three writers |
| `test_b_flags_unmapped.py` | area B: FLAG options, unmapped output in the alignment files and in Unmapped.out.mate*, read IDs, record order, `--outSAMtlen 2`, `--outSAMmode NoQS` |
| `test_b_sorting.py` | area B: coordinate sorting, BAM compression, `--outStd`, SAM header options, agreement of the three writers |
| `helpers_b.py` | area B helpers: a second synthetic genome, expected scores, header comparisons, BAM bin computation |
| `test_junctions_quant.py` | area C: SJ.out.tab, on-the-fly sjdb, 2-pass, GeneCounts, TranscriptomeSAM and Log.final.out on small hand-built read sets |
| `test_c_sj_out.py` | area C: SJ.out.tab columns for every motif, all `--outSJfilter*` options, `--outSJfilterReads`, `--outSJtype None`, read pairs counted once, `--outFilterType BySJout`, small SJ buffers, outputs independent of input chunks |
| `test_c_sjdb_twopass.py` | area C: junctions inserted at the mapping step (`--sjdbFileChrStartEnd`, `--sjdbGTF*` options, `--sjdbOverhang`, `--sjdbInsertSave`, `--sjdbScore`, `--limitSjdbInsertNsj`) and 2-pass mapping |
| `test_c_quant.py` | area C: GeneCounts overlap and strand rules, nested genes, TranscriptomeSAM coordinates and output modes, paired-end quantification, transcriptome BAM compression |
| `test_c_log_final.py` | area C: every count field of Log.final.out, too-many-mismatches and chimeric counts |
| `helpers_c.py` | area C helpers: spliced and repeat-segment reads, junction cleanliness, expected SJ.out.tab rows |
| `test_solo_chimeric.py` | area D: STARsolo CB_UMI_Simple matrices and stats, chimeric junction and WithinBAM output |
| `test_d_solo_barcodes.py` | area D: every `--soloCBmatchWLtype`, no whitelist, barcode read length and geometry, gzipped input, Barcodes.stats |
| `test_d_solo_umi.py` | area D: every `--soloUMIdedup` type and `--soloUMIfiltering`, against reference implementations of the documented rules |
| `test_d_solo_multimap.py` | area D: `--soloMultiMappers` Uniform, PropUnique, Rescue and EM; Solo outputs independent of threads and input chunks |
| `test_d_solo_features.py` | area D: Gene, GeneFull, GeneFull_Ex50pAS, GeneFull_ExonOverIntron, SJ and Velocyto features under each `--soloStrand` |
| `test_d_solo_cellfilter.py` | area D: `--soloCellFilter` CellRanger2.2, TopCells, EmptyDrops_CR and None; Summary.csv and UMIperCellSorted.txt |
| `test_d_solo_outputs.py` | area D: `--soloOutFileNames`, features.tsv field 3, `--soloCellReadStats`, barcode inside the cDNA mate |
| `test_d_solo_bam.py` | area D: STARsolo BAM tags (CR CY UR UY CB UB GX GN sS sQ sM) |
| `test_d_chimeric.py` | area D: chimeric detection in both schemes: thresholds, columns, output types, multimapping chimeras, paired-end and merged-mate chimeras, spliced segments, outputs independent of input chunks |
| `helpers_d.py` | area D helpers: 10x read builders, chimeric read builders, reference UMI collapsing, multi-gene distribution and knee filtering |
| `TARGETS.md` | functions and options on the alignment path and the tests covering each |

## The synthetic genome

`make_synthetic_genome()` builds the same genome every time (seed
`GENOME_SEED`): chr1 50 kb, chr2 40 kb, chr3 20 kb of uniform random sequence,
so `--genomeSAindexNbases` is 7 (min(14, log2(110000)/2 - 1)). Annotation
(1-based, inclusive; `GENE_SPECS` in star_helpers.py):

| Gene | Locus | Strand | Notes |
|---|---|---|---|
| G01 | chr1:1001-3800 | + | T01a 3 exons, T01b skips exon 2; GT/AG |
| G02 | chr1:6001-8500 | - | 3 exons; CT/AC on the genome |
| G03 | chr1:11001-12200 | + | GC/AG intron 11401-11800 |
| G04 | chr1:15001-16300 | + | AT/AC intron 15301-15900 |
| G05, G06 | chr1:20001-21400 (+), 21201-22600 (-) | +/- | exons overlap at 21201-21400 |
| G07 | chr1:30001-31000 | - | single exon |
| G08 | chr2:2001-5400 | - | 4 exons; CT/GC intron 2301-3000 |
| G09 | chr2:12001-17600 | + | 5 kb intron |
| G10 | chr2:30001-32400 | + | 3 exons; unannotated isoform T10skip (junction 30401-31800) |
| G11 | chr3:12001-13300 | - | GT/AT intron 12501-13000 |
| G12 | chr3:15001-16000 | + | not in the GTF; non-canonical intron 15201-15800 (T12nc) |

Repeated segments (`synth.segments`): `dup2` (1.5 kb, chr1:26001 + and
chr2:8001 reverse-complemented), `dup3` (600 bp, 3 copies), `dupmm` (1 kb, 2
copies, copy 2 with a substitution every 100 bp at offsets 50, 150, ...),
`rep12` (300 bp, 12 copies, more than `--outFilterMultimapNmax 10`). Every
intron is written so its junction cannot slide (last exon base differs from
the last intron base, first exon base from the first intron base).

## Harness API

Import from the module and use the fixtures:

```python
import star_helpers as sh
from star_helpers import ReadSet, record_truth, run_star

def test_example(genome_sjdb, sim):
	reads = sim.transcript_reads(100)
	run = run_star(genome_sjdb, reads, "--outSAMattributes NH HI AS nM")
	for rec in run.alignments():
		assert record_truth(rec) == reads.by_name[rec.qname].truth()
```

### Fixtures (conftest.py)

| Fixture | Scope | Gives |
|---|---|---|
| `star_bin` | session, autouse | path of the binary; skips the session when missing |
| `star_threads` | session | the default `--runThreadN` for this session |
| `synth_genome` | session | the `SyntheticGenome`, with `fasta_path` and `gtf_path` written |
| `genome_plain` | session | `StarGenome` built without annotations |
| `genome_sjdb` | session | `StarGenome` built with the GTF, `--sjdbOverhang 99` |
| `genome_factory` | session | `make(gtf=True, sjdb_overhang=99, extra_args=(), synth=None) -> StarGenome`, cached per argument set |
| `sim` | function | `ReadSimulator` seeded from the test's module and name |
| `tenx_whitelist` | session | `(barcodes, path)`: 24 16-bp barcodes, pairwise >= 3 mismatches apart |

### Running STAR

- `run_star(genome, reads, extra_args=(), threads=None, outdir=None, name=None, fasta=False, gz=False, name_style="plain", check=True, timeout=None) -> StarRun`
  - `reads`: a `ReadSet`; a list of `ReadSet`s (written as comma-separated files, one read group each); a list of paths passed to `--readFilesIn` as given; or `None` when `extra_args` carries `--readFilesIn`.
  - `extra_args`: list of strings or one shell-like string. Paths in it must be absolute (STAR runs inside the run directory).
  - `gz=True` writes gzipped files and adds `--readFilesCommand gzip -cd`; `fasta=True` writes FASTA; `name_style` is `plain`, `slash` (name/1) or `illumina` (name 1:N:0:1).
  - `check=True` raises `StarError` (an `AssertionError` with the argv, stderr and Log.out tail) on a non-zero exit; use `check=False` to test errors via `run.returncode`.
- `build_genome(synth, outdir, gtf=True, sjdb_overhang=99, extra_args=(), threads=None, sa_index_nbases=None) -> StarGenome` (prefer `genome_factory`).
- `OUTPUT_FORMATS`: `{"sam": [], "bam": [...Unsorted], "sorted": [...SortedByCoordinate]}`. SAM and BAM records are written by different code in STAR, so parametrize tests of SAM fields over it: `@pytest.mark.parametrize("fmt", list(sh.OUTPUT_FORMATS))` and add `sh.OUTPUT_FORMATS[fmt]` to `extra_args`.
- `find_star()`, `default_threads()`, `MAX_THREADS`.

`StarRun`: `outdir`, `argv`, `returncode`, `stdout`, `stdout_bytes`,
`stderr`, `truth` (name -> SimRead/SimPair over all input sets), and

- `alignments(name=None) -> AlignmentFile` (Aligned.out.sam, Aligned.out.bam or Aligned.sortedByCoord.out.bam, or the named file); `records(name=None)`; `transcriptome()` (Aligned.toTranscriptome.out.bam)
- `sj() -> [SJRecord]`, `gene_counts() -> {name: (unstranded, strand1, strand2)}`, `log_final() -> {field: value}`, `log_out() -> str`
- `unmapped(mate=1) -> [FastxRecord]`, `chimeric() -> [dict]`
- `solo(feature="Gene", matrix="raw", mtx=None) -> SoloMatrix`, `solo_features_stats(feature)`, `solo_summary(feature)`, `solo_barcodes_stats()`, `solo_dir(feature)`, `solo_out_names()`; all follow `--soloOutFileNames` when the run sets it
- `path(name)`, `exists(name)`

### Simulating reads

`ReadSimulator(synth, seed=1, read_len=100, const_qual=None)` (use the `sim`
fixture). Offsets into reads count from the leftmost aligned base in genome
orientation; transcript offsets count 5' to 3' along the transcript.

- `transcript_read(tx_id, offset, length=None, sense=True, mismatches=(), head="", tail="", n_offsets=(), name=None) -> SimRead`
- `genomic_read(contig, pos, strand="+", length=None, mismatches=(), insertion=None, deletion=None, head="", tail="", n_offsets=(), name=None) -> SimRead`; `insertion`/`deletion` are `(read offset, n)`
- `segment_read(segment, offset, length=None, strand="+", copy=0) -> SimRead` with `extra["loci"]` (all copies) and `extra["n_loci"]`
- `random_read(length=None) -> SimRead` (should not align; `contig` is None)
- `chimeric_read((contig, pos, strand, n), (contig, pos, strand, n)) -> SimRead` with `extra["junction"]` (expected Chimeric.out.junction columns), `extra["repeat"]`, `extra["shift_motifs"]`; `clean_chimeric_read(strands=("+", "+"), length=(50, 50))` picks pieces whose breakpoint cannot move
- `transcript_pair(tx_id, offset, insert, length=None, forward=True, mismatches1=(), mismatches2=(), adapter=ILLUMINA_ADAPTER) -> SimPair`; `genomic_pair(contig, pos, insert, strand="+", ...)`; inserts shorter than the read run into the adapter
- batches returning a `ReadSet`: `transcript_reads(n, tx_ids=None, length=None, min_block=25, sense=None, spliced=None)`, `transcript_pairs(n, tx_ids=None, length=None, insert=(200, 300), min_block=25)`, `tenx_reads(whitelist, cells, cdna_len=90, umi_len=12, reads_per_umi=1, three_prime=300, cb_errors=None)` (`cells` = `{cb: {gene_id: n_umis}}`)
- helpers: `unique_locus(length=None, avoid_genes=False) -> (contig, pos)` outside repeated segments, `deletion_is_unambiguous(contig, start, n)`, `valid_offsets(tx_id, length, min_block, region)`
- module level: `make_whitelist(n, seed=7)`, `expected_solo_counts(readset) -> {(gene_id, cb): n_umis}`, `expected_nm_md(synth, contig, pos, cigar, sam_seq) -> (NM, MD)`

`SimRead`: `name`, `seq`, `qual` (as sequenced), `contig`, `strand`, `pos`,
`cigar` (truth in SAM terms: S for clipped/adapter bases, N for introns),
`mismatches` (genomic positions), `kind`, `tx`, `gene`, `mate`, `extra`;
`truth() -> (contig, pos, cigar, strand)`, `sam_seq()`, `sam_qual()`, `end`,
`introns`. `SimPair`: `name`, `r1`, `r2` (STAR's file order; for 10x reads
`r1` is the cDNA read, `r2` the barcode read), `insert`, `extra`.
`ReadSet(records, kind=None)`: iterable, `by_name`, `paired`, `mates(i)`,
`+` concatenates, `write(directory, ...)`.

### Parsers and records

- `read_alignments(path)`, `parse_sam(path_or_text)`, `read_bam(path_or_bytes)` (pure-Python BGZF/BAM decoder) -> `AlignmentFile` with `header`, `references`, `records`, `by_name()`, `primary()`, `mapped()`, `header_lines("@RG")`
- `SamRecord`: `qname flag rname pos mapq cigar rnext pnext tlen seq qual`, `bin` and `bam_types` (the BAM bin field and the tag type codes as stored, c/C/s/S/i/I...; None for SAM records), `tags` (ordered; B arrays are lists), `tag_types` (BAM integer types reported as `i`), `tag(name)`, `strand`, `mate`, `end`, `introns`, `is_unmapped/is_reverse/is_secondary/is_supplementary/is_read1/is_read2/is_paired/is_proper_pair`, `to_line()`
- `record_truth(rec) -> (rname, pos, cigar, strand)`; `record_multiset(records, fields=None, bam_types=False)` and `first_line_diff(a, b)` for order-free comparisons (`bam_types=True` also compares the stored integer tag types of two BAM outputs)
- `parse_sj_out`, `parse_reads_per_gene`, `parse_log_final`, `parse_fastx`, `parse_chimeric_junction` (`CHIMERIC_FIELDS`), `read_solo_matrix` (`SoloMatrix.counts[(key, cb)]` with key the gene ID, or `(contig, start, end)` for the SJ feature, whose features.tsv is SJ.out.tab; `keys`, `features`, `barcodes`, `shape`), `parse_mtx`, `parse_stats`, `parse_summary_csv`
- CIGAR helpers: `cigar_ops`, `ops_to_cigar`, `cigar_ref_len`, `cigar_query_len`, `cigar_introns(pos, cigar)`, `cigar_blocks(pos, cigar)`, `cigar_softclips`; sequence helpers `revcomp`, `hamming`, `motif_code`, `MOTIF_STRAND`

## Writing tests

- Expected values come from the simulation truth or are worked by hand in the docstring. A test that only pins STAR's current output says "pins the current output" in its docstring; keep those few. Checking only that STAR exits 0 is not coverage.
- When STAR's behaviour contradicts its documentation (source/parametersDefault, doc/STARmanual.pdf, docs/STARsolo.md, CHANGES.md), leave the test out of the suite and report the discrepancy instead, so the suite only asserts documented behaviour.
- Use `sim` for reads so they do not change with test order (a module-scoped fixture seeds its own simulator from its name, never a simulator shared with tests), and the session genomes; build another genome only with `genome_factory`. Keep each test to a few STAR runs on a few hundred reads; share a module-scoped run between tests that only read its outputs.
- Shuffle read sets whose reads would otherwise be grouped (all reads of one UMI, cell or gene together), so STAR's own sorting is exercised.
- Choose reads whose truth is unambiguous: `unique_locus` for single-locus reads, `min_block >= 25` for spliced reads without annotation, `deletion_is_unambiguous` for deletions, mismatches away from read ends.
- Python 3.8 syntax, standard library only; indent with tabs.
- Tick the rows you cover in `TARGETS.md`.
