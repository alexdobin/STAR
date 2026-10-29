"""Area C: SJ.out.tab columns for every motif class, the --outSJfilter*
options, --outSJfilterReads, --outSJtype None and --outFilterType BySJout.

SJ.out.tab columns (manual 5.5): contig, first and last intron base, strand
(0 undefined, 1 +, 2 -), motif (0 non-canonical, 1 GT/AG, 2 CT/AC, 3 GC/AG,
4 CT/GC, 5 AT/AC, 6 GT/AT), annotated, unique reads, multi-mapping reads,
maximum spliced alignment overhang. The filters take four values, one per
motif class: (1) non-canonical, (2) GT/AG and CT/AC, (3) GC/AG and CT/GC,
(4) AT/AC and GT/AT; none applies to annotated junctions (parametersDefault,
"Output Filtering: Splice Junctions").
"""

import pytest

import helpers_c as hc
import star_helpers as sh
from star_helpers import ReadSet, record_multiset, record_truth, run_star

# (contig, intron, [(left, right, strand) of one read each]); left/right are
# the exonic bases before/after the intron in genome orientation.
COLUMN_SPECS = (
	("chr1", (1301, 2300), [(50, 50, "+"), (85, 15, "-")]),  # T01a, GT/AG
	("chr1", (2451, 3450), [(40, 60, "+")]),  # T01a
	("chr1", (1301, 3450), [(30, 70, "+"), (70, 30, "-"), (64, 36, "+"), (45, 55, "-")]),  # T01b
	("chr1", (6401, 7000), [(25, 75, "-")]),  # T02, CT/AC
	("chr1", (7201, 8000), [(60, 40, "+"), (33, 67, "+")]),  # T02
	("chr1", (11401, 11800), [(30, 70, "+"), (85, 15, "-"), (50, 50, "+")]),  # T03, GC/AG
	("chr1", (15301, 15900), [(20, 80, "-"), (77, 23, "+")]),  # T04, AT/AC
	("chr1", (20601, 21000), [(48, 52, "+")]),  # T05
	("chr1", (21801, 22300), [(70, 30, "-"), (35, 65, "+"), (55, 45, "-"), (61, 39, "+")]),  # T06, CT/AC
	("chr2", (2301, 3000), [(44, 56, "+"), (75, 25, "-")]),  # T08, CT/GC
	("chr2", (3151, 4000), [(66, 34, "-")]),  # T08
	("chr2", (4201, 5000), [(18, 82, "+"), (82, 18, "-"), (50, 50, "-")]),  # T08
	("chr2", (12301, 17300), [(57, 43, "+"), (43, 57, "-")]),  # T09, 5 kb
	("chr2", (30401, 30900), [(26, 74, "+")]),  # T10
	("chr2", (31301, 31800), [(40, 60, "-"), (80, 20, "+")]),  # T10
	("chr2", (30401, 31800), [(52, 48, "+"), (31, 69, "-"), (47, 53, "+")]),  # T10skip, not in the GTF
	("chr3", (12501, 13000), [(37, 63, "-"), (63, 37, "+"), (28, 72, "+"), (51, 49, "-")]),  # T11, GT/AT
	("chr3", (15201, 15800), [(50, 50, "+"), (35, 65, "-"), (62, 38, "+")]),  # T12nc, non-canonical
)


def dup2_reads(sim, multi=2):
	"""`multi` reads crossing an intron inside the dup2 segment, which align
	equally well to both copies (chr1 +, chr2 -), and one read at the chr1
	copy whose left block starts 10 bases before the segment, so only chr1
	fits it. Returns (reads, chr1 intron)."""
	a, b = hc.segment_intron(sim.g, "dup2")
	reads = [hc.segment_junction_read(sim, "dup2", a, b, 50, 50, copy=i % 2, strand="+-"[i % 2])
		for i in range(multi)]
	contig, s, e = hc.segment_intron_loci(sim.g, "dup2", a, b)[0]
	assert s - 70 == sim.g.segments["dup2"].copies[0].start - 10
	reads.append(hc.junction_read(sim, contig, (s, e), 70, 30, "+"))
	return reads, (contig, s, e)


def column_reads(sim):
	reads = []
	for contig, intron, specs in COLUMN_SPECS:
		for left, right, strand in specs:
			reads.append(hc.junction_read(sim, contig, intron, left, right, strand))
	# 200 bases over both T01a introns: blocks of 20, 150 and 30 bases
	reads.append(sim.transcript_read("T01a", 280, length=200))
	reads += dup2_reads(sim)[0]
	return reads


@pytest.mark.parametrize("annotated", [False, True], ids=["plain", "sjdb"])
def test_sj_out_columns_every_motif(genome_plain, genome_sjdb, sim, annotated):
	"""Reads across the 16 GTF introns (all six canonical motif codes), the
	unannotated T10skip (GT/AG) and T12nc (non-canonical, 3 reads, overhangs
	>= 30) introns, and an intron inside the duplicated dup2 segment (two
	reads that map to both copies, one that fits only the chr1 copy). Every
	SJ.out.tab row is worked out from the truth: strand and motif from the
	intron's genomic dinucleotides (CT/AC, CT/GC and GT/AT are - strand
	introns, non-canonical is strand 0), unique and multi-mapping read
	counts, and column 9 = the largest min(left block, right block) over the
	reads; for the 200-base read crossing two introns the blocks are the
	neighbouring exons (20|150 and 150|30), not the read ends. Without
	annotation every junction is unannotated and passes the default filters
	(overhangs >= 15 > 12, the non-canonical one has 3 unique reads and
	overhang 50 >= 30); with the GTF inserted the GTF introns have column 6 =
	1 and the others 0."""
	for contig, (s, e), _ in COLUMN_SPECS:
		assert hc.junction_is_clean(sim.g, contig, s, e, 15, 15), (contig, s, e)
	reads = column_reads(sim)
	genome = genome_sjdb if annotated else genome_plain
	run = run_star(genome, ReadSet(reads))
	expected = hc.expected_sj(hc.truth_alignments(reads), sim.g,
		hc.annotated_junctions(sim.g) if annotated else ())
	assert {r.motif for r in expected} == set(range(7))
	assert run.sj() == expected
	recs = run.alignments().by_name()
	for r in reads:
		if r.kind == "multi":
			assert sorted(record_truth(x) for x in recs[r.name]) == sorted(r.extra["loci"])
		else:
			assert [record_truth(x) for x in recs[r.name]] == [r.truth()], r.name


# Input buffer for multi-chunk runs: STAR reads at most limitIObufferSize[0]
# / mates - 2 * 651 - 2 * 50000 bytes per chunk (Parameters.cpp,
# IncludeDefine.h), 8698 bytes here, so a few hundred reads make several
# chunks that the threads map separately. With the default 30 MB buffer the
# test read sets fit in one chunk and only one thread gets reads.
SMALL_CHUNKS = ["--limitIObufferSize", "110000", "50000000"]
TIME_FIELDS = ("Started job on", "Started mapping on", "Finished on", "Mapping speed, Million of reads per hour")


def test_outputs_do_not_depend_on_input_chunks(genome_sjdb, sim):
	"""Several input chunks mapped by 4 threads (SMALL_CHUNKS) give the
	outputs of one chunk mapped by 1 thread, so the merging of per-thread
	junctions (outputSJ), statistics, gene counts and transcriptome records
	is checked: SJ.out.tab equals the truth (as in
	test_sj_out_columns_every_motif[sjdb], plus the junctions of 250
	transcript reads) and the 1-thread output; ReadsPerGene.out.tab, every
	Log.final.out field except times and speed, the alignment records and
	the transcriptome records (primary flag masked, chosen at random:
	CHANGES.md 2.5.1a) are the same."""
	reads = column_reads(sim) + list(sim.transcript_reads(250)) + \
		[sim.segment_read("dup3", 30 * i) for i in range(10)] + [sim.random_read() for _ in range(10)]
	args = ["--quantMode", "GeneCounts", "TranscriptomeSAM"]
	one = run_star(genome_sjdb, ReadSet(reads), args, threads=1)
	four = run_star(genome_sjdb, ReadSet(reads), args + SMALL_CHUNKS, threads=4)
	mapped = [r for r in reads if r.kind != "random"]
	expected = hc.expected_sj(hc.truth_alignments(mapped), sim.g, hc.annotated_junctions(sim.g))
	assert one.sj() == expected
	assert four.path("SJ.out.tab").read_bytes() == one.path("SJ.out.tab").read_bytes()
	assert four.path("ReadsPerGene.out.tab").read_bytes() == one.path("ReadsPerGene.out.tab").read_bytes()
	a, b = one.log_final(), four.log_final()
	assert [k for k in a if k not in TIME_FIELDS] == [k for k in b if k not in TIME_FIELDS]
	assert {k: v for k, v in a.items() if k not in TIME_FIELDS} == {k: v for k, v in b.items() if k not in TIME_FIELDS}
	assert a["Number of input reads"] == len(reads)
	assert record_multiset(four.records()) == record_multiset(one.records())
	tx = lambda run: sorted((r.qname, r.rname, r.pos, r.cigar, r.flag & ~sh.FLAG_SECONDARY, r.tag("NH"))
		for r in run.transcriptome())
	assert tx(four) == tx(one)


def test_sj_counts_read_pairs_once(genome_sjdb, sim):
	"""SJ.out.tab counts reads (manual 5.5, column 7 "number of uniquely
	mapping reads crossing the junction"), so a pair whose overlapping
	mates both cross the T01a junction chr1:1301-2300 adds 1, not 2. Pair 1:
	mate 1 at T01a offsets 225-324 (blocks 75|25), mate 2 at 250-349
	(50|50). Pair 2 (mate 1 antisense): the mate at 260-359 crosses
	1301-2300 (40|60), the mate at 400-499 crosses 2451-3450 (50|50). Rows by hand: 1301-2300 unique 2,
	column 9 = max(25, 50, 40) = 50; 2451-3450 unique 1, column 9 50."""
	pairs = [sim.transcript_pair("T01a", 225, 125, forward=True), sim.transcript_pair("T01a", 260, 240, forward=False)]
	run = run_star(genome_sjdb, ReadSet(pairs))
	recs = run.alignments().by_name()
	for p in pairs:
		assert sorted(record_truth(x) for x in recs[p.name]) == sorted([p.r1.truth(), p.r2.truth()])
	assert sorted(x.introns for p in pairs for x in (p.r1, p.r2)) == [[(1301, 2300)]] * 3 + [[(2451, 3450)]]
	assert run.sj() == [sh.SJRecord("chr1", 1301, 2300, 1, 1, 1, 2, 0, 50), sh.SJRecord("chr1", 2451, 3450, 1, 1, 1, 1, 0, 50)]


def test_sj_collapse_with_small_buffers(genome_plain, sim):
	"""--limitOutSJcollapsed 8 --limitOutSJoneRead 2: each thread's junction
	buffer holds 8 records, so it is collapsed every few reads and enlarged
	when 2 x limitOutSJoneRead records no longer fit after collapsing
	(ReadAlignChunk_mapChunk.cpp), and the final merge buffer of
	2 x limitOutSJcollapsed rows is enlarged for the 20 output rows
	(outputSJ.cpp). No read records more than 2 junctions (the 200-base read
	crosses 2 introns, a dup2 read 1 intron at 2 loci), and SJ.out.tab is
	the one worked out from the truth, as in
	test_sj_out_columns_every_motif[plain]."""
	reads = column_reads(sim)
	run = run_star(genome_plain, ReadSet(reads), "--limitOutSJcollapsed 8 --limitOutSJoneRead 2")
	expected = hc.expected_sj(hc.truth_alignments(reads), sim.g)
	assert len(expected) == 20
	assert run.sj() == expected
	assert "Increased the size of chunkOutSJ" in run.log_out()


# ---------------------------------------------------------------- filters

# Labelled junctions for the filter tests: label -> (contig, intron, reads as
# (left, right, strand)). A..H cover every motif class; J/K share the donor
# 1301, L/M share the acceptor 31800. The dup2 junction (I1 on chr1, I2 on
# chr2) is added separately.
FILTER_SPECS = (
	("A", "chr1", (20601, 21000), [(88, 12, "+")]),  # GT/AG, u1, overhang 12, gap 400
	("B", "chr2", (12301, 17300), [(40, 60, "+"), (70, 30, "-"), (45, 55, "+")]),  # GT/AG, u3, oh 45, gap 5000
	("C", "chr1", (6401, 7000), [(25, 75, "-"), (50, 50, "+")]),  # CT/AC, u2, oh 50, gap 600
	("D", "chr1", (11401, 11800), [(20, 80, "+")]),  # GC/AG, u1, oh 20, gap 400
	("E", "chr2", (2301, 3000), [(35, 65, "-"), (60, 40, "+")]),  # CT/GC, u2, oh 40, gap 700
	("F", "chr1", (15301, 15900), [(30, 70, "+")]),  # AT/AC, u1, oh 30, gap 600
	("G", "chr3", (12501, 13000), [(45, 55, "-"), (30, 70, "+")]),  # GT/AT, u2, oh 45, gap 500
	("H", "chr3", (15201, 15800), [(50, 50, "+"), (40, 60, "-"), (65, 35, "+")]),  # non-canonical, u3, oh 50, gap 600
	("J", "chr1", (1301, 2300), [(50, 50, "+")]),  # GT/AG, u1, oh 50, gap 1000
	("K", "chr1", (1301, 3450), [(60, 40, "-"), (30, 70, "+")]),  # GT/AG, u2, oh 40, gap 2150
	("L", "chr2", (31301, 31800), [(50, 50, "+")]),  # GT/AG, u1, oh 50, gap 500
	("M", "chr2", (30401, 31800), [(60, 40, "-")]),  # GT/AG, u1, oh 40, gap 1400; not in the GTF
)
NOVEL_IN_GTF_GENOME = {"H", "I1", "I2", "M"}

# (id, extra args, genome, labels whose rows are dropped), each worked by hand
# from the per-junction counts, overhangs and gaps above. I1 is chr1 (u1 m2
# oh 50 gap 206), I2 chr2 (u0 m2 oh 50 gap 206).
FILTER_CASES = [
	("defaults", "", "plain", set()),
	# overhang >= 30/13/21/31 by class: A 12 < 13, D 20 < 21, F 30 < 31
	("overhang", "--outSJfilterOverhangMin 30 13 21 31", "plain", {"A", "D", "F"}),
	# -1 = no output for that motif class: GC/AG (D) and CT/GC (E)
	("overhang_off", "--outSJfilterOverhangMin 30 12 -1 12", "plain", {"D", "E"}),
	# kept if unique >= U or total >= T: I1 survives on its total of 3 only
	("counts", "--outSJfilterCountUniqueMin 4 2 2 2 --outSJfilterCountTotalMin 4 3 2 3", "plain",
		{"A", "D", "F", "H", "I2", "J", "L", "M"}),
	("counts_off", "--outSJfilterCountUniqueMin -1 -1 -1 -1 --outSJfilterCountTotalMin -1 -1 -1 -1", "plain",
		{"A", "B", "C", "D", "E", "F", "G", "H", "I1", "I2", "J", "K", "L", "M"}),
	# GT/AG junctions need donors and acceptors >= 1 apart: J/K share a donor,
	# L/M an acceptor
	("dist", "--outSJfilterDistToOtherSJmin 10 1 5 10", "plain", {"J", "K", "L", "M"}),
	# 1 read: gap <= 600; 2 reads: gap <= 4999; 3+ reads: any gap. F (600, 1
	# read) sits on the limit.
	("intronmax", "--outSJfilterIntronMaxVsReadN 600 4999", "plain", {"J", "M"}),
	("intronmax_edge", "--outSJfilterIntronMaxVsReadN 599 2149", "plain", {"F", "J", "K", "M"}),
	# every filter set to remove everything: only GTF junctions are left
	("annotated_exempt", "--outSJfilterOverhangMin -1 -1 -1 -1 --outSJfilterCountUniqueMin -1 -1 -1 -1 "
		"--outSJfilterCountTotalMin -1 -1 -1 -1 --outSJfilterDistToOtherSJmin 0 1 0 0 "
		"--outSJfilterIntronMaxVsReadN 1", "sjdb", NOVEL_IN_GTF_GENOME),
]


@pytest.fixture(scope="module")
def filter_reads(synth_genome):
	"""(reads, {(contig, start, end): label}) for the filter tests."""
	sim = sh.ReadSimulator(synth_genome, seed=4101)
	reads = []
	labels = {}
	for label, contig, intron, specs in FILTER_SPECS:
		labels[(contig,) + intron] = label
		for left, right, strand in specs:
			reads.append(hc.junction_read(sim, contig, intron, left, right, strand))
	extra, chr1_intron = dup2_reads(sim)
	reads += extra
	labels[chr1_intron] = "I1"
	a, b = hc.segment_intron(synth_genome, "dup2")
	labels[hc.segment_intron_loci(synth_genome, "dup2", a, b)[1]] = "I2"
	return ReadSet(reads), labels


@pytest.fixture(scope="module")
def filter_default_run(genome_plain, filter_reads):
	return run_star(genome_plain, filter_reads[0])


def _label(row, labels):
	return labels[(row.contig, row.start, row.end)]


@pytest.mark.parametrize("case", FILTER_CASES, ids=[c[0] for c in FILTER_CASES])
def test_sj_filters(genome_plain, genome_sjdb, filter_reads, filter_default_run, case):
	"""Each --outSJfilter* setting drops exactly the junctions worked out by
	hand in FILTER_CASES from the documented rule (parametersDefault):
	OverhangMin compares column 9 with the class minimum; CountUniqueMin OR
	CountTotalMin must be met; DistToOtherSJmin is the distance to the
	nearest other donor and acceptor; IntronMaxVsReadN[n-1] caps the gap of
	junctions with n reads; -1 means no output for that class; annotated
	junctions are exempt from all of them. Rows that survive are unchanged."""
	name, args, which, dropped = case
	reads, labels = filter_reads
	assert set(labels.values()) >= dropped
	annotated = hc.annotated_junctions(genome_plain.synth) if which == "sjdb" else ()
	full = hc.expected_sj(hc.truth_alignments(reads), genome_plain.synth, annotated)
	assert {_label(r, labels) for r in full} == set(labels.values())
	if which == "sjdb":
		assert {_label(r, labels) for r in full if not r.annotated} == NOVEL_IN_GTF_GENOME
	if name == "defaults":
		run = filter_default_run
	else:
		run = run_star(genome_sjdb if which == "sjdb" else genome_plain, reads, args)
	assert run.sj() == [r for r in full if _label(r, labels) not in dropped]


def test_sj_filter_reads_unique(genome_plain, filter_reads):
	"""--outSJfilterReads Unique collects junctions from uniquely mapping
	reads only: the chr2 dup2 row (multi-mappers only) disappears and the
	chr1 row keeps its one unique read, with multi 0 and the overhang of that
	read (min(70, 30) = 30) instead of 50 from the multi-mappers."""
	reads, labels = filter_reads
	run = run_star(genome_plain, reads, "--outSJfilterReads Unique")
	unique_only = [a for a in hc.truth_alignments(reads) if a[3]]
	expected = hc.expected_sj(unique_only, genome_plain.synth)
	i1 = [r for r in expected if _label(r, labels) == "I1"]
	assert [r[6:] for r in i1] == [(1, 0, 30)]
	assert run.sj() == expected


def test_sj_type_none(genome_plain, filter_reads, filter_default_run):
	"""--outSJtype None writes no SJ.out.tab and leaves the alignments and
	the Log.final.out splice counts as with the default (Standard)."""
	reads = filter_reads[0]
	run = run_star(genome_plain, reads, "--outSJtype None")
	assert not run.exists("SJ.out.tab")
	assert filter_default_run.exists("SJ.out.tab")
	assert record_multiset(run.records()) == record_multiset(filter_default_run.records())
	a, b = run.log_final(), filter_default_run.log_final()
	splices = [k for k in a if k.startswith("Number of splices")]
	assert len(splices) == 6 and all(a[k] == b[k] for k in splices)


# ---------------------------------------------------------------- BySJout

@pytest.fixture(scope="module")
def bysjout_reads(synth_genome):
	"""Reads for the BySJout tests (unannotated genome): 12 unspliced exonic
	reads; 2 reads over the T05 intron with an 11-base overhang (junction
	fails --outSJfilterOverhangMin 12); 1 read over the T02 intron with a
	12-base overhang (passes); 2 reads over the non-canonical T12nc intron
	(fails --outSJfilterCountUniqueMin 3); 2 dup2 multi-mapping spliced reads
	(passes, total count 2 >= 1)."""
	sim = sh.ReadSimulator(synth_genome, seed=4102)
	reads = list(sim.transcript_reads(12, spliced=False))
	reads.append(hc.junction_read(sim, "chr1", (20601, 21000), 89, 11, "+", name="oh11a"))
	reads.append(hc.junction_read(sim, "chr1", (20601, 21000), 11, 89, "-", name="oh11b"))
	reads.append(hc.junction_read(sim, "chr1", (6401, 7000), 88, 12, "+", name="oh12"))
	reads.append(hc.junction_read(sim, "chr3", (15201, 15800), 50, 50, "+", name="nc1"))
	reads.append(hc.junction_read(sim, "chr3", (15201, 15800), 40, 60, "-", name="nc2"))
	reads += dup2_reads(sim)[0][:2]
	return ReadSet(reads)


@pytest.fixture(scope="module")
def bysjout_normal(genome_plain, bysjout_reads):
	return run_star(genome_plain, bysjout_reads)


def _oh11_unspliced(synth):
	"""Expected (name -> (pos, cigar)) of the two 11-base-overhang reads once
	their junction is banned: the 89-base block plus however far a local
	extension (+1 match, -1 mismatch) reaches into the intron."""
	g = synth.contigs["chr1"]
	read_a = g[20511:20600] + g[21000:21011]  # 20512-20600 + 21001-21011
	ext_a = hc.extension_length(read_a[89:], g[20600:20611])
	read_b = g[20589:20600] + g[21000:21089]  # 20590-20600 + 21001-21089
	ext_b = hc.extension_length(read_b[:11][::-1], g[20989:21000][::-1])
	return {"oh11a": (20512, "%dM%dS" % (89 + ext_a, 11 - ext_a)),
		"oh11b": (21001 - ext_b, "%dS%dM" % (11 - ext_b, 89 + ext_b))}


@pytest.mark.parametrize("fmt", list(sh.OUTPUT_FORMATS))
def test_filter_type_bysjout(genome_plain, bysjout_reads, bysjout_normal, fmt):
	"""--outFilterType BySJout keeps only reads whose junctions passed the
	SJ.out.tab filters (parametersDefault). With Normal filtering every read
	has its true alignment. With BySJout: the two T12nc reads lose their
	only junction and their best remaining alignment (about 50 bases) is
	below 0.66 of the read, so they are unmapped "too short"; the two
	11-base-overhang reads are realigned without the junction, unspliced
	with the overhang soft-clipped (_oh11_unspliced); all other reads keep
	their alignments and SJ.out.tab is unchanged. Log.final.out by hand: 19
	input reads, 15 unique (was 17), 2 too short, 2 multi; splices of unique
	reads 1 GT/AG (oh12) instead of 3 GT/AG + 2 non-canonical."""
	reads = bysjout_reads
	normal = bysjout_normal.alignments().by_name()
	for r in reads:
		if r.kind == "multi":
			assert sorted(record_truth(x) for x in normal[r.name]) == sorted(r.extra["loci"])
		else:
			assert [record_truth(x) for x in normal[r.name]] == [r.truth()], r.name
	run = run_star(genome_plain, reads, ["--outFilterType", "BySJout"] + sh.OUTPUT_FORMATS[fmt])
	recs = run.alignments().by_name()
	assert "nc1" not in recs and "nc2" not in recs
	clipped = _oh11_unspliced(genome_plain.synth)
	for name, (pos, cigar) in clipped.items():
		assert [(x.rname, x.pos, x.cigar) for x in recs[name]] == [("chr1", pos, cigar)]
	for r in reads:
		if r.name in clipped or r.name in ("nc1", "nc2"):
			continue
		assert record_multiset(recs[r.name], ("rname", "pos", "cigar", "flag")) == \
			record_multiset(normal[r.name], ("rname", "pos", "cigar", "flag")), r.name
	assert run.sj() == bysjout_normal.sj()
	a, b = hc.segment_intron(genome_plain.synth, "dup2")
	dup2 = sorted(hc.segment_intron_loci(genome_plain.synth, "dup2", a, b))
	assert [(x.contig, x.start, x.end) for x in run.sj()] == [("chr1", 6401, 7000)] + dup2
	log = run.log_final()
	mapped = 13 * 100 + sum(n for _, c in clipped.values() for n, op in sh.cigar_ops(c) if op == "M")
	assert log["Number of input reads"] == 19
	assert log["Uniquely mapped reads number"] == 15
	assert log["Number of reads unmapped: too short"] == 2
	assert log["Number of reads mapped to multiple loci"] == 2
	assert log["Number of splices: Total"] == 1
	assert log["Number of splices: GT/AG"] == 1
	assert log["Number of splices: Non-canonical"] == 0
	assert log["Average mapped length"] == hc.fmt2(mapped / 15.0)
	assert bysjout_normal.log_final()["Number of splices: Total"] == 5
