"""Area A, multimappers and seeds: --outMultimapperOrder Random with
--runRNGseed, primary choice, --winAnchorMultimapNmax, --seedPerWindowNmax,
the seed search start options, and --runThreadN 1 vs 4 on area-A options."""

import pytest

import helpers_a as ha
import star_helpers as sh
from star_helpers import ReadSet, record_truth, run_star

ATTRS = ["--outSAMattributes", "NH", "HI", "AS", "nM"]


def _repeat_reads(sim, n=36):
	"""Reads from the 2-, 3- and 12-copy segments (all copies identical) in
	turn, with a unique read after every fourth."""
	reads = []
	for i in range(n):
		seg = ("dup2", "dup3", "rep12")[i % 3]
		off = sim.rng.randint(0, sim.g.segments[seg].length - 100)
		reads.append(sim.segment_read(seg, off, strand=sim.rng.choice("+-")))
		if i % 4 == 3:
			reads.append(sim.genomic_read(*sim.unique_locus(), name="uniq%d" % i))
	return ReadSet(reads)


def _loci_in_order(aln):
	out = {}
	for r in aln:
		out.setdefault(r.qname, []).append(r)
	return out


# Old_2.4 output order of a multimapper's loci, as indices into the
# segment's copies (SEGMENT_SPECS order), for reads with the segment's
# forward sequence ('+') or its reverse complement ('-'). Taken from STAR
# 2.7.11b's output: parametersDefault calls this order "quasi-random".
OLD_ORDER = {
	("rep12", "+"): [2, 1, 3, 0, 8, 10, 11, 9, 7, 6, 4, 5],
	("rep12", "-"): [2, 1, 0, 3, 10, 9, 8, 11, 6, 7, 4, 5],
	("dup3", "+"): [2, 1, 0],
	("dup3", "-"): [1, 2, 0],
	("dup2", "+"): [0, 1],
	("dup2", "-"): [1, 0],
}


def test_old_multimapper_order(genome_plain, sim):
	"""Pins the current output: with the default --outMultimapperOrder
	Old_2.4 the loci of a multimapper come out in a fixed order that depends
	only on the read's sequence (OLD_ORDER; the same for every offset in the
	segment), HI counts 1..NH in that order and the first record is the
	primary. A rewrite that changes how windows or their alignments are
	ordered changes this order while keeping the set of loci."""
	reads = []
	for seg in ("rep12", "dup3", "dup2"):
		for i in range(6):
			off = sim.rng.randint(0, sim.g.segments[seg].length - 100)
			reads.append(sim.segment_read(seg, off, strand="+-"[i % 2], name="%s_%d" % (seg, i)))
	reads = ReadSet(reads)
	got = _loci_in_order(run_star(genome_plain, reads, ATTRS + ["--outFilterMultimapNmax", "12"]).alignments())
	assert sorted(got) == sorted(reads.by_name)
	for r in reads:
		recs = got[r.name]
		loci = r.extra["loci"]
		order = [loci.index((x.rname, x.pos, x.strand)) for x in recs]
		assert order == OLD_ORDER[r.extra["segment"], r.strand], r.name
		assert [x.tag("HI") for x in recs] == list(range(1, len(loci) + 1))
		assert [x.is_secondary for x in recs] == [False] + [True] * (len(loci) - 1)


@pytest.mark.parametrize("seed", [777, 4321])
def test_random_multimapper_order(genome_plain, sim, seed):
	"""--outMultimapperOrder Random shuffles each read's alignments with a
	Mersenne Twister seeded by runRNGseed * (thread index + 1) (ReadAlign.cpp)
	and a Fisher-Yates pass using int(u * itr + 0.5)
	(ReadAlign_multMapSelect.cpp); helpers_a.MT19937 reproduces the engine
	and libstdc++'s uniform_real_distribution. The order before shuffling is
	the Old_2.4 output order. With one thread the Random output must equal
	that order shuffled with draws taken read by read in input order (unique
	reads take no draws); HI numbers the records in output order and the
	first record is the primary. Pinned to --runThreadN 1 because each thread
	has its own generator."""
	reads = _repeat_reads(sim)
	base = ATTRS + ["--outFilterMultimapNmax", "12"]
	old = _loci_in_order(run_star(genome_plain, reads, base, threads=1).alignments())
	args = base + ["--outMultimapperOrder", "Random", "--runRNGseed", str(seed)]
	new = _loci_in_order(run_star(genome_plain, reads, args, threads=1).alignments())
	assert sorted(old) == sorted(new) == sorted(reads.by_name)
	rng = ha.MT19937(seed)
	moved = 0
	for r in reads:
		before = [(x.rname, x.pos, x.strand) for x in old[r.name]]
		after = [(x.rname, x.pos, x.strand) for x in new[r.name]]
		if r.name.startswith("uniq"):
			assert after == before == [r.truth()[:2] + (r.strand,)]
			continue
		assert sorted(before) == sorted(r.extra["loci"])
		assert after == ha.random_multimapper_order(before, len(before), rng)
		assert [x.tag("HI") for x in new[r.name]] == list(range(1, len(after) + 1))
		assert [x.is_secondary for x in new[r.name]] == [False] + [True] * (len(after) - 1)
		moved += after != before
	assert moved >= 10


def test_primary_is_best_scoring_alignment(genome_plain, sim):
	"""Copy 2 of segment dupmm has one substitution per 100 bases, so with
	--outFilterMultimapScoreRange 2 a read from copy 1 has alignments of AS
	98 (copy 1) and 96 (copy 2). OneBestScore marks only the best one
	primary (parametersDefault), with Old_2.4 order and with Random order;
	with Random the best alignment also comes first (multMapSelect moves the
	best-scoring alignments to the front before shuffling)."""
	reads = ReadSet([sim.segment_read("dupmm", 100 * sim.rng.randint(0, 8) + sim.rng.randint(5, 45),
		strand=sim.rng.choice("+-")) for _ in range(10)])
	base = ATTRS + ["--outFilterMultimapScoreRange", "2"]
	for args in (base, base + ["--outMultimapperOrder", "Random"]):
		by_name = run_star(genome_plain, reads, args, threads=1).alignments().by_name()
		assert sorted(by_name) == sorted(reads.by_name)
		for name, recs in by_name.items():
			truth = reads.by_name[name]
			primary = [x for x in recs if not x.is_secondary]
			assert len(recs) == 2 and len(primary) == 1
			assert record_truth(primary[0]) == truth.truth() and primary[0].tag("AS") == 98
			assert sorted(x.tag("AS") for x in recs) == [96, 98]
			if "Random" in args:
				assert recs[0] is primary[0]


def test_primary_is_best_when_written_second(genome_plain, sim):
	"""250-base reads from dupmm copy 1 that cover three of the sites where
	copy 2 differs and carry copy 2's base at the first: copy 1 has one
	mismatch (AS 249 - 1 - 2 = 246, length term of a 250-base span) and copy
	2 two (AS 248 - 2 - 2 = 244), so --outFilterMultimapScoreRange 2 reports
	both. OneBestScore (default) and AllBestScore (parametersDefault) both
	make exactly the copy-1 record primary. In the default Old_2.4 order at
	least one read has its AS-244 record written first (checked), so the
	primary flag is not taken from the output position."""
	c1, c2 = sim.g.segments["dupmm"].copies
	reads = []
	for i in range(8):
		off = 100 * i + 20
		s1 = sim.g.seq(c1.contig, c1.start + off, c1.start + off + 249)
		s2 = sim.g.seq(c2.contig, c2.start + off, c2.start + off + 249)
		sites = [j for j in range(250) if s1[j] != s2[j]]
		assert len(sites) == 3
		fwd = s1[:sites[0]] + s2[sites[0]] + s1[sites[0] + 1:]
		reads.append(ha.oriented("mm%d" % i, fwd, c1.contig, c1.start + off, "250M", "+-"[i % 2],
			extra={"other": (c2.contig, c2.start + off)}))
	reads = ReadSet(reads)
	base = ATTRS + ["--outFilterMultimapScoreRange", "2"]
	first_secondary = 0
	for args in (base, base + ["--outSAMprimaryFlag", "AllBestScore"]):
		recs = run_star(genome_plain, reads, args).records()
		by_name = {}
		for x in recs:
			by_name.setdefault(x.qname, []).append(x)
		assert sorted(by_name) == sorted(reads.by_name)
		for name, rs in by_name.items():
			r = reads.by_name[name]
			got = sorted((x.rname, x.pos, x.strand, x.tag("AS"), x.tag("nM"), x.is_secondary) for x in rs)
			assert got == sorted([(r.contig, r.pos, r.strand, 246, 1, False),
				r.extra["other"] + (r.strand, 244, 2, True)])
			first_secondary += rs[0].is_secondary
	assert first_secondary > 0


def test_win_anchor_multimap_nmax(genome_plain, sim):
	"""Only seeds mapping to at most --winAnchorMultimapNmax loci can anchor
	a window (ReadAlign_stitchPieces.cpp). With --outFilterMultimapNmax 20,
	reads of the 12-copy segment map to all 12 copies when the limit is 12
	and have no window at all ('unmapped: other') when it is 11; 2- and
	3-copy reads and unique reads map in both runs."""
	reads = []
	for i in range(4):
		reads.append(sim.segment_read("rep12", sim.rng.randint(0, 200), strand="+-"[i % 2], name="rep%d" % i))
		reads.append(sim.segment_read("dup3", sim.rng.randint(0, 500), strand="+-"[i % 2], name="dup3_%d" % i))
		reads.append(sim.segment_read("dup2", sim.rng.randint(0, 1400), strand="+-"[i % 2], name="dup2_%d" % i))
		reads.append(sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2], name="uniq%d" % i))
	reads = ReadSet(reads)
	for limit in (11, 12):
		run = run_star(genome_plain, reads, ATTRS + ["--outFilterMultimapNmax", "20", "--winAnchorMultimapNmax", str(limit)])
		by_name = run.alignments().by_name()
		want = [n for n in reads.by_name if limit == 12 or not n.startswith("rep")]
		assert sorted(by_name) == sorted(want)
		for name in want:
			r = reads.by_name[name]
			loci = r.extra["loci"] if "loci" in r.extra else [r.truth()[:2] + (r.strand,)]
			assert sorted((x.rname, x.pos, x.strand) for x in by_name[name]) == sorted(loci)
			assert all(x.tag("NH") == len(loci) for x in by_name[name])
		assert ha.unmapped_counts(run.log_final()) == (0, 0, 0 if limit == 12 else 4, 0)


def test_seed_per_window_nmax(genome_plain, sim):
	"""--seedPerWindowNmax 1: when a second anchor seed reaches a window that
	already holds its one allowed seed, and every stored seed is an anchor,
	STAR drops all windows of the read (ReadAlign_assignAlignToWindow.cpp,
	MARKER_TOO_MANY_ANCHORS_PER_WINDOW). Exact reads, unique or from the
	3-copy segment (copies on different chromosomes), need one seed per
	window and still map; a read with a mismatch at offset 50 (seeds on both
	sides) or spliced 50/50 (one seed per exon) is unmapped as 'other'."""
	reads = [sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2], name="exact%d" % i) for i in range(4)]
	reads += [sim.segment_read("dup3", sim.rng.randint(0, 500), strand="+-"[i % 2], name="dup3_%d" % i)
		for i in range(3)]
	reads += [sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2], mismatches=(50,), name="mm%d" % i)
		for i in range(4)]
	reads += [ha.tx_junction_read(sim, tx, 0, 50, sense=bool(i), name="spl%s%d" % (tx, i))
		for tx in ("T05", "T09") for i in range(2)]
	reads = ReadSet(reads)
	run = run_star(genome_plain, reads, ATTRS + ["--seedPerWindowNmax", "1"])
	by_name = run.alignments().by_name()
	assert sorted(by_name) == sorted(n for n in reads.by_name if n.startswith(("exact", "dup3")))
	for name, recs in by_name.items():
		r = reads.by_name[name]
		loci = r.extra["loci"] if "loci" in r.extra else [r.truth()[:2] + (r.strand,)]
		assert sorted((x.rname, x.pos, x.strand) for x in recs) == sorted(loci)
	assert ha.unmapped_counts(run.log_final()) == (0, 0, 8, 0)


def _mixed_truth_reads(sim):
	reads = [sim.genomic_read(*sim.unique_locus(120), strand="+-"[i % 2], mismatches=mm)
		for i, mm in enumerate([(), (), (25, 71), (30, 60), (10, 50, 90), (99,), (0,)])]
	reads += list(sim.transcript_reads(10, spliced=True))
	c, p = sim.unique_locus(120)
	while not sim.deletion_is_unambiguous(c, p + 45, 3):
		c, p = sim.unique_locus(120)
	reads.append(sim.genomic_read(c, p, "+", deletion=(45, 3)))
	reads.append(sim.genomic_read(*sim.unique_locus(120), strand="-", insertion=(52, 2)))
	return ReadSet(reads)


@pytest.mark.parametrize("args", [
	["--seedSearchStartLmax", "10"],
	["--seedSearchStartLmax", "30"],
	["--seedSearchStartLmax", "100"],
	["--seedSearchStartLmaxOverLread", "0.3"],
	["--seedSearchLmax", "20"],
], ids=lambda a: "_".join(a).lstrip("-"))
def test_seed_search_options_keep_true_alignments(genome_plain, sim, args):
	"""Changing where seed searches start (--seedSearchStartLmax and its
	OverLread form: the read is cut into pieces of at most that length,
	searched from both ends) or adding fixed-length searches
	(--seedSearchLmax) only changes sensitivity; reads with a few
	substitutions, splices with blocks >= 25 bases and short indels still get
	their true alignment. This catches start-point arithmetic that loses or
	corrupts seeds (ReadAlign_mapOneRead.cpp); a change of start points
	that still finds the true alignments passes (for example 2 pieces of 50
	instead of 3 of 33 for a 100-base read, which also satisfies "pieces no
	longer than this value")."""
	reads = _mixed_truth_reads(sim)
	aln = run_star(genome_plain, reads, ATTRS + args).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		r = reads.by_name[rec.qname]
		if r.mismatches and (0 in [p - r.pos for p in r.mismatches] or 99 in [p - r.pos for p in r.mismatches]):
			left = 1 if r.mismatches[0] == r.pos else 0
			right = 1 if r.mismatches[-1] == r.pos + 99 else 0
			assert record_truth(rec) == (r.contig, r.pos + left, sh.ops_to_cigar(
				[(left, "S"), (100 - left - right, "M"), (right, "S")]), r.strand)
		else:
			assert record_truth(rec) == r.truth()


def test_thread_count_area_a_options(genome_plain, sim):
	"""--runThreadN 1 and 4 give the same record multiset and Log.final.out
	statistics with area-A options switched on: adapter clipping,
	Extend5pOfRead1, multimappers up to 12 loci, unmapped reads recorded
	(single-end), and mate merging (--peOverlapNbasesMin) for pairs.
	--limitIObufferSize 106000 per mate leaves input chunks of about 20
	reads, so the four threads map different reads; with the default buffer
	one thread would load the whole input as a single chunk."""
	se = ReadSet(list(_repeat_reads(sim, 12)) + list(_mixed_truth_reads(sim))
		+ [sim.random_read() for _ in range(4)])
	pe = sim.transcript_pairs(40, insert=(120, 260))
	skip = ("Started job on", "Started mapping on", "Finished on", "Mapping speed, Million of reads per hour")
	for reads, args in (
			(se, ATTRS + ["--outSAMunmapped", "Within", "--outFilterMultimapNmax", "12", "--alignEndsType",
				"Extend5pOfRead1", "--clip3pAdapterSeq", sh.ILLUMINA_ADAPTER]),
			(pe, ATTRS + ["--peOverlapNbasesMin", "10", "--outSAMunmapped", "Within"])):
		args = args + ["--limitIObufferSize", str(106000 * (2 if reads.paired else 1)), "50000000"]
		one = run_star(genome_plain, reads, args, threads=1)
		four = run_star(genome_plain, reads, args, threads=4)
		a, b = sh.record_multiset(one.records()), sh.record_multiset(four.records())
		assert a == b, sh.first_line_diff(a, b)
		la = {k: v for k, v in one.log_final().items() if k not in skip}
		lb = {k: v for k, v in four.log_final().items() if k not in skip}
		assert la == lb
		assert la["Number of input reads"] == len(reads)


def test_transcripts_per_window_nmax(genome_plain, sim):
	"""--alignTranscriptsPerWindowNmax 1 keeps one alignment per window. A
	window covers whole genome bins (2^winBinNbits = 65536 bases here, more
	than any contig), so the copies of a segment on one chromosome and
	strand share a window: reads of the 12-copy segment (4 copies on each
	chromosome) and of the 3-copy segment (one per chromosome) each get 3
	alignments (NH 3), one per chromosome, every one at a true copy. The
	12-copy reads then pass the default --outFilterMultimapNmax 10."""
	reads = [sim.segment_read(seg, sim.rng.randint(0, 200), strand="+-"[i % 2], name="%s_%d" % (seg, i))
		for seg in ("rep12", "dup3") for i in range(3)]
	reads = ReadSet(reads)
	by_name = run_star(genome_plain, reads, ATTRS + ["--alignTranscriptsPerWindowNmax", "1"]).alignments().by_name()
	assert sorted(by_name) == sorted(reads.by_name)
	for name, recs in by_name.items():
		loci = reads.by_name[name].extra["loci"]
		assert sorted(r.rname for r in recs) == ["chr1", "chr2", "chr3"]
		assert all((r.rname, r.pos, r.strand) in loci and r.tag("NH") == 3 for r in recs)


def test_seed_multimap_nmax(genome_plain, sim):
	"""Seeds mapping to more than --seedMultimapNmax loci are not stored
	(ReadAlign_storeAligns.cpp). With 11, reads of the 12-copy segment have
	no seed and are unmapped as 'other', while 2- and 3-copy reads keep all
	their loci (--outFilterMultimapNmax 20). (Whether a seed with exactly
	that many loci is kept differs from parametersDefault's "fewer than"
	and is not tested here.)"""
	reads = [sim.segment_read(seg, sim.rng.randint(0, 200), strand="+-"[i % 2], name="%s_%d" % (seg, i))
		for seg in ("rep12", "dup3", "dup2") for i in range(3)]
	reads = ReadSet(reads)
	run = run_star(genome_plain, reads, ATTRS + ["--outFilterMultimapNmax", "20", "--seedMultimapNmax", "11"])
	by_name = run.alignments().by_name()
	assert sorted(by_name) == sorted(n for n in reads.by_name if not n.startswith("rep12"))
	for name, recs in by_name.items():
		assert sorted((r.rname, r.pos, r.strand) for r in recs) == sorted(reads.by_name[name].extra["loci"])
	assert ha.unmapped_counts(run.log_final()) == (0, 0, 3, 0)
