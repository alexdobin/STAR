"""Area D: chimeric detection, old scheme (--chimMultimapNmax 0) and the
multimapping scheme (--chimMultimapNmax > 0): thresholds, output types,
paired-end and merged-mate chimeras, Chimeric.out.junction columns.

Scores are worked from parametersDefault: a gapless segment of n matching
bases scores n + round(-0.25 * log2(n)) (scoreGenomicLengthLog2scale), e.g.
49 for 50 bases, 98 for 100; a non-GT/AG chimeric junction adds
chimScoreJunctionNonGTAG (-1) in the multimapping scheme.
"""

import subprocess

import pytest

import helpers_d as hd
import star_helpers as sh
from star_helpers import ReadSet, run_star

STRANDS = (("+", "+"), ("+", "-"), ("-", "+"), ("-", "-"))
NEW = ["--chimMultimapNmax", "10"]
EXTRA_FIELDS = ("num_chim_aln", "max_poss_aln_score", "non_chim_aln_score", "this_chim_aln_score",
	"bestall_chim_aln_score", "PEmerged_bool")


@pytest.fixture
def csim(synth_genome, request):
	"""Simulator seeded from the test's name, like conftest.sim, so a test
	draws the same reads whichever tests ran before it."""
	return hd.module_sim(synth_genome, __name__, request.node.name)


def _fixture_sim(synth, name):
	"""Simulator of one module-scoped fixture, seeded from its name."""
	return hd.module_sim(synth, __name__, name)


@pytest.fixture(scope="module")
def clean_reads(synth_genome):
	"""Four 50+50 chimeras (one per strand pair) between unique loci on
	different contigs, breakpoints that cannot move."""
	sim = _fixture_sim(synth_genome, "clean_reads")
	return ReadSet([hd.strict_chimeric_read(sim, s, name="clean%d" % i) for i, s in enumerate(STRANDS)])


@pytest.fixture(scope="module")
def gtag_read(synth_genome):
	"""50+50 chimera whose breakpoint is GT/AG (junction type 1)."""
	return hd.strict_chimeric_read(_fixture_sim(synth_genome, "gtag_read"), ("+", "-"), motif=1, name="gtag1")


def _junction_rows(run):
	if not run.exists("Chimeric.out.junction"):
		return []
	return run.chimeric()


def _cols(row):
	return [row[k] for k in sh.CHIMERIC_FIELDS]


def _extra(row):
	return [int(x) for x in row["extra"][:6]]


@pytest.mark.parametrize("scheme", ["old", "new"])
@pytest.mark.parametrize("option,value,reported", [
	("--chimSegmentMin", 48, True), ("--chimSegmentMin", 51, False),
	("--chimJunctionOverhangMin", 50, True), ("--chimJunctionOverhangMin", 51, False),
])
def test_segment_length_and_overhang(genome_plain, clean_reads, scheme, option, value, reported):
	"""50+50 chimeras: both segments (and so both junction overhangs) are
	50 bases. parametersDefault: chimSegmentMin "minimum length of chimeric
	segment", chimJunctionOverhangMin "minimum overhang for a chimeric
	junction". 48 and 50 pass, 51 rejects. (STAR 2.7.11b rejects 50-base
	segments at chimSegmentMin 49 and 50, contrary to the documented
	minimum; reported separately.)"""
	args = {"--chimSegmentMin": [], "--chimJunctionOverhangMin": ["--chimSegmentMin", "20"]}[option]
	args = args + [option, str(value)] + (NEW if scheme == "new" else [])
	rows = _junction_rows(run_star(genome_plain, clean_reads, args))
	if not reported:
		assert rows == []
		return
	assert sorted(_cols(r) for r in rows) == sorted(list(p.extra["junction"].values()) for p in clean_reads)


def test_old_scheme_junction_file_has_14_columns(genome_plain, clean_reads):
	"""Old scheme: no header, exactly the 14 columns of manual 6.4, and the
	same lines (in any order) at --runThreadN 1 and 4."""
	one = run_star(genome_plain, clean_reads, "--chimSegmentMin 20", threads=1)
	four = run_star(genome_plain, clean_reads, "--chimSegmentMin 20", threads=4)
	lines = one.path("Chimeric.out.junction").read_text().splitlines()
	assert all(len(l.split("\t")) == 14 for l in lines)
	assert sorted(lines) == sorted(four.path("Chimeric.out.junction").read_text().splitlines())


def test_new_scheme_columns(genome_plain, clean_reads, gtag_read):
	"""Multimapping scheme: a header line naming the 21 columns, then per
	chimera columns 1-14 plus num_chim_aln 1, max_poss_aln_score 100 (read
	length), non_chim_aln_score 49 (the best linear alignment, one 50-base
	segment), this/bestall_chim_aln_score 49+49-1 = 97 for a non-GT/AG
	junction and 98 for GT/AG (type 1), PEmerged_bool 0."""
	reads = ReadSet(clean_reads.records + [gtag_read])
	run = run_star(genome_plain, reads, ["--chimSegmentMin", "20"] + NEW)
	header = run.path("Chimeric.out.junction").read_text().splitlines()[0].split("\t")
	assert header == list(sh.CHIMERIC_FIELDS) + list(EXTRA_FIELDS) + ["readgrp"]
	rows = run.chimeric()
	assert sorted(r["read_name"] for r in rows) == sorted(reads.by_name)
	for row in rows:
		p = reads.by_name[row["read_name"]]
		assert _cols(row) == list(p.extra["junction"].values())
		score = 98 if p is gtag_read else 97
		assert _extra(row) == [1, 100, 49, score, score, 0]
	assert gtag_read.extra["junction"]["junction_type"] == 1


@pytest.mark.parametrize("scheme,read,option,value,reported", [
	# old scheme, non-GT/AG: summed segment score 49 + 49 = 98
	("old", "clean", "--chimScoreMin", 98, True), ("old", "clean", "--chimScoreMin", 99, False),
	("old", "clean", "--chimScoreDropMax", 2, True), ("old", "clean", "--chimScoreDropMax", 1, False),
	# multimapping scheme, GT/AG: 98, no junction penalty
	("new", "gtag", "--chimScoreMin", 98, True), ("new", "gtag", "--chimScoreMin", 99, False),
	("new", "gtag", "--chimScoreDropMax", 2, True), ("new", "gtag", "--chimScoreDropMax", 1, False),
	# multimapping scheme, non-GT/AG: 97
	("new", "clean", "--chimScoreMin", 97, True), ("new", "clean", "--chimScoreDropMax", 3, True),
	("new", "clean", "--chimScoreDropMax", 2, False),
])
def test_score_min_and_drop_max(genome_plain, clean_reads, gtag_read, scheme, read, option, value, reported):
	"""chimScoreMin "minimum total (summed) score of the chimeric segments"
	and chimScoreDropMax "max drop of chimeric score from the read length"
	at their boundaries, with the scores worked in the module docstring."""
	reads = ReadSet([gtag_read]) if read == "gtag" else ReadSet(clean_reads.records[:1])
	args = ["--chimSegmentMin", "20", option, str(value)] + (NEW if scheme == "new" else [])
	rows = _junction_rows(run_star(genome_plain, reads, args))
	assert len(rows) == (1 if reported else 0)
	if reported:
		assert _cols(rows[0]) == list(reads.records[0].extra["junction"].values())


def test_nonchim_score_drop_min(genome_plain, clean_reads):
	"""Multimapping scheme: chimeric detection runs only when the best
	linear score (49) is more than chimNonchimScoreDropMin below the read
	length (100): a drop of 51 is enough for 50 and too little for 52.
	(STAR 2.7.11b also runs detection at 51; reported separately.)"""
	reads = ReadSet(clean_reads.records[:2])
	for value, n in ((50, 2), (52, 0)):
		run = run_star(genome_plain, reads, ["--chimSegmentMin", "20", "--chimNonchimScoreDropMin", str(value)] + NEW)
		assert len(run.chimeric()) == n, value


def test_junction_non_gtag_penalty(genome_plain, clean_reads, gtag_read):
	"""chimScoreJunctionNonGTAG is added to the chimeric score of
	non-GT/AG junctions only: with -4 the multimapping scheme scores the
	non-canonical chimeras 49+49-4 = 94 and the GT/AG one 98. In the old
	scheme a -10 penalty takes the score below 100 - chimScoreDropMax 5,
	which -1 does not."""
	reads = ReadSet(clean_reads.records[:2] + [gtag_read])
	run = run_star(genome_plain, reads, ["--chimSegmentMin", "20", "--chimScoreJunctionNonGTAG", "-4"] + NEW)
	scores = {r["read_name"]: _extra(r)[3] for r in run.chimeric()}
	assert scores == {p.name: (98 if p is gtag_read else 94) for p in reads}
	old = ReadSet(clean_reads.records[:2])
	for pen, n in (("-1", 2), ("-10", 0)):
		run = run_star(genome_plain, old, ["--chimSegmentMin", "20", "--chimScoreDropMax", "5",
			"--chimScoreJunctionNonGTAG", pen])
		assert len(_junction_rows(run)) == n, pen


@pytest.fixture(scope="module")
def repeat_reads(synth_genome):
	"""50+50 chimeras whose breakpoint can slide by 1-4 bases (the last
	bases of A also continue B's genome, or the first bases of B continue
	A's), no GT/AG or CT/AC within 8 bases."""
	csim = _fixture_sim(synth_genome, "repeat_reads")
	reads = []
	while len(reads) < 6:
		ca, pa = csim.unique_locus(80, avoid_genes=True)
		cb, pb = csim.unique_locus(80, avoid_genes=True)
		if ca == cb:
			continue
		strands = STRANDS[len(reads) % 4]
		r = csim.chimeric_read((ca, pa + 10, strands[0], 50), (cb, pb + 10, strands[1], 50))
		if not 0 < sum(r.extra["repeat"]) <= 4 or any(r.extra["shift_motifs"].values()):
			continue
		reads.append(r)
	return ReadSet(reads)


@pytest.mark.parametrize("scheme", ["old", "new"])
def test_slidable_breakpoint_repeat_columns(genome_plain, repeat_reads, scheme):
	"""When the breakpoint can slide over a repeat of length L = left +
	right (manual 6.4 columns 8-9: repeat length to the left / right of the
	junction), every column is consistent with one breakpoint inside that
	range: helpers_d.shifted_junction for a shift s with -left <= s <=
	right, so the two repeat columns sum to L. STAR places the breakpoint
	at the donor end of the repeat (s = -left, left column 0): pins the
	current output."""
	run = run_star(genome_plain, repeat_reads, ["--chimSegmentMin", "20"] + (NEW if scheme == "new" else []))
	rows = run.chimeric()
	assert sorted(r["read_name"] for r in rows) == sorted(repeat_reads.by_name)
	for row in rows:
		r = repeat_reads.by_name[row["read_name"]]
		left, right = r.extra["repeat"]
		options = {s: list(hd.shifted_junction(r, s).values()) for s in range(-left, right + 1)}
		assert _cols(row) in options.values()
		assert _cols(row) == options[-left]
		assert row["repeat_left_lenA"] + row["repeat_right_lenB"] == left + right


@pytest.fixture(scope="module")
def three_piece_read(synth_genome):
	"""B1(25) + A(50) + B2(23), three unique loci on three contigs, clean
	breakpoints: chimera B1|A scores 24+49 = 73, A|B2 49+22 = 71."""
	csim = _fixture_sim(synth_genome, "three_piece_read")
	while True:
		locs = []
		for n in (25, 50, 23):
			while True:
				c, p = csim.unique_locus(n + 20, avoid_genes=True)
				if all(c != x[0] for x in locs):
					break
			locs.append((c, p + 10, csim.rng.choice("+-"), n))
		r = hd.multi_piece_read(csim, locs)
		if r:
			return r


def test_score_separation(genome_plain, three_piece_read):
	"""Old scheme: the best chimera (B1|A, 73) is reported only when the
	next one (A|B2, 71) is more than chimScoreSeparation below it:
	separation 1 reports B1|A with its columns, 3 reports nothing. (STAR
	2.7.11b also rejects at 2, a difference equal to the "minimum
	difference"; reported separately.)"""
	r = three_piece_read
	reads = ReadSet([r])
	base = ["--chimSegmentMin", "15", "--chimScoreDropMax", "40"]
	rows = _junction_rows(run_star(genome_plain, reads, base + ["--chimScoreSeparation", "1"]))
	assert [_cols(x) for x in rows] == [list(hd.junction_columns(r, 0).values())]
	assert _junction_rows(run_star(genome_plain, reads, base + ["--chimScoreSeparation", "3"])) == []


@pytest.mark.parametrize("score_range,n", [(1, 1), (2, 2)])
def test_multimap_score_range(genome_plain, three_piece_read, score_range, n):
	"""Multimapping scheme: the two chimeras of the three-piece read score
	73-1 = 72 and 71-1 = 70; chimMultimapScoreRange 1 keeps only the best,
	2 keeps both (num_chim_aln 2, bestall 72). max_poss_aln_score is the
	read length 98, non_chim_aln_score the 50-base segment's 49."""
	r = three_piece_read
	run = run_star(genome_plain, ReadSet([r]), ["--chimSegmentMin", "15", "--chimScoreDropMax", "40",
		"--chimMultimapScoreRange", str(score_range)] + NEW)
	rows = sorted(run.chimeric(), key=lambda x: -int(x["extra"][3]))
	want = [(list(hd.junction_columns(r, 0).values()), 72), (list(hd.junction_columns(r, 1).values()), 70)][:n]
	assert [(_cols(x), _extra(x)) for x in rows] == [(c, [n, 98, 49, s, 72, 0]) for c, s in want]


def test_main_segment_multimapping(genome_plain, csim):
	"""Main segment A (60 bases) in dup3, which has 3 copies, B (40) unique.
	Old scheme: chimMainSegmentMultNmax ("maximum number of multi-alignments
	for the main chimeric segment") 3 reports the chimera with A at one of
	its copies, 2 and 1 report nothing. Multimapping scheme: one line per A
	copy, num_chim_aln 3."""
	r = hd.strict_chimeric_read(csim, ("+", "-"), (60, 40), locus_a=lambda: hd.segment_locus(csim, "dup3", 60))
	(ca, pa, sa, na), b = r.extra["pieces"]
	copies = [(ca, pa, sa)] + hd.other_copy(csim, "dup3", ca, pa, sa, na)
	expected = {}
	for c, p, s in copies:
		alt = csim.chimeric_read((c, p, s, na), b, name=r.name)
		assert alt.seq == r.seq
		expected[c] = list(alt.extra["junction"].values())
	reads = ReadSet([r])
	rows = _junction_rows(run_star(genome_plain, reads, "--chimSegmentMin 20 --chimMainSegmentMultNmax 3"))
	assert len(rows) == 1 and _cols(rows[0]) == expected[rows[0]["chr_donorA"]]
	for m in (2, 1):
		assert _junction_rows(run_star(genome_plain, reads, "--chimSegmentMin 20 --chimMainSegmentMultNmax %d" % m)) == []
	rows = run_star(genome_plain, reads, ["--chimSegmentMin", "20"] + NEW).chimeric()
	assert sorted(_cols(x) for x in rows) == sorted(expected.values())
	assert {_extra(x)[0] for x in rows} == {3}


def test_multimapping_chimera_nmax(genome_plain, csim):
	"""Acceptor B (40 bases) in the exact repeat dup2 (2 copies), A (60)
	unique. Old scheme: two equally good chimeras, none reported
	(chimScoreSeparation). Multimapping scheme: chimMultimapNmax 1 is
	exceeded, 2 reports both, with B at each copy and num_chim_aln 2."""
	r = hd.strict_chimeric_read(csim, ("+", "+"), (60, 40), locus_b=lambda: hd.segment_locus(csim, "dup2", 40))
	a, (cb, pb, sb, nb) = r.extra["pieces"]
	expected = [list(r.extra["junction"].values())]
	for c, p, s in hd.other_copy(csim, "dup2", cb, pb, sb, nb):
		alt = csim.chimeric_read(a, (c, p, s, nb), name=r.name)
		assert alt.seq == r.seq
		expected.append(list(alt.extra["junction"].values()))
	reads = ReadSet([r])
	assert _junction_rows(run_star(genome_plain, reads, "--chimSegmentMin 20")) == []
	assert run_star(genome_plain, reads, "--chimSegmentMin 20 --chimMultimapNmax 1").chimeric() == []
	rows = run_star(genome_plain, reads, "--chimSegmentMin 20 --chimMultimapNmax 2").chimeric()
	assert sorted(_cols(x) for x in rows) == sorted(expected)
	assert [_extra(x) for x in rows] == [[2, 100, 59, 97, 97, 0]] * 2


@pytest.mark.parametrize("lengths", [(45, 50), (50, 45)])
def test_segment_read_gap_max(genome_plain, csim, lengths):
	"""A + 5 bases matching neither side + B, with the shorter segment
	first or second: the segments are 5 read bases apart, allowed by
	chimSegmentReadGapMax 5 and not by 4 or the default 0. The reported
	chimera joins A's contig and strand to B's."""
	r = hd.strict_chimeric_read(csim, ("+", "-"), lengths)
	pa, pb = r.extra["pieces"]
	_, seq_a, after_a = hd.piece_context(csim.g, pa, 10)
	before_b, seq_b, _ = hd.piece_context(csim.g, pb, 10)
	gap = "".join(next(x for x in "ACGT" if x not in (after_a[k], before_b[k + 5])) for k in range(5))
	seq = seq_a + gap + seq_b
	reads = ReadSet([sh.SimRead(r.name, seq, csim.qual(len(seq)), kind="chimeric")])
	for gmax, n in ((0, 0), (4, 0), (5, 1)):
		rows = _junction_rows(run_star(genome_plain, reads, "--chimSegmentMin 20 --chimSegmentReadGapMax %d" % gmax))
		assert len(rows) == n, gmax
	row = rows[0]
	assert (row["chr_donorA"], row["strand_donorA"], row["chr_acceptorB"], row["strand_acceptorB"]) == (
		pa[0], pa[2], pb[0], pb[2])


@pytest.mark.parametrize("scheme", ["old", "new"])
@pytest.mark.parametrize("layout", ["opposite_strands", "same_strand_far", "backsplice"])
def test_intrachromosomal_chimeras(genome_plain, csim, layout, scheme):
	"""Pieces on one contig >= 5 kb apart (manual 6: chimeric segments are
	on different strands, or far from each other). On opposite strands they
	are chimeric with default settings. On the same strand, in genomic
	order, only once --alignIntronMax/--alignMatesGapMax 1000 forbid
	joining them linearly. On the same strand with piece B upstream of
	piece A (a back-splice) no linear alignment exists, and with the
	default window size (2^winBinNbits x winAnchorDistNbins = 589,824
	bases, parametersDefault) both segments lie in one alignment window,
	whereas the other chimeras of this module span two windows."""
	if layout == "opposite_strands":
		r = hd.strict_chimeric_read(csim, ("-", "+"), same_contig=True)
	elif layout == "same_strand_far":
		r = hd.strict_chimeric_read(csim, ("+", "+"), same_contig=True)
	else:
		while True:
			r = hd.strict_chimeric_read(csim, ("+", "+"), same_contig=True)
			(_, pa, _, na), (_, pb, _, nb) = r.extra["pieces"]
			if pb + nb + 5000 <= pa:
				break
	extra = ["--alignIntronMax", "1000", "--alignMatesGapMax", "1000"] if layout == "same_strand_far" else []
	extra += NEW if scheme == "new" else []
	rows = run_star(genome_plain, ReadSet([r]), ["--chimSegmentMin", "20"] + extra).chimeric()
	assert [_cols(x) for x in rows] == [list(r.extra["junction"].values())]


def test_spliced_donor_segment(genome_plain, csim):
	"""Old scheme, segment A itself spliced: 50 bases before and 30 after
	the T01a intron chr1:1301-2300 (chr1:1251-1300 and 2301-2330, +
	strand), followed by 40 bases B at a unique gene-free locus on another
	contig, with a clean breakpoint after 2330 (120-base read). Columns by
	hand (manual 6.4): donor breakpoint 2331, the base after A's last
	block; start_alnA 1251 and cigar_alnA 50M1000N30M40S; B's columns as
	for an unspliced chimera with an 80-base A."""
	g = csim.g
	a_last = ("chr1", 2301, "+", 30)
	while True:
		cb, pb = csim.unique_locus(60, avoid_genes=True)
		b = (cb, pb + 10, csim.rng.choice("+-"), 40)
		if cb != "chr1" and hd.junction_is_clean(csim, a_last, b):
			break
	seq = g.seq("chr1", 1251, 1300) + g.seq("chr1", 2301, 2330) + hd.piece_context(g, b, 0)[1]
	read = sh.SimRead("spliced_a", seq, csim.qual(len(seq)), kind="chimeric")
	want = csim.chimeric_read(a_last, b, name="spliced_a").extra["junction"]
	want["start_alnA"], want["cigar_alnA"] = 1251, "50M1000N30M40S"
	want["cigar_alnB"] = "80S40M" if b[2] == "+" else "40M80S"
	assert (want["brkpt_donorA"], want["junction_type"], want["repeat_left_lenA"], want["repeat_right_lenB"]) == (2331, 0, 0, 0)
	rows = run_star(genome_plain, ReadSet([read]), ["--chimSegmentMin", "20"]).chimeric()
	assert [_cols(x) for x in rows] == [list(want.values())]


@pytest.fixture(scope="module")
def chim_pairs(synth_genome):
	"""Pairs with mate 1 at locus A and mate 2 at locus B (different
	contigs, all strand pairs): the junction lies between the mates."""
	csim = _fixture_sim(synth_genome, "chim_pairs")
	pairs = []
	for sa, sb in STRANDS:
		while True:
			ca, pa = csim.unique_locus(120, avoid_genes=True)
			cb, pb = csim.unique_locus(120, avoid_genes=True)
			if ca != cb:
				break
		pairs.append(hd.chimeric_pair(csim, (ca, pa + 10, sa), (cb, pb + 10, sb)))
	return ReadSet(pairs)


@pytest.mark.parametrize("scheme", ["old", "new"])
def test_pe_encompassing_junctions(genome_plain, chim_pairs, scheme):
	"""Manual 6.4, encompassing junctions: type -1, column 2 the base after
	mate 1's 3' end, column 5 the base before mate 2's 5' end (mate 2
	reverse-complemented), repeats 0, each segment's start and CIGAR 100M.
	Multimapping scheme: max_poss_aln_score 200, non_chim 98 (one mate),
	chimeric score 98+98 = 196 (no junction penalty: the junction is not
	sequenced)."""
	run = run_star(genome_plain, chim_pairs, ["--chimSegmentMin", "20"] + (NEW if scheme == "new" else []))
	rows = run.chimeric()
	assert sorted(r["read_name"] for r in rows) == sorted(chim_pairs.by_name)
	for row in rows:
		p = chim_pairs.by_name[row["read_name"]]
		(ca, pa, sa), (cb, pb, sb) = p.extra["pieces"]
		want = list(p.extra["junction"].values()) + [pa, "100M", pb, "100M"]
		assert _cols(row) == want
		if scheme == "new":
			assert _extra(row) == [1, 200, 98, 196, 196, 0]


@pytest.fixture(scope="module")
def merged_pair(synth_genome):
	"""Fragment A(75)+B(75) read as 100-base mates that overlap by 50:
	mate 1 = A + B[:25], mate 2 = reverse complement of A[50:] + B."""
	csim = _fixture_sim(synth_genome, "merged_pair")
	r = hd.strict_chimeric_read(csim, ("+", "-"), (75, 75))
	name = "merged1"
	m1 = sh.SimRead(name, r.seq[:100], csim.qual(100))
	m2 = sh.SimRead(name, sh.revcomp(r.seq[50:]), csim.qual(100))
	m1.mate, m2.mate = 1, 2
	return r, ReadSet([sh.SimPair(name, m1, m2)])


def test_pe_overlap_merged_chimera(genome_plain, merged_pair):
	"""Manual 10: with --peOverlapNbasesMin the overlapping mates are merged
	and chimeric detection runs on the merged sequence. Multimapping scheme:
	one line whose columns 1-14 are those of the 150-base single-end
	chimera, PEmerged_bool 1, max_poss_aln_score 150, non_chim 73 (75
	bases), chimeric score 73+73-1 = 145."""
	r, reads = merged_pair
	run = run_star(genome_plain, reads, ["--chimSegmentMin", "20", "--peOverlapNbasesMin", "20"] + NEW)
	rows = run.chimeric()
	assert len(rows) == 1
	want = list(r.extra["junction"].values())
	want[sh.CHIMERIC_FIELDS.index("read_name")] = "merged1"
	assert _cols(rows[0]) == want
	assert _extra(rows[0]) == [1, 150, 73, 145, 145, 1]


def test_pe_overlap_merged_chimera_within_bam(genome_plain, merged_pair):
	"""Old scheme with --peOverlapNbasesMin and --chimOutType WithinBAM: the
	merged chimera is written back as a pair: mate 1's A part as a
	supplementary, hard-clipped record at A (75 aligned bases), mate 1's B
	part (25 bases) and mate 2 (75 bases, reverse strand) at B's start."""
	r, reads = merged_pair
	(ca, pa, sa, na), (cb, pb, sb, nb) = r.extra["pieces"]
	run = run_star(genome_plain, reads, "--chimSegmentMin 20 --peOverlapNbasesMin 20 --chimOutType WithinBAM "
		"--outSAMtype BAM Unsorted")
	recs = run.records()
	sup = [x for x in recs if x.is_supplementary]
	assert len(recs) == 3 and len(sup) == 1
	s = sup[0]
	assert (s.mate, s.rname, s.pos, s.strand, sh.cigar_ref_len(s.cigar)) == (1, ca, pa, sa, 75)
	assert "H" in s.cigar
	rest = sorted((x.mate, x.rname, x.strand, sh.cigar_ref_len(x.cigar)) for x in recs if not x.is_supplementary)
	assert rest == [(1, cb, sb, 25), (2, cb, "+" if sb == "-" else "-", 75)]


def test_junction_format_comment_lines(genome_plain, clean_reads, csim):
	"""--chimOutJunctionFormat 1 appends '# <version>   <command line>' and
	'# Nreads N\\tNreadsUnique U\\tNreadsMulti M' (parametersDefault); with
	4 chimeric reads (not mapped linearly), 3 unique reads and one read
	from the 3-copy repeat dup3 these are 8, 3 and 1. Format 0 (default)
	writes no comment lines."""
	extra = [csim.genomic_read(*csim.unique_locus(100, avoid_genes=True), name="lin%d" % i) for i in range(3)]
	extra.append(csim.segment_read("dup3", 200, name="multi1"))
	reads = ReadSet(clean_reads.records + extra)
	run = run_star(genome_plain, reads, "--chimSegmentMin 20 --chimOutJunctionFormat 1")
	lines = run.path("Chimeric.out.junction").read_text().splitlines()
	assert [l for l in lines if not l.startswith("#")] and len(run.chimeric()) == 4
	comments = [l for l in lines if l.startswith("#")]
	assert lines[-2:] == comments
	version = subprocess.run([str(sh.find_star()), "--version"], stdout=subprocess.PIPE).stdout.decode().strip()
	assert comments[0].startswith("# %s   " % version) and comments[0].endswith("--chimOutJunctionFormat 1")
	assert comments[1] == "# Nreads 8\tNreadsUnique 3\tNreadsMulti 1"
	plain = run_star(genome_plain, reads, "--chimSegmentMin 20")
	assert not any(l.startswith("#") for l in plain.path("Chimeric.out.junction").read_text().splitlines())


def test_separate_sam_old(genome_plain, clean_reads):
	"""--chimOutType Junctions SeparateSAMold writes Chimeric.out.junction
	and Chimeric.out.sam, the latter with both segments of each chimera:
	the positions, strands and CIGARs of columns 11-14, NH 2, HI 1 and 2,
	exactly one primary record per read; Aligned.out.sam has none."""
	run = run_star(genome_plain, clean_reads, "--chimSegmentMin 20 --chimOutType Junctions SeparateSAMold")
	assert len(run.chimeric()) == len(clean_reads)
	by_name = run.alignments("Chimeric.out.sam").by_name()
	assert sorted(by_name) == sorted(clean_reads.by_name)
	for name, recs in by_name.items():
		j = clean_reads.by_name[name].extra["junction"]
		(ca, pa, sa, _), (cb, pb, sb, _) = clean_reads.by_name[name].extra["pieces"]
		assert sorted((x.rname, x.pos, x.strand, x.cigar) for x in recs) == sorted([
			(ca, j["start_alnA"], sa, j["cigar_alnA"]), (cb, j["start_alnB"], sb, j["cigar_alnB"])])
		assert sorted((x.tag("NH"), x.tag("HI")) for x in recs) == [(2, 1), (2, 2)]
		assert sum(not x.is_secondary for x in recs) == 1
	assert run.records() == [] or all(x.is_unmapped for x in run.records())


def test_within_bam_soft_clip(genome_plain, clean_reads):
	"""--chimOutType WithinBAM SoftClip: as with HardClip, two records per
	chimeric read at the two segments, one supplementary, but the
	supplementary record keeps the whole read soft-clipped (no H, 100
	bases of SEQ)."""
	run = run_star(genome_plain, clean_reads, "--chimSegmentMin 20 --chimOutType WithinBAM SoftClip "
		"--outSAMtype BAM Unsorted")
	by_name = run.alignments().by_name()
	assert sorted(by_name) == sorted(clean_reads.by_name)
	for name, recs in by_name.items():
		j = clean_reads.by_name[name].extra["junction"]
		assert sorted((x.pos, x.cigar) for x in recs) == sorted([(j["start_alnA"], j["cigar_alnA"]),
			(j["start_alnB"], j["cigar_alnB"])])
		sup = [x for x in recs if x.is_supplementary]
		assert len(sup) == 1 and "H" not in sup[0].cigar and len(sup[0].seq) == 100


@pytest.fixture(scope="module")
def n_genome(genome_factory, synth_genome, tmp_path_factory):
	"""The synthetic genome with chr3:17001-17010 replaced by N."""
	synth = sh.make_synthetic_genome()
	s = synth.contigs["chr3"]
	synth.contigs["chr3"] = s[:17000] + "N" * 10 + s[17010:]
	synth.write(tmp_path_factory.mktemp("synth_n"))
	return synth, genome_factory(gtf=False, synth=synth)


def test_chim_filter_genomic_n(n_genome):
	"""Piece A ends 20 bases before the N run, so the genome continuing A
	inside the read's B part has Ns. parametersDefault chimFilter:
	banGenomicN (default) "Ns are not allowed in the genome sequence around
	the chimeric junction" rejects it; None reports it."""
	synth, genome = n_genome
	sim = sh.ReadSimulator(synth, seed=17)
	a = ("chr3", 16931, "+", 50)
	r = hd.strict_chimeric_read(sim, ("+", "-"), locus_a=lambda: a[:2])
	assert "N" in hd.piece_context(synth, r.extra["pieces"][0], 50)[2]
	reads = ReadSet([r])
	assert _junction_rows(run_star(genome, reads, "--chimSegmentMin 20")) == []
	rows = run_star(genome, reads, "--chimSegmentMin 20 --chimFilter None").chimeric()
	assert [_cols(x) for x in rows] == [list(r.extra["junction"].values())]


@pytest.mark.parametrize("scheme", ["old", "new"])
def test_read_group_column(genome_plain, clean_reads, scheme):
	"""Two input files with one read group each (--outSAMattrRGline ID:grpA
	, ID:grpB) and RG in --outSAMattributes: every Chimeric.out.junction
	line ends with its read's group, as column 21 'readgrp' after the six
	score columns in the multimapping scheme (manual 6.4) and right after
	column 14 in the old scheme (pins the current output)."""
	a, b = ReadSet(clean_reads.records[:2]), ReadSet(clean_reads.records[2:])
	args = ["--chimSegmentMin", "20", "--outSAMattrRGline", "ID:grpA", ",", "ID:grpB",
		"--outSAMattributes", "NH", "HI", "AS", "nM", "RG"] + (NEW if scheme == "new" else [])
	rows = run_star(genome_plain, [a, b], args).chimeric()
	assert sorted(r["read_name"] for r in rows) == sorted(clean_reads.by_name)
	for row in rows:
		p = clean_reads.by_name[row["read_name"]]
		grp = "grpA" if p.name in a.by_name else "grpB"
		assert _cols(row) == list(p.extra["junction"].values())
		assert row["extra"] == (["1", "100", "49", "97", "97", "0", grp] if scheme == "new" else [grp])


# limitIObufferSize[0] - 2 * 651 - 2 * 50000 = 8698 bytes per input chunk
# (Parameters.cpp, IncludeDefine.h): the reads below make several chunks.
SMALL_CHUNKS = ["--limitIObufferSize", "110000", "50000000"]


@pytest.mark.parametrize("scheme", ["old", "new"])
def test_chimeric_output_does_not_depend_on_input_chunks(genome_plain, csim, scheme):
	"""24 clean chimeras among 150 transcript reads, in random order, mapped
	as several input chunks by 4 threads (SMALL_CHUNKS): the
	Chimeric.out.junction lines are those of one chunk mapped by 1 thread
	(in any order), columns 1-14 of each equal the truth, and the
	'# Nreads' line (--chimOutJunctionFormat 1) and the Log.final.out
	chimeric count are the same."""
	chim = [csim.clean_chimeric_read(s) for s in STRANDS for _ in range(6)]
	reads = chim + list(csim.transcript_reads(150))
	csim.rng.shuffle(reads)
	args = ["--chimSegmentMin", "20", "--chimOutJunctionFormat", "1"] + (NEW if scheme == "new" else [])
	one = run_star(genome_plain, ReadSet(reads), args, threads=1)
	four = run_star(genome_plain, ReadSet(reads), args + SMALL_CHUNKS, threads=4)
	lines = lambda run: [l for l in run.path("Chimeric.out.junction").read_text().splitlines() if not l.startswith("#")]
	assert sorted(lines(four)) == sorted(lines(one))
	truth = {r.name: list(r.extra["junction"].values()) for r in chim}
	assert {r["read_name"]: _cols(r) for r in one.chimeric()} == truth
	nreads = lambda run: [l for l in run.path("Chimeric.out.junction").read_text().splitlines() if l.startswith("# Nreads")]
	assert nreads(four) == nreads(one) and len(nreads(one)) == 1
	assert one.log_final()["Number of chimeric reads"] == four.log_final()["Number of chimeric reads"] == len(chim)
