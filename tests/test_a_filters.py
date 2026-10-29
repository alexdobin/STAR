"""Area A, output filters (ReadAlign_mappedFilter.cpp, stitchWindowAligns.cpp):
the Log.final.out unmapped categories, --outFilterScoreMin(OverLread),
--outFilterMatchNmin(OverLread), the mismatch filters,
--outFilterIntronMotifs and --outFilterIntronStrands.

The OverLread filters compare with (int)(ratio * (Lread - 1)); for
single-end reads parametersDefault says "normalized to read length", which
differs by one base, so
single-end boundaries here use only scores on which both readings agree."""

import pytest

import helpers_a as ha
import star_helpers as sh
from star_helpers import ReadSet, SimPair, record_truth, run_star

ATTRS = ["--outSAMattributes", "NH", "HI", "AS", "nM"]


def _pct(n, total):
	return round(100.0 * n / total, 2)


def test_unmapped_categories(genome_plain, sim):
	"""One run (--alignEndsType EndToEnd --outFilterMismatchNmax 3) with
	reads of every outcome, counted in Log.final.out as in mappedFilter():
	8 unique and 4 three-copy reads map; 6 reads with an N every 11 bases
	have no N-free piece of seedSplitMin bases, hence no window: 'other'; 6
	random reads are extended end to end over random bases, so their score
	is negative, clamped to 0 (stitchWindowAligns.cpp) and a best score of
	0 means no window (stitchPieces.cpp): 'other'; 6 reads of 50 Ns + 50
	genome bases score 48 < int(0.66 * 99): 'too short'; 6 reads with
	mismatches at offsets 1, 2, 97, 98 are forced end to end with 4
	mismatches > 3: 'too many mismatches'; 6 reads of the 12-copy segment
	exceed outFilterMultimapNmax 10: 'too many loci'. Each unmapped record
	(--outSAMunmapped Within) carries the matching uT code (manual: 0 other,
	1 too short, 2 too many mismatches, 3 too many loci). Percentages are of
	42 input reads."""
	reads = [sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2], name="uniq%d" % i) for i in range(8)]
	reads += [sim.segment_read("dup3", sim.rng.randint(0, 500), strand="+-"[i % 2], name="multi%d" % i)
		for i in range(4)]
	reads += [sim.genomic_read(*sim.unique_locus(), n_offsets=tuple(range(10, 100, 11)), name="nopiece%d" % i)
		for i in range(6)]
	reads += [sim.random_read(name="random%d" % i) for i in range(6)]
	for i in range(6):
		c, p = sim.unique_locus()
		reads.append(ha.n_head_read("nhead%d" % i, sim.g, c, p, 50, "+-"[i % 2]))
	reads += [sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2], mismatches=(1, 2, 97, 98), name="mm%d" % i)
		for i in range(6)]
	reads += [sim.segment_read("rep12", sim.rng.randint(0, 200), name="rep%d" % i) for i in range(6)]
	reads = ReadSet(reads)
	run = run_star(genome_plain, reads, ATTRS + ["--alignEndsType", "EndToEnd", "--outFilterMismatchNmax", "3",
		"--outSAMunmapped", "Within"])
	aln = run.alignments()
	reason = {"nopiece": "0", "random": "0", "nhead": "1", "mm": "2", "rep": "3"}
	unmapped = [r for r in aln if r.is_unmapped]
	assert sorted(r.qname for r in unmapped) == sorted(n for n in reads.by_name if not n.startswith(("uniq", "multi")))
	for rec in unmapped:
		assert rec.tag("uT") == reason[rec.qname.rstrip("0123456789")], rec.qname
	by_name = sh.AlignmentFile(aln.header, aln.references, aln.mapped()).by_name()
	assert sorted(by_name) == sorted(n for n in reads.by_name if n.startswith(("uniq", "multi")))
	for name, recs in by_name.items():
		r = reads.by_name[name]
		if name.startswith("uniq"):
			assert [record_truth(x) for x in recs] == [r.truth()]
		else:
			assert sorted((x.rname, x.pos, x.strand) for x in recs) == sorted(r.extra["loci"])
	log = run.log_final()
	expect = {
		"Number of input reads": 42,
		"Uniquely mapped reads number": 8,
		"Uniquely mapped reads %": _pct(8, 42),
		"Number of reads mapped to multiple loci": 4,
		"% of reads mapped to multiple loci": _pct(4, 42),
		"Number of reads mapped to too many loci": 6,
		"% of reads mapped to too many loci": _pct(6, 42),
		"Number of reads unmapped: too many mismatches": 6,
		"% of reads unmapped: too many mismatches": _pct(6, 42),
		"Number of reads unmapped: too short": 6,
		"% of reads unmapped: too short": _pct(6, 42),
		"Number of reads unmapped: other": 12,
		"% of reads unmapped: other": _pct(12, 42),
	}
	assert {k: log[k] for k in expect} == expect


def _n_head_reads(sim, ks, n=2):
	reads = []
	for k in ks:
		for i in range(n):
			c, p = sim.unique_locus(130)
			reads.append(ha.n_head_read("k%02d_%d" % (k, i), sim.g, c, p + k, k, "+-"[i % 2]))
	return ReadSet(reads)


# Each case: extra args and {k: mapped}. Reads are k Ns + (100 - k) genome
# bases, so nMatch = 100 - k and AS = 100 - k - 2 (length term; 98 - k).
FILTER_CASES = [
	("default", [], {32: True, 34: False}),
	("ScoreMin78", ["--outFilterScoreMin", "78"], {20: True, 21: False}),
	("MatchNmin80", ["--outFilterMatchNmin", "80"], {20: True, 21: False}),
	("ScoreMinOverLread0.8", ["--outFilterScoreMinOverLread", "0.8"], {18: True, 20: False}),
	("MatchNminOverLread0.8", ["--outFilterMatchNminOverLread", "0.8"], {20: True, 22: False}),
]


@pytest.mark.parametrize("case", FILTER_CASES, ids=[c[0] for c in FILTER_CASES])
def test_score_and_match_filters(genome_plain, sim, case):
	"""Reads with k leading Ns align as kS(100-k)M (Ns never seed or score).
	default 0.66/0.66: AS 66 (k 32) passes, AS 64 (k 34) is too short.
	--outFilterScoreMin 78: AS >= 78 passes (k 20), 77 fails (k 21).
	--outFilterMatchNmin 80: 80 matched bases pass, 79 fail.
	--outFilterScoreMinOverLread 0.8: AS 80 passes, 78 fails.
	--outFilterMatchNminOverLread 0.8: 80 matched bases pass, 78 fail.
	In each non-default case the other ratio filters are set to 0 so only the
	named one applies. Failing reads count as 'too short'."""
	name, args, ks = case
	if name != "default":
		zero = {"--outFilterScoreMinOverLread": "0", "--outFilterMatchNminOverLread": "0"}
		for a in args[::2]:
			zero.pop(a, None)
		args = args + [x for kv in zero.items() for x in kv]
	reads = _n_head_reads(sim, sorted(ks))
	run = run_star(genome_plain, reads, ATTRS + args)
	got = {r.qname: r for r in run.alignments()}
	for r in reads:
		k = int(r.name[1:3])
		if ks[k]:
			assert record_truth(got[r.name]) == r.truth()
			assert got[r.name].tag("AS") == 98 - k
		else:
			assert r.name not in got
	assert ha.unmapped_counts(run.log_final()) == (0, sum(not v for v in ks.values()) * 2, 0, 0)


def test_score_filter_paired_end_sum_of_mates(genome_plain, sim):
	"""For pairs the ratio filters use the sum of the mates' lengths
	(parametersDefault): int(0.66 * 200) = 132. Mate 2 carries k Ns at its
	outer end, so AS = 200 - k - 2 (length term of a 300-base span): k 66
	gives 132 and maps, k 67 gives 131 and is too short."""
	pairs = []
	for i, k in enumerate((66, 66, 67, 67)):
		c, p = sim.unique_locus(400)
		name = "pk%d_%d" % (k, i)
		r1 = sim.genomic_read(c, p, "+", name=name)
		fwd = sim.g.seq(c, p + 200, p + 299 - k) + "N" * k
		r2 = ha.oriented(name, fwd, c, p + 200, "%dM%dS" % (100 - k, k), "-")
		r1.mate, r2.mate = 1, 2
		pairs.append(SimPair(name, r1, r2))
	pairs = ReadSet(pairs)
	run = run_star(genome_plain, pairs, ATTRS)
	by_name = run.alignments().by_name()
	assert sorted(by_name) == sorted(n for n in pairs.by_name if n.startswith("pk66"))
	for name, recs in by_name.items():
		p = pairs.by_name[name]
		m = {r.mate: r for r in recs}
		assert record_truth(m[1]) == p.r1.truth() and record_truth(m[2]) == p.r2.truth()
		assert m[1].tag("AS") == 132
	assert ha.unmapped_counts(run.log_final()) == (0, 2, 0, 0)


def _central_mm_reads(sim, ks, n=2):
	"""Reads with k mismatches at offsets 36, 40, ... (inside 36..60), so a
	partial alignment that avoids them is too short."""
	reads = []
	for k in ks:
		for i in range(n):
			reads.append(sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2],
				mismatches=tuple(36 + 4 * j for j in range(k)), name="mm%d_%d" % (k, i)))
	return ReadSet(reads)


@pytest.mark.parametrize("args,limit", [
	(["--outFilterMismatchNmax", "3"], 3),
	(["--outFilterMismatchNoverLmax", "0.04"], 4),
	(["--outFilterMismatchNoverReadLmax", "0.05"], 5),
])
def test_mismatch_filters(genome_plain, sim, args, limit):
	"""Reads with 0-6 mismatches in their centre. A read is output only if
	its mismatches are within --outFilterMismatchNmax, within
	--outFilterMismatchNoverLmax x mapped length (100) and within
	--outFilterMismatchNoverReadLmax x read length. STAR extends seeds over
	mismatches only up to these limits (extendAlign.cpp), so a read over the
	limit keeps a partial alignment below 0.66 of the read and is counted as
	'too short'. Reads within the limit align end to end with nM = k."""
	reads = _central_mm_reads(sim, range(7))
	run = run_star(genome_plain, reads, ATTRS + args)
	got = {r.qname: r for r in run.alignments()}
	for r in reads:
		k = len(r.mismatches)
		if k <= limit:
			assert record_truth(got[r.name]) == r.truth() and got[r.name].tag("nM") == k
		else:
			assert r.name not in got
	assert ha.unmapped_counts(run.log_final()) == (0, 2 * (6 - limit), 0, 0)


def test_mismatch_budget_paired_end(genome_plain, sim):
	"""--outFilterMismatchNoverReadLmax 0.05 gives a pair of 100-base mates a
	budget of int(0.05 * 200) = 10 mismatches (ReadAlign_oneRead.cpp), so
	pairs with 4 central mismatches in each mate align end to end with nM 8."""
	pairs = ReadSet([sim.genomic_pair(*sim.unique_locus(400), insert=300, strand="+-"[i % 2],
		mismatches1=(36, 44, 52, 60), mismatches2=(38, 46, 54, 62)) for i in range(6)])
	aln = run_star(genome_plain, pairs, ATTRS + ["--outFilterMismatchNoverReadLmax", "0.05"]).alignments()
	by_name = aln.by_name()
	assert sorted(by_name) == sorted(pairs.by_name)
	for name, recs in by_name.items():
		p = pairs.by_name[name]
		m = {r.mate: r for r in recs}
		assert record_truth(m[1]) == p.r1.truth() and record_truth(m[2]) == p.r2.truth()
		assert m[1].tag("nM") == 8


def _motif_filter_reads(sim):
	reads = []
	for left in (40, 50, 60):
		reads.append(ha.tx_junction_read(sim, "T12nc", 0, left, sense=left != 50, name="nc%d" % left))
	for tx in ("T01a", "T03", "T04"):
		reads.append(ha.tx_junction_read(sim, tx, 0, 50, name="can_" + tx))
	reads += [sim.genomic_read(*sim.unique_locus(), name="unspl%d" % i) for i in range(3)]
	return ReadSet(reads)


def test_intron_motif_filter(genome_plain, sim, tmp_path):
	"""--outFilterIntronMotifs RemoveNoncanonical drops alignments with a
	non-canonical junction; RemoveNoncanonicalUnannotated drops them only
	when the junction is not annotated (parametersDefault). The T12nc
	junction (non-canonical) is unannotated in genome_plain and annotated in
	the second pair of runs with --sjdbFileChrStartEnd. A dropped T12nc read
	has no other alignment within 0.66 of its length (unmapped, too short);
	canonical and unspliced reads are unaffected."""
	reads = _motif_filter_reads(sim)
	for r in reads:
		if r.name.startswith("nc"):
			(s1, e1), (s2, e2) = sh.cigar_blocks(r.pos, r.cigar)
			assert ha.unspliced_fallback(sim.g.contigs[r.contig], r.sam_seq(), s1, e1 - s1 + 1, s2)[0][0] < 65
	sj = tmp_path / "t12nc.sj"
	sj.write_text("chr3\t15201\t15800\t+\n")
	annot = ["--sjdbFileChrStartEnd", str(sj), "--sjdbOverhang", "99"]
	for args, dropped in (
			(["--outFilterIntronMotifs", "RemoveNoncanonical"], True),
			(["--outFilterIntronMotifs", "RemoveNoncanonicalUnannotated"], True),
			(["--outFilterIntronMotifs", "RemoveNoncanonical"] + annot, True),
			(["--outFilterIntronMotifs", "RemoveNoncanonicalUnannotated"] + annot, False)):
		run = run_star(genome_plain, reads, ATTRS + args)
		got = run.alignments().by_name()
		want = [n for n in reads.by_name if not (dropped and n.startswith("nc"))]
		assert sorted(got) == sorted(want), args
		for name in want:
			assert [record_truth(x) for x in got[name]] == [reads.by_name[name].truth()]
		assert ha.unmapped_counts(run.log_final()) == (0, 3 if dropped else 0, 0, 0)


def _two_junction_read(sim, name, motifs, blocks=(36, 28, 36)):
	"""Unspliced-genome read of three blocks joined across two de novo
	introns (200-1000 bases) whose + strand motifs are `motifs`, e.g.
	(("GT", "AG"), ("CT", "AC")); junctions cannot slide and the scan keeps
	them in place."""
	g = sim.g
	b1, b2, b3 = blocks
	while True:
		c, p = sim.unique_locus(60, avoid_genes=True)
		seq = g.contigs[c]
		pos = [p]
		ok = True
		for j, (don, acc) in enumerate(motifs):
			d = pos[-1] + (b1 if j == 0 else b2) - 1
			if seq[d:d + 2] != don:
				ok = False
				break
			cands = [a for a in range(d + 201, min(d + 1001, len(seq) - 200))
				if seq[a - 3:a - 1] == acc and seq[d - 1] != seq[a - 2] and seq[d] != seq[a - 1]]
			if not cands:
				ok = False
				break
			pos.append(sim.rng.choice(cands))
		if not ok:
			continue
		(p1, a1, a2) = pos
		if any(_overlaps_segment(g, c, s, s + 40) for s in (a1, a2)):
			continue
		fwd = seq[p1 - 1:p1 - 1 + b1] + seq[a1 - 1:a1 - 1 + b2] + seq[a2 - 1:a2 - 1 + b3]
		d1, d2 = p1 + b1 - 1, a1 + b2 - 1
		if ha.junction_placement(seq, fwd[:b1 + b2], b1, d1, a1)[0] != 0 or \
				ha.junction_placement(seq, fwd[b1:], b2, d2, a2)[0] != 0:
			continue
		cig = "%dM%dN%dM%dN%dM" % (b1, a1 - d1 - 1, b2, a2 - d2 - 1, b3)
		return ha.oriented(name, fwd, c, p1, cig, sim.rng.choice("+-"))


def _overlaps_segment(synth, contig, start, end):
	return any(k.contig == contig and not (end < k.start or start > k.end)
		for seg in synth.segments.values() for k in seg.copies)


def test_intron_strand_filter(genome_plain, sim):
	"""A read spliced at a + strand junction (GT/AG) and a - strand junction
	(CT/AC) has inconsistent intron strands: the default
	--outFilterIntronStrands RemoveInconsistentStrands drops the alignment
	(the best partial alignment, 64 bases with one junction, is too short);
	None keeps it. Reads whose two junctions are both GT/AG align in both
	runs."""
	bad = [_two_junction_read(sim, "mixed%d" % i, (("GT", "AG"), ("CT", "AC"))) for i in range(3)]
	bad += [_two_junction_read(sim, "mixedr%d" % i, (("CT", "AC"), ("GT", "AG"))) for i in range(2)]
	good = [_two_junction_read(sim, "same%d" % i, (("GT", "AG"), ("GT", "AG"))) for i in range(3)]
	reads = ReadSet(bad + good)
	run = run_star(genome_plain, reads, ATTRS)
	got = run.alignments().by_name()
	assert sorted(got) == sorted(r.name for r in good)
	for name, recs in got.items():
		assert [record_truth(x) for x in recs] == [reads.by_name[name].truth()]
	assert ha.unmapped_counts(run.log_final()) == (0, len(bad), 0, 0)
	got = run_star(genome_plain, reads, ATTRS + ["--outFilterIntronStrands", "None"]).alignments().by_name()
	assert sorted(got) == sorted(reads.by_name)
	for name, recs in got.items():
		assert [record_truth(x) for x in recs] == [reads.by_name[name].truth()]
