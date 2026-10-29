"""Area A, alignment scores: junction motif penalties (--scoreGap*), indel
penalties, --scoreGenomicLengthLog2scale, --scoreStitchSJshift,
--alignSJstitchMismatchNmax and --alignInsertionFlush. AS values are worked
from the scoring rules: +1 per matched base, -1 per mismatch, the gap
penalties from parametersDefault and the length term in helpers_a.length_score.
Junction positions come from helpers_a.junction_placement, a transcription of
the junction scan in stitchAlignToTranscript.cpp."""

import pytest

import helpers_a as ha
import star_helpers as sh
from star_helpers import ReadSet, record_truth, run_star

ATTRS = ["--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD"]

# (transcript, intron index): every motif class, both strands, introns of
# 400 to 5000 bases
MOTIF_INTRONS = [("T01a", 0), ("T02", 0), ("T03", 0), ("T04", 0), ("T08", 0), ("T11", 0),
	("T12nc", 0), ("T09", 0), ("T10", 1)]

CUSTOM = {"gap": -1, "Noncan": -10, "GCAG": -3, "ATAC": -9}
PEN_KEY = {0: "Noncan", 3: "GCAG", 4: "GCAG", 5: "ATAC", 6: "ATAC"}


def _motif_reads(sim):
	reads = []
	for tx, k in MOTIF_INTRONS:
		for left in (40, 50, 60):
			reads.append(ha.tx_junction_read(sim, tx, k, left, sense=left != 50, name="%s_%d" % (tx, left)))
	return ReadSet(reads)


def _junction(synth, read):
	"""(donor end, acceptor start, motif code, intron length) of a one-intron read."""
	(s1, e1), (s2, e2) = sh.cigar_blocks(read.pos, read.cigar)
	g = synth.contigs[read.contig]
	code, _ = ha.motif_and_penalty(g[e1:e1 + 2], g[s2 - 3:s2 - 1], ha.DEFAULT_PENS)
	return e1, s2, code, s2 - e1 - 1


@pytest.mark.parametrize("custom", [False, True])
def test_junction_motif_penalties(genome_plain, sim, custom):
	"""Spliced reads over every intron motif, without annotation. AS = 100 +
	length term of the span (100 + intron) + scoreGap + the motif penalty:
	0 for GT/AG and CT/AC, scoreGapGCAG for GC/AG and CT/GC, scoreGapATAC for
	AT/AC and GT/AT, scoreGapNoncan otherwise. Defaults 0/-4/-8/-8 give e.g.
	97 for the 1000-base GT/AG intron of T01a and 90 for the 600-base AT/AC
	intron of T04. The custom run uses scoreGap -1, Noncan -10, GCAG -3,
	ATAC -9, chosen so the junction scan still keeps every junction in place
	(checked with junction_placement)."""
	reads = _motif_reads(sim)
	pens = {k: CUSTOM[k] for k in ("Noncan", "GCAG", "ATAC")} if custom else ha.DEFAULT_PENS
	gap = CUSTOM["gap"] if custom else 0
	args = ATTRS + (["--scoreGap", str(gap), "--scoreGapNoncan", str(pens["Noncan"]),
		"--scoreGapGCAG", str(pens["GCAG"]), "--scoreGapATAC", str(pens["ATAC"])] if custom else [])
	aln = run_star(genome_plain, reads, args).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	codes = set()
	for rec in aln:
		r = reads.by_name[rec.qname]
		d, a, code, n = _junction(sim.g, r)
		left = sh.cigar_ops(r.cigar)[0][0]
		assert ha.junction_placement(sim.g.contigs[r.contig], r.sam_seq(), left, d, a, pens=pens) == (0, code, 0)
		assert record_truth(rec) == r.truth()
		pen = pens[PEN_KEY[code]] if code in PEN_KEY else 0
		assert rec.tag("AS") == 100 + ha.length_score(100 + n) + gap + pen, (rec.qname, code)
		codes.add(code)
	assert codes == {0, 1, 2, 3, 4, 5, 6}


def test_indel_penalties(genome_plain, sim):
	"""--scoreDelOpen -3 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -3:
	a 3-base deletion scores 100 - 3 - 3 - 2 (length term) = 92 and a 2-base
	insertion 98 - 1 - 6 - 2 = 89, at the simulated CIGAR (pilot
	test_alignment_score_hand_worked has the defaults)."""
	reads = []
	while len(reads) < 12:
		c, p = sim.unique_locus(120)
		strand = "+-"[len(reads) % 2]
		if len(reads) % 4 < 2:
			if sim.deletion_is_unambiguous(c, p + 45, 3):
				reads.append(sim.genomic_read(c, p, strand, deletion=(45, 3), name="del%d" % len(reads)))
		else:
			reads.append(sim.genomic_read(c, p, strand, insertion=(52, 2), name="ins%d" % len(reads)))
	reads = ReadSet(reads)
	args = ATTRS + ["--scoreDelOpen", "-3", "--scoreDelBase", "-1", "--scoreInsOpen", "-1", "--scoreInsBase", "-3"]
	aln = run_star(genome_plain, reads, args).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		assert record_truth(rec) == reads.by_name[rec.qname].truth()
		assert rec.tag("AS") == (92 if rec.qname.startswith("del") else 89)


@pytest.mark.parametrize("scale", [0, -1])
def test_genomic_length_log2scale(genome_plain, sim, scale):
	"""The length term is int(ceil(log2(span) * scale - 0.5)) and is skipped
	when the scale is 0 (stitchWindowAligns.cpp). Scale 0: 100M scores 100,
	every spliced read 100 + motif penalty. Scale -1: 100M scores 100 - 7 =
	93 (ceil(-6.64 - 0.5)), the 1000-base T01a intron 100 - 10 = 90, the
	5000-base T09 intron 100 - 12 = 88, the 400-base GC/AG T03 intron 100 - 9
	- 4 = 87."""
	reads = ReadSet([sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2]) for i in range(4)]
		+ [ha.tx_junction_read(sim, tx, 0, 50, sense=bool(i), name="%s_%d" % (tx, i))
			for tx in ("T01a", "T03", "T09") for i in range(2)])
	aln = run_star(genome_plain, reads, ATTRS + ["--scoreGenomicLengthLog2scale", str(scale)]).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		r = reads.by_name[rec.qname]
		assert record_truth(rec) == r.truth()
		span = sh.cigar_ref_len(r.cigar)
		pen = -4 if r.tx == "T03" else 0
		assert rec.tag("AS") == 100 + pen + ha.length_score(span, scale)
	if scale == -1:
		assert sorted({r.tag("AS") for r in aln}) == [87, 88, 90, 93]


def _in_segment(synth, contig, start, end):
	return any(c.contig == contig and not (end < c.start or start > c.end)
		for seg in synth.segments.values() for c in seg.copies)


def _shift_junction_reads(sim, n_each=3):
	"""De novo junctions (50 + 50 bases, intron 300-3000) that the scan keeps
	in place as non-canonical with scoreStitchSJshift 1 but moves k bases left
	to a GT/AG or CT/AC motif with scoreStitchSJshift 3, at the cost of 2 (or
	3) read bases that then mismatch."""
	g = sim.g
	reads = []
	count = {2: 0, 3: 0}
	while min(count.values()) < n_each:
		c, p = sim.unique_locus(60, avoid_genes=True)
		d = p + 49
		a = d + 1 + sim.rng.randint(300, 3000)
		seq = g.contigs[c]
		if a + 60 > len(seq) or _in_segment(g, c, a, a + 49) or seq[d - 1] == seq[a - 2] or seq[d] == seq[a - 1]:
			continue
		fwd = seq[p - 1:d] + seq[a - 1:a + 49]
		j1 = ha.junction_placement(seq, fwd, 50, d, a, shift=1)
		j3 = ha.junction_placement(seq, fwd, 50, d, a, shift=3)
		if j1 != (0, 0, 0) or j3[0] >= 0 or j3[1] not in (1, 2) or j3[2] not in count or count[j3[2]] >= n_each:
			continue
		count[j3[2]] += 1
		r = ha.oriented("shift%d" % len(reads), fwd, c, p, "50M%dN50M" % (a - d - 1), "+-"[len(reads) % 2],
			extra={"j3": j3})
		reads.append(r)
	return ReadSet(reads)


def test_stitch_sj_shift_and_mismatch_limit(genome_plain, sim):
	"""--scoreStitchSJshift bounds how far left of the seed boundary the
	junction scan starts (moving left while lost matches minus the shift
	stay >= 0). Default 1: the scan cannot reach the canonical motif k bases
	left, so the true non-canonical junction is kept: 50M<I>N50M, AS = 100 -
	8 + length term. Shift 3: the junction moves to the GT/AG (CT/AC) motif,
	(50-k)M<I>N(50+k)M with the 2 or 3 displaced bases as mismatches, AS =
	100 - 2 * mismatches + length term (motif penalty 0). Adding
	--alignSJstitchMismatchNmax 0 2 0 0 allows at most 2 mismatches when
	stitching GT/AG and CT/AC junctions: the 3-mismatch reads lose the
	junction and their 50-base remainder is too short, so they are unmapped."""
	reads = _shift_junction_reads(sim)
	aln = run_star(genome_plain, reads, ATTRS).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		r = reads.by_name[rec.qname]
		assert record_truth(rec) == r.truth()
		assert rec.tag("AS") == 100 - 8 + ha.length_score(sh.cigar_ref_len(r.cigar)) and rec.tag("nM") == 0
	expect = {}
	for r in reads:
		jr, code, mm = r.extra["j3"]
		n = sh.cigar_ops(r.cigar)[1][0]
		expect[r.name] = (r.contig, r.pos, "%dM%dN%dM" % (50 + jr, n, 50 - jr), r.strand, mm,
			100 - 2 * mm + ha.length_score(100 + n))
	aln = run_star(genome_plain, reads, ATTRS + ["--scoreStitchSJshift", "3"]).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		assert record_truth(rec) + (rec.tag("nM"), rec.tag("AS")) == expect[rec.qname]
	run = run_star(genome_plain, reads, ATTRS + ["--scoreStitchSJshift", "3", "--alignSJstitchMismatchNmax", "0", "2", "0", "0"])
	aln = run.alignments()
	assert sorted(aln.by_name()) == sorted(n for n, e in expect.items() if e[4] == 2)
	for rec in aln:
		assert record_truth(rec) + (rec.tag("nM"), rec.tag("AS")) == expect[rec.qname]
	assert ha.unmapped_counts(run.log_final()) == (0, sum(e[4] == 3 for e in expect.values()), 0, 0)


def test_stitch_mismatch_limit_atac(genome_plain, sim):
	"""With --scoreGapATAC -4 the junction scan prefers, for the
	non-canonical T12nc intron, an AT/AC motif 3 bases to the right at the
	cost of one mismatch (-1 - 4 > -8; junction_placement). The default
	--alignSJstitchMismatchNmax allows 0 mismatches for AT/AC junctions, so
	the stitch is rejected and the reads are unmapped (too short); with
	--alignSJstitchMismatchNmax 0 -1 0 1 they align as 53M600N47M (40/60
	and 60/40 splits shift likewise) with nM 1 and AS = 100 - 2 - 4 - 2 =
	92."""
	reads = ReadSet([ha.tx_junction_read(sim, "T12nc", 0, left, sense=bool(i % 2), name="nc%d_%d" % (left, i))
		for left in (40, 50, 60) for i in range(2)])
	pens = {"Noncan": -8, "GCAG": -4, "ATAC": -4}
	for r in reads:
		d, a, code, n = _junction(sim.g, r)
		left = sh.cigar_ops(r.cigar)[0][0]
		assert ha.junction_placement(sim.g.contigs[r.contig], r.sam_seq(), left, d, a, pens=pens) == (3, 5, 1)
	run = run_star(genome_plain, reads, ATTRS + ["--scoreGapATAC", "-4"])
	assert len(run.alignments()) == 0
	assert ha.unmapped_counts(run.log_final()) == (0, len(reads), 0, 0)
	aln = run_star(genome_plain, reads, ATTRS + ["--scoreGapATAC", "-4", "--alignSJstitchMismatchNmax", "0", "-1", "0", "1"]).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		r = reads.by_name[rec.qname]
		left = sh.cigar_ops(r.cigar)[0][0]
		assert record_truth(rec) == (r.contig, r.pos, "%dM600N%dM" % (left + 3, 97 - left), r.strand)
		assert (rec.tag("nM"), rec.tag("AS")) == (1, 92)


def _homopolymer_insert_reads(sim, n=6, mismatch_before=(), prefix="hp"):
	"""Reads with one extra base inserted inside a run of 3 identical genome
	bases near the read centre; the run cannot be extended (the bases on both
	sides differ from it). 'run' = (first, last) genome position of the run.
	mismatch_before: distances, used in turn, at which the base that far left
	of the run is also substituted (by a base other than the run's), giving
	the alignment one mismatch ('mm' in extra)."""
	g = sim.g
	reads = []
	while len(reads) < n:
		c, p = sim.unique_locus(120)
		seq = g.contigs[c]
		for off in range(44, 56):
			q = p + off  # 1-based first base of the candidate run
			b = seq[q - 1]
			if seq[q:q + 2] == b * 2 and seq[q - 2] != b and seq[q + 2] != b:
				fwd = seq[p - 1:q + 2] + b + seq[q + 2:p - 1 + 99]
				mm = 0
				if mismatch_before:
					i = off - mismatch_before[len(reads) % len(mismatch_before)]
					sub = ha.mismatched(fwd[i])
					if sub == b:
						break
					fwd = fwd[:i] + sub + fwd[i + 1:]
					mm = 1
				reads.append(ha.oriented("%s%d" % (prefix, len(reads)), fwd, c, p, None, "+-"[len(reads) % 2],
					extra={"run": (q, q + 2), "mm": mm}))
				break
	return ReadSet(reads)


def test_insertion_flush_right(genome_plain, sim):
	"""--alignInsertionFlush Right moves an insertion whose position is
	ambiguous to its rightmost position (parametersDefault): an extra base
	in a 3-base run is reported after the run, i.e. the CIGAR is
	(run end - POS + 1)M1I...M, with AS 99 - 2 - 2 - 2 = 93; reads that also
	carry a substitution 2 to 6 bases left of the run score 98 - 1 - 4 - 2 =
	91 with nM 1. Without a substitution the default (None: "insertions are
	not flushed") also reports the insertion after the run; with one, it
	places at least one insertion further left (checked), so these reads
	tell the two settings apart."""
	reads = _homopolymer_insert_reads(sim) + _homopolymer_insert_reads(sim, 10, range(2, 7), prefix="hpmm")
	aln = run_star(genome_plain, reads, ATTRS + ["--alignInsertionFlush", "Right"]).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		r = reads.by_name[rec.qname]
		m1 = r.extra["run"][1] - r.pos + 1
		assert (rec.rname, rec.pos, rec.cigar, rec.strand) == (r.contig, r.pos, "%dM1I%dM" % (m1, 99 - m1), r.strand)
		assert (rec.tag("AS"), rec.tag("nM")) == (93 - 2 * r.extra["mm"], r.extra["mm"])
	right = {rec.qname: rec.cigar for rec in aln}
	default = {rec.qname: rec.cigar for rec in run_star(genome_plain, reads, ATTRS).alignments()}
	assert sorted(default) == sorted(reads.by_name)
	assert any(default[n] != right[n] for n in default if n.startswith("hpmm"))
