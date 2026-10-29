"""Area A: read input, seed search, stitching and extension, scoring, output
filters and multimapper handling. Expected values come from the simulation
truth (star_helpers.SimRead) or are worked by hand in the docstrings."""

import math

import pytest

import star_helpers as sh
from star_helpers import ReadSet, record_truth, run_star


def _assert_truth(aln, reads):
	"""Every input read has exactly one record, at its simulated alignment."""
	by_name = aln.by_name()
	assert sorted(by_name) == sorted(reads.by_name)
	for name, recs in by_name.items():
		truth = reads.by_name[name]
		assert len(recs) == 1, recs
		rec = recs[0]
		assert record_truth(rec) == truth.truth(), (rec, truth)
		assert rec.seq == truth.sam_seq()
		assert rec.qual == truth.sam_qual()


def test_se_reads_align_to_their_origin(genome_sjdb, sim):
	"""Unspliced and spliced single-end reads from annotated transcripts align
	uniquely (NH 1, MAPQ 255) with the simulated position, strand and CIGAR;
	SEQ and QUAL are reverse-complemented for reverse-strand reads."""
	reads = sim.transcript_reads(240) + sim.transcript_reads(60, spliced=True)
	assert sum("N" in r.cigar for r in reads) >= 60
	assert sum(r.strand == "-" for r in reads) >= 100
	aln = run_star(genome_sjdb, reads).alignments()
	_assert_truth(aln, reads)
	for rec in aln:
		assert rec.mapq == 255
		assert rec.tag("NH") == 1 and rec.tag("HI") == 1
		assert rec.flag in (0, 16)


def test_spliced_reads_found_without_annotation(genome_plain, sim):
	"""Without a splice-junction database, spliced reads with at least 25
	bases on each side of every junction still get the true CIGAR, including
	GC/AG, AT/AC and minus-strand introns and the 5 kb intron of G09."""
	reads = sim.transcript_reads(200, spliced=True)
	motifs = {sh.motif_code(sim.g.seq(r.contig, s, s + 1), sim.g.seq(r.contig, e - 1, e))
		for r in reads for s, e in r.introns}
	assert motifs >= {1, 2, 3, 4, 5, 6}
	aln = run_star(genome_plain, reads).alignments()
	_assert_truth(aln, reads)


@pytest.mark.parametrize("fmt", list(sh.OUTPUT_FORMATS))
def test_substitutions_and_indels(genome_plain, sim, fmt):
	"""Reads with two substitutions, a 3-base deletion or a 2-base insertion
	(placed where they cannot slide) align with the simulated CIGAR; nM, NM
	and MD equal the values computed from the truth (expected_nm_md), in SAM,
	unsorted BAM and sorted BAM output (separate code paths in STAR)."""
	reads = []
	while len(reads) < 60:
		contig, pos = sim.unique_locus(120)
		kind = len(reads) % 3
		strand = "+" if len(reads) % 2 else "-"
		if kind == 0:
			r = sim.genomic_read(contig, pos, strand, mismatches=(25, 71))
		elif kind == 1:
			if not sim.deletion_is_unambiguous(contig, pos + 45, 3):
				continue
			r = sim.genomic_read(contig, pos, strand, deletion=(45, 3))
		else:
			r = sim.genomic_read(contig, pos, strand, insertion=(52, 2))
		reads.append(r)
	reads = ReadSet(reads)
	args = ["--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD"] + sh.OUTPUT_FORMATS[fmt]
	aln = run_star(genome_plain, reads, args).alignments()
	_assert_truth(aln, reads)
	for rec in aln:
		truth = reads.by_name[rec.qname]
		nm, md = sh.expected_nm_md(sim.g, truth.contig, truth.pos, truth.cigar, truth.sam_seq())
		assert rec.tag("NM") == nm, rec
		assert rec.tag("MD") == md, rec
		assert rec.tag("nM") == len(truth.mismatches)


def test_alignment_score_hand_worked(genome_plain, sim):
	"""AS = matches - mismatches + gap penalties + the genomic-length term
	int(ceil(log2(span) * -0.25 - 0.5)) from scoreGenomicLengthLog2scale.
	100M exact: 100 - 2 = 98. Two mismatches: 98 - 4 = 94. 3-base deletion
	(scoreDelOpen -2, scoreDelBase -2 x 3): 100 - 8 - 2 = 90 (span 103 still
	gives -2). 2-base insertion (scoreInsOpen -2, scoreInsBase -2 x 2): 98
	matched bases - 6 - 2 = 90."""
	contig, pos = sim.unique_locus(120)
	while not sim.deletion_is_unambiguous(contig, pos + 45, 3):
		contig, pos = sim.unique_locus(120)
	reads = ReadSet([
		sim.genomic_read(contig, pos, "+", name="exact"),
		sim.genomic_read(contig, pos, "-", mismatches=(30, 60), name="mm2"),
		sim.genomic_read(contig, pos, "+", deletion=(45, 3), name="del3"),
		sim.genomic_read(contig, pos, "-", insertion=(50, 2), name="ins2"),
	])
	aln = run_star(genome_plain, reads).alignments()
	_assert_truth(aln, reads)
	got = {r.qname: r.tag("AS") for r in aln}
	assert got == {"exact": 98, "mm2": 94, "del3": 90, "ins2": 90}


def _mapq(n):
	return 255 if n == 1 else int(-10 * math.log10(1 - 1.0 / n))


def test_multimappers_report_all_loci(genome_plain, sim):
	"""Reads from the 2-copy (one copy reverse-complemented) and 3-copy
	segments get one record per copy: NH = copies, HI = 1..NH, MAPQ =
	int(-10*log10(1-1/NH)) (3 for NH 2, 1 for NH 3; STAR manual 5.2.1), the
	loci and strands of the copies, and exactly one record without 0x100."""
	reads = []
	for seg, n in (("dup2", 15), ("dup3", 15)):
		for i in range(n):
			off = sim.rng.randint(0, sim.g.segments[seg].length - 100)
			reads.append(sim.segment_read(seg, off, strand=sim.rng.choice("+-")))
	reads = ReadSet(reads)
	aln = run_star(genome_plain, reads).alignments()
	by_name = aln.by_name()
	assert sorted(by_name) == sorted(reads.by_name)
	for name, recs in by_name.items():
		truth = reads.by_name[name]
		n = truth.extra["n_loci"]
		assert len(recs) == n
		assert sorted((r.rname, r.pos, r.strand) for r in recs) == sorted(truth.extra["loci"])
		assert sorted(r.tag("HI") for r in recs) == list(range(1, n + 1))
		assert all(r.tag("NH") == n and r.mapq == _mapq(n) and r.cigar == "100M" for r in recs)
		assert sum(not r.is_secondary for r in recs) == 1


def test_too_many_loci(genome_plain, sim):
	"""Reads from the 12-copy segment exceed outFilterMultimapNmax 10: no
	records, and Log.final.out counts them as mapped to too many loci. With
	--outFilterMultimapNmax 12 each gets 12 records with NH 12 and MAPQ 0."""
	reads = ReadSet([sim.segment_read("rep12", sim.rng.randint(0, 200)) for _ in range(12)])
	run = run_star(genome_plain, reads)
	assert len(run.alignments()) == 0
	log = run.log_final()
	assert log["Number of reads mapped to too many loci"] == 12
	assert log["Uniquely mapped reads number"] == 0
	assert log["Number of reads mapped to multiple loci"] == 0

	aln = run_star(genome_plain, reads, "--outFilterMultimapNmax 12").alignments()
	by_name = aln.by_name()
	assert len(by_name) == 12
	for name, recs in by_name.items():
		assert len(recs) == 12
		assert sorted((r.rname, r.pos, r.strand) for r in recs) == sorted(reads.by_name[name].extra["loci"])
		assert all(r.tag("NH") == 12 and r.mapq == 0 for r in recs)


def test_multimap_score_range(genome_plain, sim):
	"""Copy 2 of segment dupmm differs from copy 1 by one substitution every
	100 bases, so a 100-base read from copy 1 scores 98 there and 96 on copy
	2. With the default outFilterMultimapScoreRange 1 the read is unique at
	copy 1; with range 2 both loci are reported (NH 2), copy 2 with nM 1."""
	reads = []
	for i in range(10):
		off = 100 * sim.rng.randint(0, 8) + sim.rng.randint(5, 45)  # one substitution at off..off+99
		reads.append(sim.segment_read("dupmm", off, strand=sim.rng.choice("+-")))
	reads = ReadSet(reads)
	aln = run_star(genome_plain, reads).alignments()
	_assert_truth(aln, reads)
	assert all(r.tag("NH") == 1 and r.tag("AS") == 98 for r in aln)

	aln = run_star(genome_plain, reads, "--outFilterMultimapScoreRange 2").alignments()
	for name, recs in aln.by_name().items():
		truth = reads.by_name[name]
		assert sorted((r.rname, r.pos, r.strand) for r in recs) == sorted(truth.extra["loci"])
		assert all(r.tag("NH") == 2 for r in recs)
		scores = {(r.rname, r.pos): (r.tag("AS"), r.tag("nM")) for r in recs}
		assert scores[(truth.contig, truth.pos)] == (98, 0)
		other = [v for k, v in scores.items() if k != (truth.contig, truth.pos)]
		assert other == [(96, 1)]


@pytest.mark.parametrize("fasta,gz,name_style", [
	(False, True, "plain"),
	(True, False, "plain"),
	(False, False, "slash"),
	(False, False, "illumina"),
])
def test_read_input_formats(genome_sjdb, sim, fasta, gz, name_style):
	"""Gzipped FASTQ through --readFilesCommand, FASTA, and read names with
	/1 or an Illumina comment after a space all give the true alignments; the
	QNAME is the name without '/1' or the comment (readNameSeparator '/'), and
	FASTA input gives QUAL '*'."""
	reads = sim.transcript_reads(60)
	aln = run_star(genome_sjdb, reads, fasta=fasta, gz=gz, name_style=name_style).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		truth = reads.by_name[rec.qname]
		assert record_truth(rec) == truth.truth()
		assert rec.seq == truth.sam_seq()
		assert rec.qual == ("*" if fasta else truth.sam_qual())


def test_clip_fixed_bases(genome_plain, sim):
	"""--clip5pNbases 4 --clip3pNbases 6 soft-clips 4 bases at the read's 5'
	end and 6 at its 3' end: a + strand read at P becomes 4S90M6S at P+4, a
	- strand read at P becomes 6S90M4S at P+6."""
	reads = ReadSet([sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2]) for i in range(20)])
	aln = run_star(genome_plain, reads, "--clip5pNbases 4 --clip3pNbases 6").alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		truth = reads.by_name[rec.qname]
		if truth.strand == "+":
			assert (rec.pos, rec.cigar, rec.strand) == (truth.pos + 4, "4S90M6S", "+")
		else:
			assert (rec.pos, rec.cigar, rec.strand) == (truth.pos + 6, "6S90M4S", "-")


def test_clip_3p_adapter_and_polya(genome_plain, sim):
	"""A 60-base insert followed by 40 adapter bases aligns as 60M40S (+) or
	40S60M (-) at the insert's position with --clip3pAdapterSeq; the 60
	aligned bases pass the 0.66 filters because they are normalised to the
	clipped length. --clip3pAdapterSeq polyA clips a 30-base poly-A tail the
	same way (reads chosen so the genomic part does not end in A)."""
	adapter = sh.ILLUMINA_ADAPTER
	reads = []
	for i in range(16):
		contig, pos = sim.unique_locus()
		reads.append(sim.genomic_read(contig, pos, "+-"[i % 2], length=60, tail=adapter + adapter[:6]))
	reads = ReadSet(reads)
	aln = run_star(genome_plain, reads, ["--clip3pAdapterSeq", adapter]).alignments()
	_assert_truth(aln, reads)

	polya = []
	while len(polya) < 16:
		contig, pos = sim.unique_locus()
		r = sim.genomic_read(contig, pos, "+-"[len(polya) % 2], length=70, tail="A" * 30)
		if "A" not in r.seq[64:70]:
			polya.append(r)
	polya = ReadSet(polya)
	aln = run_star(genome_plain, polya, "--clip3pAdapterSeq polyA").alignments()
	_assert_truth(aln, polya)


@pytest.mark.parametrize("fmt", list(sh.OUTPUT_FORMATS))
def test_paired_end_reads(genome_sjdb, sim, fmt):
	"""Pairs from transcripts align as proper pairs: FLAG 99/147 or 83/163,
	RNEXT '=', PNEXT the mate's POS, TLEN +-(end of the - mate - start of the
	+ mate + 1) with + on the + mate (outSAMtlen 1), both mates at their true
	alignment, NH 1 and MC the mate's CIGAR; in SAM and both BAM outputs."""
	pairs = sim.transcript_pairs(150)
	args = ["--outSAMattributes", "NH", "HI", "AS", "nM", "MC"] + sh.OUTPUT_FORMATS[fmt]
	aln = run_star(genome_sjdb, pairs, args).alignments()
	by_name = aln.by_name()
	assert sorted(by_name) == sorted(pairs.by_name)
	for name, recs in by_name.items():
		p = pairs.by_name[name]
		assert sorted(r.mate for r in recs) == [1, 2]
		m = {r.mate: r for r in recs}
		plus, minus = (p.r1, p.r2) if p.r1.strand == "+" else (p.r2, p.r1)
		tlen = minus.end - plus.pos + 1
		for i, truth, other in ((1, p.r1, p.r2), (2, p.r2, p.r1)):
			rec = m[i]
			assert record_truth(rec) == truth.truth()
			assert rec.flag in ((99, 147) if p.r1.strand == "+" else (83, 163))[i - 1:i]
			assert (rec.rnext, rec.pnext) == ("=", other.pos)
			assert rec.tlen == (tlen if truth is plus else -tlen)
			assert rec.tag("NH") == 1 and rec.tag("MC") == other.cigar
