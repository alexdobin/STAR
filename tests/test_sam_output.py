"""Area B: SAM/BAM output, attributes, flags, ordering and sorting, read
groups, unmapped-read output and --runThreadN equivalence."""

import pytest

import helpers_b as hb
import star_helpers as sh
from star_helpers import ReadSet, record_multiset, record_truth, run_star

ATTRS = "--outSAMattributes NH HI AS nM NM MD jM jI".split()


def _mixed_reads(sim):
	"""Unique spliced/unspliced reads, 2- and 3-copy multimappers, reads with
	substitutions, 12-copy reads (too many loci) and random reads."""
	reads = sim.transcript_reads(150) + sim.transcript_reads(40, spliced=True)
	extra = [sim.segment_read(s, sim.rng.randint(0, 200), strand=sim.rng.choice("+-"))
		for s in ("dup2", "dup3", "rep12") for _ in range(10)]
	extra += [sim.genomic_read(*sim.unique_locus(), strand="-", mismatches=(30, 64)) for _ in range(10)]
	extra += [sim.random_read() for _ in range(10)]
	return reads + ReadSet(extra)


def test_sam_unsorted_and_sorted_bam_agree(genome_sjdb, sim):
	"""The same reads written as SAM, unsorted BAM and coordinate-sorted BAM
	(--outSAMtype BAM Unsorted SortedByCoordinate) give the same multiset of
	records; the sorted BAM is ordered by (reference index, POS) with
	unmapped reads last and its header says SO:coordinate. The BAM side is
	decoded by star_helpers.read_bam."""
	reads = _mixed_reads(sim)
	args = ATTRS + ["--outSAMunmapped", "Within"]
	sam = run_star(genome_sjdb, reads, args).alignments()
	bam = run_star(genome_sjdb, reads, args + ["--outSAMtype", "BAM", "Unsorted", "SortedByCoordinate"])
	unsorted = bam.alignments("Aligned.out.bam")
	srt = bam.alignments("Aligned.sortedByCoord.out.bam")
	assert len(sam) > len(reads)  # multimappers add records
	assert record_multiset(unsorted) == record_multiset(sam)
	assert record_multiset(srt) == record_multiset(sam)
	ref_index = {name: i for i, (name, _) in enumerate(srt.references)}
	keys = [(ref_index[r.rname], r.pos) if not r.is_unmapped else (len(ref_index), 0) for r in srt]
	assert keys == sorted(keys)
	assert srt.header[0].split("\t")[:3] == ["@HD", "VN:1.4", "SO:coordinate"]
	assert srt.references == [(n, len(s)) for n, s in sim.g.contigs.items()]


def test_unmapped_within_reports_reason(genome_plain, sim):
	"""--outSAMunmapped Within writes one record per unmapped read: FLAG 4,
	RNAME '*', POS 0, CIGAR '*', SEQ/QUAL as in the input and uT giving the
	reason (manual 5.4): 3 for reads mapping to too many loci; the counts of
	uT 0, 1, 2, 3 equal Log.final.out's unmapped other / too short / too many
	mismatches / too many loci."""
	reads = sim.transcript_reads(40) + ReadSet(
		[sim.segment_read("rep12", sim.rng.randint(0, 200)) for _ in range(8)]
		+ [sim.random_read() for _ in range(12)])
	run = run_star(genome_plain, reads, "--outSAMunmapped Within")
	aln = run.alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	unmapped = [r for r in aln if r.is_unmapped]
	assert len(unmapped) == 20
	for rec in unmapped:
		truth = reads.by_name[rec.qname]
		assert (rec.flag, rec.rname, rec.pos, rec.cigar, rec.mapq) == (4, "*", 0, "*", 0)
		assert (rec.seq, rec.qual) == (truth.seq, truth.qual)
		assert rec.tag("uT") == "3" if truth.kind == "multi" else rec.tag("uT") in ("0", "1", "2")
	log = run.log_final()
	counts = {k: sum(r.tag("uT") == k for r in unmapped) for k in "0123"}
	assert counts["3"] == 8
	assert counts["0"] == log["Number of reads unmapped: other"]
	assert counts["1"] == log["Number of reads unmapped: too short"]
	assert counts["2"] == log["Number of reads unmapped: too many mismatches"]
	assert counts["3"] == log["Number of reads mapped to too many loci"]


def test_unmapped_fastx_output(genome_plain, sim):
	"""--outReadsUnmapped Fastx writes the unmapped reads, and only those, to
	Unmapped.out.mate1/2 with the input sequence and quality. For pairs where
	neither mate maps the name carries the mate-status tag 00 (manual 5.4)."""
	se = sim.transcript_reads(30) + ReadSet([sim.random_read() for _ in range(10)])
	run = run_star(genome_plain, se, "--outReadsUnmapped Fastx")
	out = run.unmapped(1)
	assert sorted(r.name for r in out) == sorted(r.name for r in se if r.kind == "random")
	for r in out:
		truth = se.by_name[r.name]
		assert (r.seq, r.qual) == (truth.seq, truth.qual)

	rnd = [sh.SimPair(n, sim.random_read(name=n), sim.random_read(name=n)) for n in ("rndA", "rndB", "rndC")]
	pe = sim.transcript_pairs(20) + ReadSet(rnd)
	run = run_star(genome_plain, pe, "--outReadsUnmapped Fastx")
	assert [r.name for r in run.unmapped(1)] == [r.name for r in run.unmapped(2)]
	for mate in (1, 2):
		out = run.unmapped(mate)
		hb.assert_fastx_order([r.name for r in out], ["rndA", "rndB", "rndC"])
		for r in out:
			truth = getattr(pe.by_name[r.name], "r%d" % mate)
			assert (r.seq, r.qual) == (truth.seq, truth.qual)
			assert r.description.split()[-1] == "00"


def test_attribute_order_and_junction_tags(genome_sjdb, sim):
	"""Attributes come out in the order requested. jM lists the motif of each
	junction plus 20 for annotated ones, jI the 1-based intron start/end
	pairs; unspliced reads get jM -1 and jI -1. NM/MD match expected_nm_md."""
	reads = sim.transcript_reads(80, spliced=True) + sim.transcript_reads(20, spliced=False)
	order = ["MD", "jI", "NH", "AS", "jM", "NM", "HI", "nM"]
	aln = run_star(genome_sjdb, reads, ["--outSAMattributes"] + order).alignments()
	annotated = {(j.contig, j.start, j.end): j.motif for j in sim.g.junctions()}
	assert len(aln) == len(reads)
	for rec in aln:
		truth = reads.by_name[rec.qname]
		assert record_truth(rec) == truth.truth()
		assert list(rec.tags) == order
		assert rec.tag_types["jM"] == "B:c" and rec.tag_types["jI"] == "B:i"
		if truth.introns:
			assert rec.tag("jM") == [annotated[(truth.contig, s, e)] + 20 for s, e in truth.introns]
			assert rec.tag("jI") == [x for se in truth.introns for x in se]
		else:
			assert rec.tag("jM") == [-1] and rec.tag("jI") == [-1]
		assert (rec.tag("NM"), rec.tag("MD")) == sh.expected_nm_md(
			sim.g, truth.contig, truth.pos, truth.cigar, truth.sam_seq())


def _dup3_reads(sim, n=12):
	return ReadSet([sim.segment_read("dup3", sim.rng.randint(0, 500), strand=sim.rng.choice("+-"))
		for _ in range(n)])


def test_primary_flag_all_best_score(genome_plain, sim):
	"""The three copies of dup3 score the same, so with --outSAMprimaryFlag
	AllBestScore all three records of each read are primary (no 0x100), while
	the default OneBestScore leaves exactly one primary."""
	reads = _dup3_reads(sim)
	one = run_star(genome_plain, reads).alignments().by_name()
	allbest = run_star(genome_plain, reads, "--outSAMprimaryFlag AllBestScore").alignments().by_name()
	assert sorted(one) == sorted(allbest) == sorted(reads.by_name)
	for name in one:
		assert len(one[name]) == len(allbest[name]) == 3
		assert sum(not r.is_secondary for r in one[name]) == 1
		assert all(not r.is_secondary for r in allbest[name])


def test_mult_nmax_limits_records_not_nh(genome_plain, sim):
	"""--outSAMmultNmax 1 writes exactly one line per multimapper; NH still
	reports the number of loci (manual 5.2.1), and the line is primary and at
	one of the true loci."""
	reads = _dup3_reads(sim)
	aln = run_star(genome_plain, reads, "--outSAMmultNmax 1").alignments()
	by_name = aln.by_name()
	assert sorted(by_name) == sorted(reads.by_name)
	for name, recs in by_name.items():
		assert len(recs) == 1
		rec = recs[0]
		assert rec.tag("NH") == 3 and not rec.is_secondary
		assert (rec.rname, rec.pos, rec.strand) in reads.by_name[name].extra["loci"]


def test_read_groups_from_comma_separated_inputs(genome_sjdb, sim):
	"""Two read files given comma-separated with --outSAMattrRGline ID:grpA
	SM:a , ID:grpB SM:b: the header has both @RG lines and every record
	carries RG:Z of the file its read came from."""
	a = sim.transcript_reads(30)
	b = sim.transcript_reads(30)
	aln = run_star(genome_sjdb, [a, b], "--outSAMattrRGline ID:grpA SM:a , ID:grpB SM:b").alignments()
	assert aln.header_lines("@RG") == ["@RG\tID:grpA\tSM:a", "@RG\tID:grpB\tSM:b"]
	assert len(aln) == 60
	for rec in aln:
		assert rec.tag("RG") == ("grpA" if rec.qname in a.by_name else "grpB")


def test_strand_field_intron_motif(genome_sjdb, sim):
	"""--outSAMstrandField intronMotif adds XS:A with the strand of the intron
	motif to spliced alignments (+ for GT/AG, GC/AG, AT/AC; - for their
	reverse complements) whatever the read's own strand, and no XS to
	unspliced alignments."""
	reads = sim.transcript_reads(60, spliced=True) + sim.transcript_reads(20, spliced=False)
	aln = run_star(genome_sjdb, reads, "--outSAMstrandField intronMotif").alignments()
	assert len(aln) == len(reads)
	for rec in aln:
		truth = reads.by_name[rec.qname]
		assert record_truth(rec) == truth.truth()
		if truth.introns:
			tx = sim.g.transcripts[truth.tx]
			assert rec.tag("XS") == tx.strand
		else:
			assert "XS" not in rec.tags


def _comparable_log(log):
	skip = ("Started job on", "Started mapping on", "Finished on", "Mapping speed, Million of reads per hour")
	return {k: v for k, v in log.items() if k not in skip}


def test_thread_count_does_not_change_output(genome_sjdb, sim):
	"""--runThreadN 1 and 4 give the same multiset of SAM records (unmapped
	reads included), the same SJ.out.tab and the same Log.final.out
	statistics, for single-end and paired-end input. --limitIObufferSize
	106000 per mate leaves input chunks of about 20 reads, so the four
	threads map different reads; with the default buffer one thread would
	load the whole input as a single chunk."""
	se = _mixed_reads(sim)
	pe = sim.transcript_pairs(120)
	args = ATTRS + ["--outSAMunmapped", "Within"]
	for reads in (se, pe):
		chunks = ["--limitIObufferSize", str(106000 * (2 if reads.paired else 1)), "50000000"]
		one = run_star(genome_sjdb, reads, args + chunks, threads=1)
		four = run_star(genome_sjdb, reads, args + chunks, threads=4)
		a, b = record_multiset(one.records()), record_multiset(four.records())
		assert a == b, sh.first_line_diff(a, b)
		assert one.sj() == four.sj()
		assert _comparable_log(one.log_final()) == _comparable_log(four.log_final())


def _protruding_pairs(sim, n=4):
	"""Pairs whose + mate (P+10..P+109) starts 10 bases after the - mate
	(P..P+99), so the two TLEN methods differ."""
	pairs = []
	for i in range(n):
		c, p = sim.unique_locus(150)
		name = "prot%d" % i
		r1 = sim.genomic_read(c, p + 10, "+", name=name)
		r2 = sim.genomic_read(c, p, "-", name=name)
		r1.mate, r2.mate = 1, 2
		pairs.append(sh.SimPair(name, r1, r2))
	return ReadSet(pairs)


@pytest.mark.parametrize("fmt", list(sh.OUTPUT_FORMATS))
def test_tlen_methods_for_protruding_mates(genome_plain, sim, fmt):
	"""Mates protruding by 10 bases (allowed by --alignEndsProtrude 20
	ConcordantPair). outSAMtlen 1 (default): + mate start to - mate end = 90,
	+ on the + mate. outSAMtlen 2: leftmost to rightmost base of either mate =
	110, + on the mate with the leftmost base (the - mate). Method 2 is
	checked for BAM output only: the SAM writer
	(ReadAlign_outputTranscriptSAM.cpp) does not read --outSAMtlen."""
	pairs = _protruding_pairs(sim)
	base = ["--alignEndsProtrude", "20", "ConcordantPair"] + sh.OUTPUT_FORMATS[fmt]
	recs = run_star(genome_plain, pairs, base).records()
	assert len(recs) == 8
	for rec in recs:
		truth = getattr(pairs.by_name[rec.qname], "r%d" % rec.mate)
		assert record_truth(rec) == truth.truth()
		assert rec.tlen == (-90 if rec.is_reverse else 90)
	if fmt != "sam":
		recs = run_star(genome_plain, pairs, base + ["--outSAMtlen", "2"]).records()
		assert len(recs) == 8
		for rec in recs:
			assert rec.tlen == (110 if rec.is_reverse else -110)
