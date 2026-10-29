"""Area B: FLAG options (--outSAMflagOR/AND, --outSAMprimaryFlag, the Illumina
filter bit), unmapped reads in the alignment files (--outSAMunmapped None /
Within / Within KeepPairs) and in Unmapped.out.mate1/2 (--outReadsUnmapped
Fastx), --outSAMreadID, --outSAMorder, --outSAMtlen 2, --outSAMmode NoQS /
None and --outSAMtype None, in the SAM, unsorted-BAM and sorted-BAM writers."""

import pytest

import helpers_b as hb
import star_helpers as sh
from star_helpers import ReadSet, record_truth, run_star

FORMATS = list(sh.OUTPUT_FORMATS)
# lets a single mate (100 of 201 bases) pass the score and match filters
SINGLE_MATE_OK = ["--outFilterScoreMinOverLread", "0.3", "--outFilterMatchNminOverLread", "0.3"]


def _flag(strand, secondary=False):
	return (0x10 if strand == "-" else 0) | (0x100 if secondary else 0)


@pytest.mark.parametrize("fmt", FORMATS)
def test_flag_and_or(genome_plain, sim, fmt):
	"""--outSAMflagAND 65263 (clears 0x10 and 0x100) and --outSAMflagOR 1040
	(sets 0x400 and 0x10): parametersDefault says FLAG = (FLAG & AND) | OR,
	AND applied before OR. Every mapped record equals the record of a run
	without the options except for FLAG, which is (old & 65263) | 1040, so
	0x10 is set on every record (OR after AND) and 0x100 on none. The old
	FLAGs are checked against the truth (strand bit; one primary record per
	read)."""
	reads = ReadSet([sim.genomic_read(*sim.unique_locus(), strand=sim.rng.choice("+-")) for _ in range(20)]
		+ [sim.segment_read("dup3", sim.rng.randint(0, 400), strand=sim.rng.choice("+-")) for _ in range(5)])
	ref = run_star(genome_plain, reads, sh.OUTPUT_FORMATS[fmt]).alignments()
	new = run_star(genome_plain, reads, ["--outSAMflagAND", "65263", "--outSAMflagOR", "1040"]
		+ sh.OUTPUT_FORMATS[fmt]).alignments()
	for name, recs in ref.by_name().items():
		truth = reads.by_name[name]
		assert sum(not r.is_secondary for r in recs) == 1
		for r in recs:
			if truth.kind != "multi":
				assert r.flag == _flag(truth.strand)
			assert r.flag & ~0x110 == 0
	assert all(r.flag & 0x410 == 0x410 and not r.is_secondary for r in new)
	for r in ref:
		r.flag = (r.flag & 65263) | 1040
	assert sh.record_multiset(new) == sh.record_multiset(ref)


@pytest.mark.parametrize("fmt", FORMATS)
def test_primary_flag_best_score(genome_plain, sim, fmt):
	"""Reads from the 3-copy segment (three equal scores) and from the 2-copy
	segment dupmm whose copies differ at one base per 100, with
	--outFilterMultimapScoreRange 2 so the copy with one mismatch (AS 96,
	the exact copy has AS 98) is reported too. OneBestScore (default): one
	primary record per read, at a best-scoring locus; AllBestScore: every
	record with the best AS is primary and the AS-96 record is secondary
	(parametersDefault, manual 5.2.1). MAPQ follows NH (3 for 2, 1 for 3)."""
	reads = [sim.segment_read("dup3", sim.rng.randint(0, 500), strand=sim.rng.choice("+-")) for _ in range(5)]
	for i, off in enumerate((0, 110, 220, 330, 440, 560)):
		reads.append(sim.segment_read("dupmm", off, copy=i % 2, strand=sim.rng.choice("+-")))
	reads = ReadSet(reads)
	for mode in ("OneBestScore", "AllBestScore"):
		args = ["--outFilterMultimapScoreRange", "2", "--outSAMprimaryFlag", mode] + sh.OUTPUT_FORMATS[fmt]
		by_name = run_star(genome_plain, reads, args).alignments().by_name()
		assert sorted(by_name) == sorted(reads.by_name)
		for name, recs in by_name.items():
			truth = reads.by_name[name]
			n = truth.extra["n_loci"]
			assert len(recs) == n
			assert sorted((r.rname, r.pos, r.strand) for r in recs) == sorted(truth.extra["loci"])
			assert all(r.tag("NH") == n and r.mapq == {2: 3, 3: 1}[n] for r in recs)
			best = {(r.rname, r.pos) for r in recs if r.tag("AS") == 98}
			if truth.extra["segment"] == "dupmm":
				assert best == {(truth.contig, truth.pos)}
				assert sorted(r.tag("AS") for r in recs) == [96, 98]
			else:
				assert len(best) == 3
			primary = {(r.rname, r.pos) for r in recs if not r.is_secondary}
			if mode == "OneBestScore":
				assert len(primary) == 1 and primary <= best
			else:
				assert primary == best


@pytest.mark.parametrize("fmt", FORMATS)
def test_mult_nmax_writes_top_score_first(genome_plain, sim, fmt):
	"""--outSAMmultNmax N writes min(NH, N) records per read and, when N is
	not -1, "the top scoring alignment will be output first" (manual 5.2.1,
	parametersDefault); NH still counts all loci. dupmm reads (AS 98 and 96
	with --outFilterMultimapScoreRange 2) with N = 1: one record, the AS-98
	locus, primary, HI 1. 3-copy reads with N = 2: two records, HI 1 and 2,
	the first primary and the second secondary."""
	reads = ReadSet([sim.segment_read("dupmm", off, copy=i % 2, strand=sim.rng.choice("+-"))
		for i, off in enumerate((0, 110, 220, 330))]
		+ [sim.segment_read("dup3", sim.rng.randint(0, 500), strand=sim.rng.choice("+-")) for _ in range(4)])
	for nmax in (1, 2):
		args = ["--outFilterMultimapScoreRange", "2", "--outSAMmultNmax", str(nmax)] + sh.OUTPUT_FORMATS[fmt]
		aln = run_star(genome_plain, reads, args).alignments()
		if fmt != "sorted":
			assert hb.contiguous_by_name(aln.records)
		by_name = aln.by_name()
		assert sorted(by_name) == sorted(reads.by_name)
		for name, recs in by_name.items():
			truth = reads.by_name[name]
			n = truth.extra["n_loci"]
			assert len(recs) == min(n, nmax)
			assert all(r.tag("NH") == n for r in recs)
			assert {(r.rname, r.pos, r.strand) for r in recs} <= set(truth.extra["loci"])
			if fmt == "sorted":
				recs = sorted(recs, key=lambda r: r.tag("HI"))
			assert [r.tag("HI") for r in recs] == list(range(1, len(recs) + 1))
			assert [r.is_secondary for r in recs] == [False] + [True] * (len(recs) - 1)
			assert recs[0].tag("AS") == 98
			if truth.extra["segment"] == "dupmm" and nmax == 1:
				assert (recs[0].rname, recs[0].pos) == (truth.contig, truth.pos)


def _unmapped_mode_reads(sim):
	"""Proper pairs; pairs where one mate is unique and the other all N (each
	mate, each strand); pairs where one mate has 2 or 3 loci and the other is
	all N; pairs of two all-N mates."""
	pairs = list(sim.transcript_pairs(10))
	for i, (mate, strand) in enumerate(((1, "+"), (1, "-"), (2, "+"), (2, "-"))):
		name = "half%d" % i
		c, p = sim.unique_locus()
		pairs.append(hb.half_pair(name, sim.genomic_read(c, p, strand, name=name), mate))
	for i, (seg, mate) in enumerate((("dup3", 1), ("dup2", 2), ("dup3", 2))):
		name = "halfmulti%d" % i
		r = sim.segment_read(seg, sim.rng.randint(0, 400), strand=sim.rng.choice("+-"), name=name)
		pairs.append(hb.half_pair(name, r, mate))
	pairs += [hb.unmappable_pair("unmapped%d" % i) for i in range(2)]
	return ReadSet(pairs)


def _check_unmapped_mate(rec, mate, mapped_rec, input_read):
	"""SAM fields of the unmapped mate of a half-mapped pair, given the
	mapped record it refers to (SAM spec 1.4: 0x1, 0x4, 0x40/0x80, 0x20 if
	the mate is reversed; RNEXT/PNEXT the mate's position)."""
	assert rec.flag & ~0x100 == 0x1 | 0x4 | (0x40 if mate == 1 else 0x80) | (0x20 if mapped_rec.is_reverse else 0)
	assert (rec.rname, rec.pos, rec.mapq, rec.cigar, rec.tlen) == ("*", 0, 0, "*", 0)
	assert (rec.rnext, rec.pnext) == (mapped_rec.rname, mapped_rec.pos)
	assert (rec.seq, rec.qual) == (input_read.seq, input_read.qual)
	assert list(rec.tags) == ["NH", "HI", "AS", "nM", "uT"]
	assert (rec.tag("NH"), rec.tag("HI"), rec.tag("AS"), rec.tag("nM"), rec.tag("uT")) == (
		0, 0, mapped_rec.tag("AS"), mapped_rec.tag("nM"), "4")


@pytest.mark.parametrize("mode,fmt", [(m, f) for m in ("None", "Within") for f in FORMATS]
	+ [("Within KeepPairs", "sam"), ("Within KeepPairs", "bam")])
def test_unmapped_reads_in_alignment_files(genome_plain, sim, mode, fmt):
	"""--outSAMunmapped (parametersDefault, manual 5.4) with paired reads.
	None: no record has 0x4. Within: pairs with no mapped mate give two
	records (FLAG 77 and 141, '*' fields, uT 0); a half-mapped pair gives one
	record for its unmapped mate (uT 4, AS/nM of the mapped mate, RNEXT/PNEXT
	and 0x20 from the primary alignment). Within KeepPairs: one unmapped-mate
	record after each alignment of the mapped mate, 0x100 when that alignment
	is secondary; for unique reads this is the same single record as Within
	(it "only affects multi-mapping reads"). KeepPairs is not checked for the
	sorted BAM, which writes one unmapped record. In every
	mode the mapped mate of a half-mapped pair has 0x8, RNEXT '*', TLEN 0,
	and proper pairs are unchanged. In the sorted BAM the unmapped records
	come after all mapped ones."""
	reads = _unmapped_mode_reads(sim)
	args = SINGLE_MATE_OK + ["--outSAMunmapped"] + mode.split() + sh.OUTPUT_FORMATS[fmt]
	aln = run_star(genome_plain, reads, args).alignments()
	if fmt != "sorted":
		assert hb.contiguous_by_name(aln.records)
	else:
		flags = [r.is_unmapped for r in aln]
		assert flags == sorted(flags)
	by_name = aln.by_name()
	expected_names = [p.name for p in reads if mode != "None" or p.kind != "unmappable"]
	assert sorted(by_name) == sorted(expected_names)
	for name, recs in by_name.items():
		pair = reads.by_name[name]
		if pair.kind == "unmappable":
			assert [(r.flag, r.rname, r.pos, r.cigar, r.rnext, r.pnext, r.tlen) for r in recs] == [
				(77, "*", 0, "*", "*", 0, 0), (141, "*", 0, "*", "*", 0, 0)]
			for r in recs:
				assert (r.tag("NH"), r.tag("AS"), r.tag("uT")) == (0, 0, "0")
			continue
		if pair.kind != "half":
			assert len(recs) == 2 and not any(r.flag & 0xC for r in recs)
			assert all(record_truth(r) == hb.mate_truth(pair, r.mate).truth() for r in recs)
			continue
		mapped_mate = 1 if pair.r1.kind != "unmappable" else 2
		truth = hb.mate_truth(pair, mapped_mate)
		n = truth.extra.get("n_loci", 1)
		mapped = [r for r in recs if not r.is_unmapped]
		unmapped = [r for r in recs if r.is_unmapped]
		assert len(mapped) == n
		loci = truth.extra.get("loci", [(truth.contig, truth.pos, truth.strand)])
		assert sorted((r.rname, r.pos, r.strand) for r in mapped) == sorted(loci)
		for r in mapped:
			assert r.flag & ~0x110 == 0x1 | 0x8 | (0x40 if mapped_mate == 1 else 0x80)
			assert (r.rnext, r.pnext, r.tlen, r.cigar) == ("*", 0, 0, "100M")
		primary = [r for r in mapped if not r.is_secondary]
		assert len(primary) == 1
		other = hb.mate_truth(pair, 3 - mapped_mate)
		if mode == "None":
			assert not unmapped
		elif mode == "Within":
			assert len(unmapped) == 1
			_check_unmapped_mate(unmapped[0], 3 - mapped_mate, primary[0], other)
			assert not unmapped[0].is_secondary
			if fmt != "sorted":
				assert recs[-1] is unmapped[0]
		else:
			assert len(unmapped) == n and len(recs) == 2 * n
			for m, u in zip(recs[::2], recs[1::2]):
				assert not m.is_unmapped and u.is_unmapped
				_check_unmapped_mate(u, 3 - mapped_mate, m, other)
				assert u.is_secondary == m.is_secondary


def test_unmapped_fastx_pairs_and_names(genome_plain, sim):
	"""--outReadsUnmapped Fastx (manual 5.4): for pairs, both mates of every
	pair that is not mapped as a whole go to Unmapped.out.mate1/2 with the
	mate status 00 (neither mapped), 10 (mate 1 mapped) or 01 (mate 2
	mapped) at the end of the name line, after the mate index and filter
	flag "0:N:" / "1:N:"; the name is the first word without /1 /2.
	Sequence and quality are the input's. Proper pairs are not written. For
	FASTA single-end input the file is FASTA and the line has no status."""
	pairs = sim.transcript_pairs(6) + ReadSet([
		hb.half_pair("m1map", sim.genomic_read(*sim.unique_locus(), strand="+", name="m1map"), 1),
		hb.half_pair("m2map", sim.genomic_read(*sim.unique_locus(), strand="-", name="m2map"), 2),
		hb.unmappable_pair("nomap"),
	])
	run = run_star(genome_plain, pairs, SINGLE_MATE_OK + ["--outReadsUnmapped", "Fastx"], name_style="slash")
	status = {"m1map": "10", "m2map": "01", "nomap": "00"}
	assert [r.name for r in run.unmapped(1)] == [r.name for r in run.unmapped(2)]
	for mate in (1, 2):
		out = run.unmapped(mate)
		hb.assert_fastx_order([r.name for r in out], ["m1map", "m2map", "nomap"])
		for r in out:
			truth = hb.mate_truth(pairs.by_name[r.name], mate)
			assert r.description.split() == ["%d:N:" % (mate - 1), status[r.name]]
			assert (r.seq, r.qual) == (truth.seq, truth.qual)

	se = sim.transcript_reads(5) + ReadSet([sim.random_read() for _ in range(3)])
	run = run_star(genome_plain, se, ["--outReadsUnmapped", "Fastx"], fasta=True)
	text = run.path("Unmapped.out.mate1").read_text()
	assert text.startswith(">")
	out = run.unmapped(1)
	hb.assert_fastx_order([r.name for r in out], [r.name for r in se if r.kind == "random"])
	for r in out:
		assert r.description.split() == ["0:N:"] and r.qual is None and r.seq == se.by_name[r.name].seq


@pytest.mark.parametrize("fmt", FORMATS)
def test_illumina_filter_flag(genome_plain, sim, tmp_path, fmt):
	"""A read whose name line has 'Y' as the third character of its second
	field (Illumina "1:Y:0:..." = did not pass filter) gets FLAG 0x200
	(CHANGES.md: "Implemented 0x200 SAM flag for reads that did not pass
	Illumina filtering"), on mapped and on unmapped (--outSAMunmapped Within)
	records; 'N' reads do not. In Unmapped.out.mate1 the filter character
	follows the mate index: "0:Y:" / "0:N:"."""
	mapped = [sim.genomic_read(*sim.unique_locus(), strand=sim.rng.choice("+-")) for _ in range(8)]
	unmapped = [sim.random_read() for _ in range(4)]
	reads = mapped + unmapped
	not_passed = {r.name for i, r in enumerate(reads) if i % 2 == 0}
	path = hb.write_fastq(tmp_path / "reads.fq", [
		("%s 1:%s:0:ACGTAC" % (r.name, "Y" if r.name in not_passed else "N"), r.seq, r.qual) for r in reads])
	run = run_star(genome_plain, [path], ["--outSAMunmapped", "Within", "--outReadsUnmapped", "Fastx"]
		+ sh.OUTPUT_FORMATS[fmt])
	recs = run.records()
	assert sorted(r.qname for r in recs) == sorted(r.name for r in reads)
	truth = {r.name: r for r in reads}
	for rec in recs:
		t = truth[rec.qname]
		qc = 0x200 if rec.qname in not_passed else 0
		if t.kind == "random":
			assert rec.flag == 0x4 | qc
		else:
			assert record_truth(rec) == t.truth() and rec.flag == _flag(t.strand) | qc
	out = run.unmapped(1)
	hb.assert_fastx_order([r.name for r in out], [r.name for r in unmapped])
	for r in out:
		assert r.description.split() == ["0:%s:" % ("Y" if r.name in not_passed else "N")]


@pytest.mark.parametrize("fmt", FORMATS)
def test_read_id_number(genome_sjdb, sim, fmt):
	"""--outSAMreadID Number (parametersDefault: "read number (index) in the
	FASTx file"): QNAME is the index of the read in the input, both mates
	get it, as do unmapped records and the Unmapped.out.mate1/2 names. That
	the index is 1-based and runs on across comma-separated files pins the
	current output; parametersDefault does not say."""
	a = sim.transcript_pairs(6) + ReadSet([hb.unmappable_pair("unmA")])
	b = sim.transcript_pairs(5) + ReadSet([hb.unmappable_pair("unmB")])
	order = [p for p in a] + [p for p in b]
	args = ["--outSAMreadID", "Number", "--outSAMunmapped", "Within", "--outReadsUnmapped", "Fastx"]
	run = run_star(genome_sjdb, [a, b], args + sh.OUTPUT_FORMATS[fmt])
	by_name = run.alignments().by_name()
	assert sorted(by_name, key=int) == [str(i) for i in range(1, len(order) + 1)]
	for qname, recs in by_name.items():
		pair = order[int(qname) - 1]
		assert sorted(r.mate for r in recs) == [1, 2]
		for r in recs:
			if pair.kind == "unmappable":
				assert r.is_unmapped and r.seq == "N" * 100
			else:
				assert record_truth(r) == hb.mate_truth(pair, r.mate).truth()
	unmapped_ids = [str(i + 1) for i, p in enumerate(order) if p.kind == "unmappable"]
	assert [r.name for r in run.unmapped(1)] == [r.name for r in run.unmapped(2)]
	hb.assert_fastx_order([r.name for r in run.unmapped(1)], unmapped_ids)


def _mates_adjacent(records):
	"""Records come as consecutive mate pairs of one alignment (same QNAME
	and HI, mates 1 and 2), the leftmost mate first."""
	assert len(records) % 2 == 0
	for x, y in zip(records[::2], records[1::2]):
		assert x.qname == y.qname and x.tag("HI") == y.tag("HI")
		assert {x.mate, y.mate} == {1, 2}
		assert x.pos <= y.pos


def test_sam_order_keep_input_order(genome_sjdb, sim):
	"""--outSAMorder PairedKeepInputOrder keeps the input order (parametersDefault)
	with 4 threads mapping separate input chunks (--limitIObufferSize 220000
	gives chunks of about 40 pairs): the records equal those of a 1-thread
	run line for line. With the default Paired order the same records come
	out (as a multiset) and, in both orders, the two mates of each alignment
	are adjacent ("one mate after the other"), leftmost mate first (pins the
	current output). The test passes threads explicitly: the option only has
	its own code path with more than one thread."""
	pairs = sim.transcript_pairs(300)
	args = ["--limitIObufferSize", "220000", "50000000"]
	ref = run_star(genome_sjdb, pairs, args, threads=1).records()
	keep = run_star(genome_sjdb, pairs, args + ["--outSAMorder", "PairedKeepInputOrder"], threads=4).records()
	paired = run_star(genome_sjdb, pairs, args + ["--outSAMorder", "Paired"], threads=4).records()
	assert [r.to_line() for r in keep] == [r.to_line() for r in ref]
	assert [r.qname for r in keep[::2]] == [p.name for p in pairs]
	assert sh.record_multiset(paired) == sh.record_multiset(ref)
	_mates_adjacent(keep)
	_mates_adjacent(paired)


@pytest.mark.parametrize("fmt", ["bam", "sorted"])
def test_tlen_method2_for_non_protruding_pairs(genome_plain, sim, fmt):
	"""For mates that do not protrude past each other the two --outSAMtlen
	methods agree (parametersDefault: they differ only "for overlapping mates
	with protruding ends"): with --outSAMtlen 2 TLEN is the leftmost to the
	rightmost aligned base, + on the leftmost (the + strand) mate, for
	inserts of 250, 150 and 100 bases (the last: both mates cover the same
	bases). SAM output is not checked: its writer ignores --outSAMtlen."""
	pairs = []
	for insert in (250, 250, 150, 150, 100, 100):
		c, p = sim.unique_locus(260)
		pairs.append(sim.genomic_pair(c, p, insert, strand=sim.rng.choice("+-")))
	pairs = ReadSet(pairs)
	recs = run_star(genome_plain, pairs, ["--outSAMtlen", "2"] + sh.OUTPUT_FORMATS[fmt]).records()
	assert len(recs) == 2 * len(pairs)
	for rec in recs:
		pair = pairs.by_name[rec.qname]
		assert record_truth(rec) == hb.mate_truth(pair, rec.mate).truth()
		plus, minus = (pair.r1, pair.r2) if pair.r1.strand == "+" else (pair.r2, pair.r1)
		tlen = minus.end - plus.pos + 1
		assert tlen == pair.insert
		assert rec.tlen == (tlen if rec.strand == "+" else -tlen)


@pytest.mark.parametrize("fmt", FORMATS)
def test_concordant_protruding_pairs_are_proper(genome_plain, sim, fmt):
	"""--alignEndsProtrude 20 ConcordantPair "report[s] alignments with
	non-zero protrusion as concordant pairs" (parametersDefault): pairs whose
	+ mate starts 10 bases after the - mate (and ends 10 bases after it) get
	the proper-pair bit, FLAG 99/147 when mate 1 is the + mate and 83/163
	when mate 1 is the - mate, with RNEXT '=' and PNEXT the mate's POS."""
	pairs = []
	for i in range(6):
		c, p = sim.unique_locus(150)
		name = "prot%d" % i
		plus = sim.genomic_read(c, p + 10, "+", name=name)
		minus = sim.genomic_read(c, p, "-", name=name)
		r1, r2 = (plus, minus) if i % 2 == 0 else (minus, plus)
		r1.mate, r2.mate = 1, 2
		pairs.append(sh.SimPair(name, r1, r2))
	pairs = ReadSet(pairs)
	args = ["--alignEndsProtrude", "20", "ConcordantPair"] + sh.OUTPUT_FORMATS[fmt]
	recs = run_star(genome_plain, pairs, args).records()
	assert len(recs) == 12
	for rec in recs:
		pair = pairs.by_name[rec.qname]
		truth, other = (pair.r1, pair.r2) if rec.mate == 1 else (pair.r2, pair.r1)
		assert record_truth(rec) == truth.truth()
		assert rec.flag == ((99, 147) if pair.r1.strand == "+" else (83, 163))[rec.mate - 1]
		assert (rec.rnext, rec.pnext) == ("=", other.pos)


@pytest.mark.parametrize("fmt", FORMATS)
def test_sam_mode_noqs(genome_sjdb, sim, fmt):
	"""--outSAMmode NoQS: "full SAM but without quality scores"
	(parametersDefault): QUAL is '*' (BAM: 0xFF bytes) on every mapped record
	and all other fields equal those of a Full run. Unmapped records are
	left out: the SAM writer keeps their QUAL."""
	reads = sim.transcript_reads(30)
	full = run_star(genome_sjdb, reads, sh.OUTPUT_FORMATS[fmt]).records()
	noqs = run_star(genome_sjdb, reads, ["--outSAMmode", "NoQS"] + sh.OUTPUT_FORMATS[fmt]).records()
	assert all(r.qual == "*" for r in noqs)
	assert all(r.qual == reads.by_name[r.qname].sam_qual() for r in full)
	for r in full:
		r.qual = "*"
	assert sh.record_multiset(noqs) == sh.record_multiset(full)


@pytest.mark.parametrize("fmt", FORMATS)
def test_qs_conversion_add(genome_plain, sim, fmt):
	"""--outQSconversionAdd -15 "add[s] this number to the quality score"
	(parametersDefault): QUAL of every record, mapped or unmapped, and of
	Unmapped.out.mate1 is the input quality minus 15 (inputs are Phred 20-40,
	so no value reaches the ends of the printable range)."""
	reads = ReadSet([sim.genomic_read(*sim.unique_locus(), strand=sim.rng.choice("+-")) for _ in range(10)]
		+ [sim.random_read() for _ in range(3)])
	shift = lambda q: "".join(chr(ord(c) - 15) for c in q)
	args = ["--outQSconversionAdd", "-15", "--outSAMunmapped", "Within", "--outReadsUnmapped", "Fastx"]
	run = run_star(genome_plain, reads, args + sh.OUTPUT_FORMATS[fmt])
	recs = run.records()
	assert len(recs) == len(reads)
	for rec in recs:
		truth = reads.by_name[rec.qname]
		assert rec.qual == shift(truth.sam_qual() if not rec.is_unmapped else truth.qual)
	out = run.unmapped(1)
	assert len(out) == 3 and all(r.qual == shift(reads.by_name[r.name].qual) for r in out)


@pytest.mark.parametrize("option", [["--outSAMmode", "None"], ["--outSAMtype", "None"]])
def test_no_alignment_output(genome_sjdb, sim, option):
	"""--outSAMmode None ("no SAM output") and --outSAMtype None ("no SAM/BAM
	output") write no Aligned.* file; SJ.out.tab and the Log.final.out read
	counts equal those of a normal run."""
	reads = sim.transcript_reads(40, spliced=True) + sim.transcript_reads(20)
	ref = run_star(genome_sjdb, reads)
	run = run_star(genome_sjdb, reads, option)
	assert not [p.name for p in run.outdir.iterdir() if p.name.startswith("Aligned")]
	assert run.sj() == ref.sj() and len(ref.sj()) > 0
	for key in ("Number of input reads", "Uniquely mapped reads number", "Number of splices: Total"):
		assert run.log_final()[key] == ref.log_final()[key]
