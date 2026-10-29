"""Area B: SAM attribute values (NH HI AS nM NM MD jM jI MC ch XS RG), the
--outSAMattributes presets, --outSAMattrIHstart, the MAPQ ladder with
--outSAMmapqUnique, --outSAMstrandField intronMotif and --outSAMattrRGline,
each checked in the SAM, unsorted-BAM and sorted-BAM writers, which are
separate code in STAR (ReadAlign_outputTranscriptSAM.cpp and
ReadAlign_alignBAM.cpp)."""

import pytest

import helpers_b as hb
import star_helpers as sh
from star_helpers import ReadSet, record_truth, run_star

FORMATS = list(sh.OUTPUT_FORMATS)
BAM_FORMATS = ["bam", "sorted"]
ATTRS = ["NH", "HI", "AS", "nM", "NM", "MD", "jM", "jI", "MC"]
INT_TAGS = ("NH", "HI", "AS", "nM", "NM")


@pytest.fixture(scope="session")
def tiny_genome(genome_factory, tmp_path_factory):
	"""STAR index of helpers_b's tiny genome (no annotation)."""
	return genome_factory(gtf=False, synth=hb.tiny_synth(tmp_path_factory))


def _annotated(synth):
	return {(j.contig, j.start, j.end) for j in synth.junctions()}


def _se_reads(sim):
	"""Spliced and unspliced transcript reads, reads with two substitutions,
	a 3-base deletion or a 2-base insertion at unique loci, and reads from
	the 2- and 3-copy segments."""
	reads = list(sim.transcript_reads(50)) + list(sim.transcript_reads(20, spliced=True))
	for i in range(8):
		c, p = sim.unique_locus(120)
		reads.append(sim.genomic_read(c, p, sim.rng.choice("+-"), mismatches=(30, 64)))
	n = 0
	while n < 8:
		c, p = sim.unique_locus(120)
		if sim.deletion_is_unambiguous(c, p + 50, 3):
			reads.append(sim.genomic_read(c, p, sim.rng.choice("+-"), deletion=(50, 3)))
			n += 1
	for i in range(8):
		c, p = sim.unique_locus(120)
		reads.append(sim.genomic_read(c, p, sim.rng.choice("+-"), insertion=(50, 2)))
	for seg in ("dup2", "dup3"):
		for i in range(5):
			reads.append(sim.segment_read(seg, sim.rng.randint(0, 400), strand=sim.rng.choice("+-")))
	return ReadSet(reads)


def _pe_reads(sim):
	"""Transcript pairs (some spliced) and genomic pairs with substitutions
	in both mates."""
	pairs = list(sim.transcript_pairs(40))
	for i in range(10):
		c, p = sim.unique_locus(260)
		pairs.append(sim.genomic_pair(c, p, 250, strand=sim.rng.choice("+-"), mismatches1=(30,),
			mismatches2=(25, 70)))
	return ReadSet(pairs)


def _check_tags(rec, truth, mates, synth, annotated, expected_order):
	"""Every attribute of one mapped record against the truth of its mate
	(`truth`) and of its whole read (`mates`)."""
	assert list(rec.tags) == expected_order
	for t in INT_TAGS:
		assert rec.tag_types[t] == "i"
	assert rec.tag_types["MD"] == "Z" and rec.tag_types["jM"] == "B:c" and rec.tag_types["jI"] == "B:i"
	assert rec.tag("AS") == hb.expected_score(synth, mates, annotated)
	assert rec.tag("nM") == sum(len(m.mismatches) for m in mates)
	seq = truth.seq if rec.strand == "+" else sh.revcomp(truth.seq)
	assert rec.seq == seq
	assert (rec.tag("NM"), rec.tag("MD")) == sh.expected_nm_md(synth, rec.rname, rec.pos, truth.cigar, seq)
	introns = sh.cigar_introns(rec.pos, rec.cigar)
	if introns:
		motifs = [sh.motif_code(synth.seq(rec.rname, s, s + 1), synth.seq(rec.rname, e - 1, e))
			+ (20 if (rec.rname, s, e) in annotated else 0) for s, e in introns]
		assert rec.tag("jM") == motifs
		assert rec.tag("jI") == [x for se in introns for x in se]
	else:
		assert rec.tag("jM") == [-1] and rec.tag("jI") == [-1]


@pytest.mark.parametrize("fmt", FORMATS)
def test_attribute_values(genome_sjdb, sim, fmt):
	"""--outSAMattributes NH HI AS nM NM MD jM jI MC (+ ch for BAM) on
	single-end reads with substitutions, indels, annotated junctions and
	2-3 loci, and on pairs. Each value against the truth: NH = number of
	loci, HI = 1..NH, AS = helpers_b.expected_score (manual 5.2.2: +1/-1 for
	matches/mismatches plus gap penalties; for PE the total of both mates), nM
	= substitutions in the read or pair, NM/MD per mate from the SAM
	specification (expected_nm_md), jM = motif + 20 for annotated junctions,
	jI = 1-based intron start/end (-1 for unspliced), MC = the other mate's
	CIGAR (pairs only), ch absent (it marks chimeric alignments only). Tags
	come in the requested order and with the specification's types."""
	g = sim.g
	annotated = _annotated(g)
	attrs = ATTRS + (["ch"] if fmt != "sam" else [])
	args = ["--outSAMattributes"] + attrs + sh.OUTPUT_FORMATS[fmt]

	se = _se_reads(sim)
	by_name = run_star(genome_sjdb, se, args).alignments().by_name()
	assert sorted(by_name) == sorted(se.by_name)
	order = [a for a in ATTRS if a != "MC"]
	for name, recs in by_name.items():
		truth = se.by_name[name]
		n = truth.extra.get("n_loci", 1)
		assert len(recs) == n
		assert sorted(r.tag("HI") for r in recs) == list(range(1, n + 1))
		if n == 1:
			assert record_truth(recs[0]) == truth.truth()
		else:
			assert sorted((r.rname, r.pos, r.strand) for r in recs) == sorted(truth.extra["loci"])
		for rec in recs:
			assert rec.tag("NH") == n
			_check_tags(rec, truth, [truth], g, annotated, order)

	pe = _pe_reads(sim)
	by_name = run_star(genome_sjdb, pe, args).alignments().by_name()
	assert sorted(by_name) == sorted(pe.by_name)
	for name, recs in by_name.items():
		pair = pe.by_name[name]
		assert sorted(r.mate for r in recs) == [1, 2]
		for rec in recs:
			truth, other = (pair.r1, pair.r2) if rec.mate == 1 else (pair.r2, pair.r1)
			assert record_truth(rec) == truth.truth()
			assert (rec.tag("NH"), rec.tag("HI")) == (1, 1)
			_check_tags(rec, truth, [pair.r1, pair.r2], g, annotated, ATTRS)
			assert rec.tag("MC") == other.cigar


@pytest.mark.parametrize("fmt", BAM_FORMATS)
def test_bam_integer_tag_types(genome_sjdb, sim, fmt):
	"""In BAM, single integer tags may use any of the types cCsSiI, "chosen
	according to the field value's magnitude" (SAM specification 4.2.4);
	STAR writes the smallest unsigned type that holds the value (pins the
	current output: star_helpers.read_bam reports every integer type as
	'i', so only the raw bytes show it). Single-end reads give values below
	256 (type C); pairs of 150-base mates give AS near 298 (type S). jM and
	jI are B arrays of c and i (manual 5.2.2) and MD is Z."""
	pe = sim.transcript_pairs(20, length=150, insert=(320, 400))
	for reads in (_se_reads(sim), pe):
		run = run_star(genome_sjdb, reads, ["--outSAMattributes"] + ATTRS[:-1] + sh.OUTPUT_FORMATS[fmt])
		raw = run.alignments().records
		types = hb.bam_aux_types(run.main_alignment_path().read_bytes())
		assert len(types) == len(raw)
		seen = set()
		for rec, (qname, tags) in zip(raw, types):
			assert qname == rec.qname
			assert [t for t, _, _ in tags] == ATTRS[:-1]
			for tag, typ, val in tags:
				if tag in INT_TAGS:
					assert val == rec.tag(tag)
					assert typ == hb.smallest_uint_type(val), (qname, tag, typ, val)
					seen.add(typ)
				else:
					assert typ == {"MD": "Z", "jM": "Bc", "jI": "Bi"}[tag]
		assert seen == ({"C", "S"} if reads.paired else {"C"})


PRESETS = {
	"None": [],
	"Standard": ["NH", "HI", "AS", "nM"],
	"default": ["NH", "HI", "AS", "nM"],
	"All": ["NH", "HI", "AS", "nM", "NM", "MD", "jM", "jI", "MC"],
}


@pytest.mark.parametrize("preset,fmt", [(p, f) for p in ("None", "Standard") for f in FORMATS]
	+ [("default", "sam"), ("default", "sorted"), ("All", "bam"), ("All", "sorted")])
def test_attribute_presets(genome_sjdb, sim, preset, fmt):
	"""Presets (parametersDefault, manual 5.2.2): None writes no attributes,
	Standard and the default write NH HI AS nM, All writes NH HI AS nM NM MD
	jM jI MC ch, where ch appears only on chimeric alignments, so proper
	pairs get NH HI AS nM NM MD jM jI MC. All is checked for BAM only (with
	SAM output STAR rejects ch). Unmapped records written by
	--outSAMunmapped Within carry NH HI AS nM uT whatever the preset (pins the
	current output; the manual documents only uT)."""
	pairs = sim.transcript_pairs(20) + ReadSet([hb.unmappable_pair("unm%d" % i) for i in range(2)])
	args = ["--outSAMunmapped", "Within"] + sh.OUTPUT_FORMATS[fmt]
	if preset != "default":
		args += ["--outSAMattributes", preset]
	aln = run_star(genome_sjdb, pairs, args).alignments()
	assert len(aln) == 2 * len(pairs)
	for rec in aln:
		if rec.qname.startswith("unm"):
			assert list(rec.tags) == ["NH", "HI", "AS", "nM", "uT"]
			assert (rec.tag("NH"), rec.tag("HI"), rec.tag("uT")) == (0, 0, "0")
		else:
			assert list(rec.tags) == PRESETS[preset]


def _ladder_reads(sim):
	reads = [sim.genomic_read(*sim.unique_locus(), strand=sim.rng.choice("+-")) for _ in range(3)]
	for seg in ("dup2", "dup3", "dup4", "dup5", "dup6"):
		for _ in range(2):
			reads.append(sim.segment_read(seg, sim.rng.randint(0, 50), strand=sim.rng.choice("+-")))
	return ReadSet(reads)


@pytest.mark.parametrize("fmt", FORMATS)
def test_mapq_ladder_and_hi_start(tiny_genome, tmp_path_factory, fmt):
	"""Reads with 1 to 6 loci (helpers_b tiny genome). MAPQ is
	--outSAMmapqUnique for unique reads and int(-10*log10(1-1/NH)) for
	multimappers (manual 5.2.1): 3 for NH 2, 1 for NH 3 and 4, 0 for NH 5 and
	6. HI counts from --outSAMattrIHstart (parametersDefault: start value of
	the HI attribute). Run with the defaults (255, HI from 1) and with
	--outSAMmapqUnique 60 --outSAMattrIHstart 0. Each record is at one of the
	read's loci and the loci are all reported."""
	g = hb.tiny_synth(tmp_path_factory)
	sim = sh.ReadSimulator(g, seed=2401)
	reads = _ladder_reads(sim)
	ladder = {2: 3, 3: 1, 4: 1, 5: 0, 6: 0}
	for unique_q, hi0, extra in ((255, 1, []), (60, 0, ["--outSAMmapqUnique", "60", "--outSAMattrIHstart", "0"])):
		by_name = run_star(tiny_genome, reads, extra + sh.OUTPUT_FORMATS[fmt]).alignments().by_name()
		assert sorted(by_name) == sorted(reads.by_name)
		for name, recs in by_name.items():
			truth = reads.by_name[name]
			n = truth.extra.get("n_loci", 1)
			loci = truth.extra.get("loci", [(truth.contig, truth.pos, truth.strand)])
			assert sorted((r.rname, r.pos, r.strand) for r in recs) == sorted(loci)
			assert sorted(r.tag("HI") for r in recs) == list(range(hi0, hi0 + n))
			for rec in recs:
				assert rec.tag("NH") == n
				assert rec.mapq == (unique_q if n == 1 else ladder[n])


def _strand_reads(sim):
	return ReadSet([
		sim.transcript_read("Tmix", 66, name="mix"),
		sim.transcript_read("Tnc", 50, name="noncanonical"),
		sim.transcript_read("Tgtag", 50, sense=False, name="gtag"),
		sim.transcript_read("Tctac", 50, name="ctac"),
		sim.genomic_read(*sim.unique_locus(), strand="+", name="plain_fwd"),
		sim.genomic_read(*sim.unique_locus(), strand="-", name="plain_rev"),
	])


@pytest.mark.parametrize("fmt", FORMATS)
def test_strand_field_intron_motif(tiny_genome, tmp_path_factory, fmt):
	"""--outSAMstrandField intronMotif (parametersDefault: "reads with
	inconsistent and/or non-canonical introns are filtered out"; manual
	5.2.3: alignments with undefined strand are suppressed, XS is added to
	all spliced alignments). Reads: Tmix crosses a GT/AG and a CT/AC intron
	(blocks 34/32/34), Tnc a non-canonical intron (50/50), Tgtag and Tctac
	one canonical intron each, plus two unspliced reads. Both runs use
	--outFilterIntronStrands None and --scoreGapNoncan 0 so that, without
	intronMotif, the mixed and non-canonical reads align at their truth
	(jM 1,2 and 0). With intronMotif those two lose their spliced alignment
	and every other alignment scores below 0.66 of the read length, so they
	are unmapped (uT 1); Tgtag gets XS:A:+ although the read is on the -
	strand, Tctac XS:A:- although it is on +, and unspliced reads get no XS."""
	g = hb.tiny_synth(tmp_path_factory)
	reads = _strand_reads(sh.ReadSimulator(g, seed=2402))
	base = ["--outFilterIntronStrands", "None", "--scoreGapNoncan", "0", "--outSAMunmapped", "Within",
		"--outSAMattributes", "NH", "HI", "AS", "nM", "jM"] + sh.OUTPUT_FORMATS[fmt]

	plain = {r.qname: r for r in run_star(tiny_genome, reads, base).records()}
	assert sorted(plain) == sorted(reads.by_name)
	for name, rec in plain.items():
		assert record_truth(rec) == reads.by_name[name].truth()
		assert "XS" not in rec.tags
	assert plain["mix"].tag("jM") == [1, 2]
	assert plain["noncanonical"].tag("jM") == [0]

	motif = {r.qname: r for r in run_star(tiny_genome, reads, base + ["--outSAMstrandField", "intronMotif"]).records()}
	assert sorted(motif) == sorted(reads.by_name)
	for name in ("mix", "noncanonical"):
		rec = motif[name]
		assert rec.is_unmapped and rec.tag("uT") == "1"
	expected_xs = {"gtag": "+", "ctac": "-"}
	for name in ("gtag", "ctac", "plain_fwd", "plain_rev"):
		rec = motif[name]
		assert record_truth(rec) == reads.by_name[name].truth()
		assert rec.tag("XS") == expected_xs.get(name)
		if name in expected_xs:
			assert list(rec.tags) == ["NH", "HI", "AS", "nM", "jM", "XS"] and rec.tag_types["XS"] == "A"
	assert reads.by_name["gtag"].strand == "-" and reads.by_name["ctac"].strand == "+"


@pytest.mark.parametrize("fmt", FORMATS)
def test_read_groups(genome_sjdb, sim, fmt):
	"""--outSAMattrRGline (parametersDefault): three comma-separated RG lines
	for three comma-separated read files, one with extra fields and a quoted
	value with spaces, give three @RG lines (fields tab-separated) and RG:Z
	of each read's file on every record, unmapped ones (--outSAMunmapped
	Within) included. One RG line for two files puts its ID on all records;
	RG named in --outSAMattributes is written at that position, otherwise it
	is appended after the requested attributes (Parameters_samAttributes.cpp)."""
	a = sim.transcript_reads(10)
	b = sim.transcript_reads(10) + ReadSet([sim.random_read() for _ in range(3)])
	c = sim.transcript_reads(10)
	rg = ["--outSAMattrRGline", "ID:grpA", "SM:alpha", "DS:first group", ",", "ID:grpB", ",", "ID:grpC", "PL:ILLUMINA"]
	aln = run_star(genome_sjdb, [a, b, c], rg + ["--outSAMunmapped", "Within"] + sh.OUTPUT_FORMATS[fmt]).alignments()
	assert aln.header_lines("@RG") == ["@RG\tID:grpA\tSM:alpha\tDS:first group", "@RG\tID:grpB",
		"@RG\tID:grpC\tPL:ILLUMINA"]
	assert len(aln) == 33
	for rec in aln:
		group = "grpA" if rec.qname in a.by_name else "grpB" if rec.qname in b.by_name else "grpC"
		assert rec.tag("RG") == group and rec.tag_types["RG"] == "Z"
		assert list(rec.tags) == ["NH", "HI", "AS", "nM"] + (["uT"] if rec.is_unmapped else []) + ["RG"]

	aln = run_star(genome_sjdb, [a, c], ["--outSAMattrRGline", "ID:only", "--outSAMattributes", "NH", "RG", "HI"]
		+ sh.OUTPUT_FORMATS[fmt]).alignments()
	assert aln.header_lines("@RG") == ["@RG\tID:only"]
	assert len(aln) == 20
	for rec in aln:
		assert list(rec.tags) == ["NH", "RG", "HI"] and rec.tag("RG") == "only"
