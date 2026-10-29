"""Area C: --quantMode GeneCounts (ReadsPerGene.out.tab) and TranscriptomeSAM
(Aligned.toTranscriptome.out.bam): overlap boundaries, strandedness
columns, paired-end reads, transcript coordinates and strand flipping,
--quantTranscriptomeSAMoutput, --quantTranscriptomeBAMcompression, and both
quantifications in one run.

ReadsPerGene.out.tab (manual 8): a read is counted for a gene if it overlaps
(1 nt or more) the exons of one and only one gene; column 2 unstranded,
column 3 read (mate 1) on the gene's strand, column 4 opposite strand.
"""

import pytest

import helpers_c as hc
import star_helpers as sh
from star_helpers import ReadSet, revcomp, run_star

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
FLAG_SECONDARY = sh.FLAG_SECONDARY


def test_gene_counts_overlap_boundaries(genome_sjdb, sim):
	"""ReadsPerGene.out.tab for reads placed by hand around exon edges.
	G07 (-, single exon chr1:30001-31000): a + read ending on 30001 overlaps
	1 base (antisense: columns 2 and 4), a + read ending on 30000 overlaps
	nothing; a - read starting on 31000 overlaps 1 base (sense: columns 2
	and 3), one starting on 31001 nothing. A read inside the first T10
	intron (in the G10 gene body, no exon) is no feature. A + read over the
	G05 exon 1 / intron boundary counts for G05. Two sense T01b reads over
	the exon 1-exon 3 junction count for G01, one antisense spliced T08 read
	for G08 (column 4). Two dup3 reads (3 loci) and one dup2 read (2 loci)
	are N_multimapping. N_unmapped counts
	a random read and a rep12 read (12 loci > --outFilterMultimapNmax 10);
	that it includes "mapped to too many loci" reads pins the current
	output (Transcriptome::quantsOutput sums the four Log.final.out unmapped
	and too-many-loci counts)."""
	reads = [sim.genomic_read("chr1", 29902, "+", name="g07_first_base"),
		sim.genomic_read("chr1", 29901, "+", name="g07_before"),
		sim.genomic_read("chr1", 31000, "-", name="g07_last_base"),
		sim.genomic_read("chr1", 31001, "-", name="g07_after"),
		sim.genomic_read("chr2", 30551, "+", name="g10_intron"),
		sim.genomic_read("chr1", 20551, "+", name="g05_exon_intron"),
		sim.transcript_read("T01b", 250), sim.transcript_read("T01b", 270),
		sim.transcript_read("T08", 350, sense=False),
		sim.segment_read("dup3", 50), sim.segment_read("dup3", 300), sim.segment_read("dup2", 500),
		sim.segment_read("rep12", 100), sim.random_read()]
	assert all("N" in r.cigar for r in reads[6:9])
	run = run_star(genome_sjdb, ReadSet(reads), "--quantMode GeneCounts")
	counts = run.gene_counts()
	expected = {"N_unmapped": (2, 2, 2), "N_multimapping": (3, 3, 3), "N_noFeature": (3, 5, 7),
		"N_ambiguous": (0, 0, 0), "G01": (2, 2, 0), "G05": (1, 1, 0), "G07": (2, 1, 1), "G08": (1, 0, 1)}
	assert list(counts)[:4] == ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"]
	assert list(counts)[4:] == [g.gene_id for g in sim.g.annotated_genes()]
	for k, v in counts.items():
		assert v == expected.get(k, (0, 0, 0)), k


# Genes added for test_nested_annotation (chr1:42001-46000 holds no other
# annotation or repeat): GX1 (+) one 4 kb exon; GX2 (-) two exons inside
# GX1's exon; GX3 (+) one exon inside GX1's exon. Transcripts sorted by start
# are TX1, TX2, TX3, so a read past TX3's end is only found in TX1 through
# the running maximum of transcript ends, and exons sorted by start likewise
# need the running maximum of exon ends.
NESTED_GENES = (
	("GX1", "TX1", "+", [(42001, 46000)]),
	("GX2", "TX2", "-", [(43001, 43200), (43601, 43800)]),
	("GX3", "TX3", "+", [(44001, 44300)]),
)


@pytest.fixture(scope="module")
def nested_genome(genome_factory, synth_genome):
	lines = []
	for gid, tid, strand, exons in NESTED_GENES:
		attr = 'gene_id "%s"; transcript_id "%s";' % (gid, tid)
		lines.append(["chr1", "x", "transcript", exons[0][0], exons[-1][1], ".", strand, ".", attr])
		for s, e in (exons if strand == "+" else exons[::-1]):
			lines.append(["chr1", "x", "exon", s, e, ".", strand, ".", attr])
	path = synth_genome.gtf_path.parent / "annotation_nested.gtf"
	hc.write_gtf(path, synth_genome.gtf_text().splitlines() + lines)
	return genome_factory(gtf=str(path))


def test_nested_annotation(nested_genome, sim):
	"""GeneCounts and TranscriptomeSAM with genes nested inside another
	gene's exon (NESTED_GENES). Reads, worked by hand (columns: unstranded,
	read strand = gene strand, opposite):
	long (+, 45001-45100) and long_as (-, 45501-45600) lie in GX1 only, past
	TX3's end: GX1 in columns 2+3 / 2+4, TX1 at offsets 3000 / 3500;
	in_gx3 (+, 44101-44200): GX1 and GX3, ambiguous in columns 2 and 3, no
	feature in column 4; TX1 offset 2100 and TX3 offset 100 (NH 2);
	in_gx2 (-, 43051-43150): ambiguous unstranded, GX2 in column 3, GX1 in
	column 4; TX1 offset 1050 antisense and TX2 (- strand, 5' end at
	43800) offset 200 + 50 = 250 sense;
	between (+, 43301-43400): inside GX2's intron, so GX1 only; TX1 offset
	1300 and no TX2 record."""
	names = (("long", 45001, "+"), ("long_as", 45501, "-"), ("in_gx3", 44101, "+"), ("in_gx2", 43051, "-"),
		("between", 43301, "+"))
	reads = ReadSet([sim.genomic_read("chr1", p, s, name=n) for n, p, s in names])
	run = run_star(nested_genome, reads, "--quantMode GeneCounts TranscriptomeSAM")
	recs = run.alignments().by_name()
	for r in reads:
		assert [sh.record_truth(x) for x in recs[r.name]] == [r.truth()], r.name
	counts = run.gene_counts()
	expected = {"N_noFeature": (0, 1, 3), "N_ambiguous": (2, 1, 0), "GX1": (3, 2, 2), "GX2": (0, 1, 0)}
	assert set(counts) >= {"GX1", "GX2", "GX3"}
	for k, v in counts.items():
		assert v == expected.get(k, (0, 0, 0)), k
	got, nh = _tx_records(run)
	want = {"long": [("TX1", 3001, "100M", 0)], "long_as": [("TX1", 3501, "100M", 16)],
		"in_gx3": [("TX1", 2101, "100M", 0), ("TX3", 101, "100M", 0)],
		"in_gx2": [("TX1", 1051, "100M", 16), ("TX2", 251, "100M", 0)],
		"between": [("TX1", 1301, "100M", 0)]}
	assert got == want
	assert {k: v[0] for k, v in nh.items()} == {k: len(v) for k, v in want.items()}


def _tx_records(run):
	"""{read name: sorted [(rname, pos, cigar, flag without 0x100)]} and
	{read name: [NH values]} from Aligned.toTranscriptome.out.bam, checking
	that exactly one alignment per read is primary (CHANGES.md 2.5.1a: all
	but one randomly selected alignment are non-primary)."""
	out, nh = {}, {}
	for name, recs in run.transcriptome().by_name().items():
		out[name] = sorted((r.rname, r.pos, r.cigar, r.flag & ~FLAG_SECONDARY) for r in recs)
		nh[name] = [r.tag("NH") for r in recs]
		mates = max(1, len({r.mate for r in recs}))
		assert sum(1 for r in recs if not r.is_secondary) == mates, name
	return out, nh


def _clip_read(sim, tx, o, n, sense, name, mismatches=()):
	"""Read covering transcript bases [o, o+100-n) whose last n bases (as
	sequenced) mismatch the transcript bases that would follow, so the
	alignment soft-clips exactly n bases at the read's 3' end. mismatches:
	offsets of extra substitutions inside the aligned part."""
	t = sim.g.transcript_seq(tx)
	if sense:
		tail = "".join(COMP[b] for b in t[o + 100 - n:o + 100])
	else:
		tail = t[o - n:o][::-1]
	return sim.transcript_read(tx, o, length=100 - n, sense=sense, tail=tail, name=name, mismatches=mismatches)


MODES = [None, "BanSingleEnd_BanIndels_ExtendSoftclip", "BanSingleEnd", "BanSingleEnd_ExtendSoftclip"]


@pytest.mark.parametrize("mode", MODES, ids=["default"] + MODES[1:])
def test_transcriptome_sam_output_modes(genome_sjdb, sim, mode):
	"""Aligned.toTranscriptome.out.bam for single-end reads under each
	--quantTranscriptomeSAMoutput (parametersDefault): the default
	BanSingleEnd_BanIndels_ExtendSoftclip drops alignments with indels and
	extends soft-clipped ends (dropping the alignment if the extension has
	more mismatches than --outFilterMismatchNmax 10: 6 internal plus 4
	extended mismatches are kept, 11 extended are not); BanSingleEnd keeps
	indels and soft clips; BanSingleEnd_ExtendSoftclip keeps indels and
	extends soft clips. In transcript coordinates POS is the 5'-most
	transcript base + 1, sense reads are forward and antisense reverse, and
	the CIGAR is reversed for - strand transcripts (T02): a + read
	40M2D60M at chr1:8101 becomes 60M2D40M at T02 offset 8500-8202. Reads
	whose 3' end mismatches over 4 or 11 bases (_clip_read) align as
	96M4S/89M11S on the genome; extended, the sense read keeps POS and the
	antisense read starts 4 bases earlier. (The 4-base tails are shorter
	than --seedMapMin 5, so they are never mapped as seeds; a 5-base tail
	sits on the seedMapMin boundary, where STAR 2.7.11b differs from the
	documentation.) A read in exon 3 of G01 maps to
	T01a and T01b (NH 2); reads over the T01b and T01a-only junctions map
	to one isoform; reads over the unannotated T10skip junction, inside an
	intron, across an exon end, or intergenic have no record. SEQ is the
	read in transcript orientation."""
	g = sim.g
	reads = {}

	def add(r, records):
		reads[r.name] = (r, records)

	extend = mode != "BanSingleEnd"
	indels = mode in ("BanSingleEnd", "BanSingleEnd_ExtendSoftclip")
	for tx, o in (("T02", 470), ("T03", 350)):
		r = sim.transcript_read(tx, o)
		add(r, [(tx, o + 1, "100M", 0)])
	r = sim.transcript_read("T08", 350, sense=False)
	add(r, [("T08", 351, "100M", 16)])
	add(sim.transcript_read("T01a", 500, name="g01_exon3"), [("T01a", 501, "100M", 0), ("T01b", 351, "100M", 0)])
	add(sim.transcript_read("T01b", 250, sense=False, name="t01b_only"), [("T01b", 251, "100M", 16)])
	add(sim.transcript_read("T01a", 400, name="t01a_only"), [("T01a", 401, "100M", 0)])
	add(sim.transcript_read("T10skip", 350, name="t10skip"), [])
	add(sim.genomic_read("chr2", 30551, "+", name="intronic"), [])
	add(sim.genomic_read("chr1", 20551, "+", name="exon_end"), [])
	c, p = sim.unique_locus(avoid_genes=True)
	add(sim.genomic_read(c, p, "-", name="intergenic"), [])
	for tx, o in (("T03", 100), ("T02", 200)):
		for sense in (True, False):
			fl = 0 if sense else 16
			r = _clip_read(sim, tx, o, 4, sense, "clip4_%s_%d" % (tx, fl))
			if extend:
				add(r, [(tx, o + 1 if sense else o - 3, "100M", fl)])
			else:
				add(r, [(tx, o + 1, "96M4S" if sense else "4S96M", fl)])
			r = _clip_read(sim, tx, o + 30, 11, sense, "clip11_%s_%d" % (tx, fl))
			add(r, [] if extend else [(tx, o + 31, "89M11S" if sense else "11S89M", fl)])
			r = _clip_read(sim, tx, o + 60, 4, sense, "clip4mm6_%s_%d" % (tx, fl), mismatches=(13, 26, 39, 52, 65, 78))
			if extend:
				add(r, [(tx, o + 61 if sense else o + 57, "100M", fl)])
			else:
				add(r, [(tx, o + 61, "96M4S" if sense else "4S96M", fl)])
	assert sim.deletion_is_unambiguous("chr1", 11141, 2) and sim.deletion_is_unambiguous("chr1", 8141, 2)
	r = sim.genomic_read("chr1", 11101, "+", deletion=(40, 2), name="del_T03")
	add(r, [("T03", hc.tx_offset_of(g, "T03", 11101) + 1, "40M2D60M", 0)] if indels else [])
	r = sim.genomic_read("chr1", 8101, "+", deletion=(40, 2), name="del_T02")
	add(r, [("T02", hc.tx_offset_of(g, "T02", 8202) + 1, "60M2D40M", 16)] if indels else [])
	r = sim.genomic_read("chr1", 11201, "-", insertion=(40, 3), name="ins_T03")
	add(r, [("T03", hc.tx_offset_of(g, "T03", 11201) + 1, "40M3I57M", 16)] if indels else [])
	r = sim.genomic_read("chr1", 8201, "-", insertion=(40, 3), name="ins_T02")
	add(r, [("T02", hc.tx_offset_of(g, "T02", 8297) + 1, "57M3I40M", 0)] if indels else [])

	rs = ReadSet([v[0] for v in reads.values()])
	args = "--quantMode TranscriptomeSAM" + ("" if mode is None else " --quantTranscriptomeSAMoutput " + mode)
	run = run_star(genome_sjdb, rs, args)
	genomic = run.alignments().by_name()
	for name, (r, _) in reads.items():
		if r.kind != "random":
			assert [sh.record_truth(x) for x in genomic[name]] == [r.truth()], name
	got, nh = _tx_records(run)
	for name, (r, records) in reads.items():
		assert got.get(name, []) == sorted(records), name
		if records:
			assert nh[name] == [len(records)] * len(records), name
	for rec in run.transcriptome():
		r = reads[rec.qname][0]
		assert rec.seq == (revcomp(r.seq) if rec.is_reverse else r.seq), rec.qname
	tseq = g.transcript_seq("T02")
	rec = [x for x in run.transcriptome() if x.qname == rs.records[0].name][0]
	assert rec.seq == tseq[470:570]


def test_gene_counts_and_transcriptome_paired_end(genome_sjdb, sim):
	"""--quantMode TranscriptomeSAM GeneCounts on read pairs (manual 8: both
	quantifications from one run). Gene counts use the strand of mate 1 and
	both mates' overlaps (manual 8: "Both ends of the paired-end read are
	checked"): T02 pairs with mate 1 sense / antisense go to column 3 / 4; a
	pair with mate 1 intergenic and mate 2 in G07 counts for G07; a pair
	with mate 1 in G05 only and mate 2 in G06 only is ambiguous unstranded,
	G05 in column 3 and G06 in column 4. In transcript coordinates each mate
	has POS = its 5'-most transcript base + 1, FLAG paired/proper/first or
	second with the sense mate forward and TLEN = +-insert. A pair inside
	exon 1 of G01 has one alignment per isoform (NH 2). A pair whose mate 2
	is all N aligns as a single-end alignment (mate 1 only, with the
	score/match minimum lowered to 0.3 of the pair length): it is counted
	for G03 but has no transcriptomic record under the default BanSingleEnd
	rule; nor do the G07 and G05/G06 pairs, which do not lie inside one
	transcript."""
	pf = sim.transcript_pair("T02", 100, 300, forward=True, name="pf")
	pr = sim.transcript_pair("T02", 150, 250, forward=False, name="pr")
	p3 = sim.transcript_pair("T03", 350, 300, forward=True, name="p3")
	iso = sim.transcript_pair("T01a", 20, 200, forward=True, name="iso")
	half = sim.transcript_pair("T03", 100, 300, forward=True, name="half")
	half.r2.seq = "N" * 100
	g07 = sim.genomic_pair("chr1", 29801, 300, "+", name="g07half")
	g0506 = sim.genomic_pair("chr1", 20401, 1100, "+", name="g05g06")
	pairs = [pf, pr, p3, iso, half, g07, g0506]
	run = run_star(genome_sjdb, ReadSet(pairs), "--quantMode TranscriptomeSAM GeneCounts "
		"--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3")
	genomic = run.alignments().by_name()
	assert [(x.flag, x.pos) for x in genomic["half"]] == [(73, 11101)]
	counts = run.gene_counts()
	expected = {"N_noFeature": (0, 2, 4), "N_ambiguous": (1, 0, 0), "G01": (1, 1, 0), "G02": (2, 1, 1),
		"G03": (2, 2, 0), "G05": (0, 1, 0), "G06": (0, 0, 1), "G07": (1, 0, 1)}
	for k, v in counts.items():
		assert v == expected.get(k, (0, 0, 0)), k
	recs = {}
	for r in run.transcriptome():
		recs.setdefault(r.qname, []).append(r)
	assert sorted(recs) == ["iso", "p3", "pf", "pr"]
	fields = lambda rs: sorted((r.rname, r.pos, r.flag & ~FLAG_SECONDARY, r.pnext, r.tlen, r.tag("NH")) for r in rs)
	# 99 = paired, proper, mate reverse, first; 147 = paired, proper, reverse, second;
	# 83 = paired, proper, reverse, first; 163 = paired, proper, mate reverse, second
	assert fields(recs["pf"]) == [("T02", 101, 99, 301, 300, 1), ("T02", 301, 147, 101, -300, 1)]
	assert fields(recs["pr"]) == [("T02", 151, 163, 301, 250, 1), ("T02", 301, 83, 151, -250, 1)]
	assert fields(recs["p3"]) == [("T03", 351, 99, 551, 300, 1), ("T03", 551, 147, 351, -300, 1)]
	assert fields(recs["iso"]) == [("T01a", 21, 99, 121, 200, 2), ("T01a", 121, 147, 21, -200, 2),
		("T01b", 21, 99, 121, 200, 2), ("T01b", 121, 147, 21, -200, 2)]
	assert sum(1 for r in recs["iso"] if not r.is_secondary) == 2
	assert len({r.rname for r in recs["iso"] if not r.is_secondary}) == 1


@pytest.fixture(scope="module")
def compression_reads(synth_genome):
	return sh.ReadSimulator(synth_genome, seed=4301).transcript_reads(150)


def _tx_multiset(run):
	return sorted((r.qname, r.rname, r.pos, r.cigar, r.flag & ~FLAG_SECONDARY, r.pnext, r.tlen, r.seq, r.qual,
		tuple(r.tags.items())) for r in run.transcriptome())


@pytest.fixture(scope="module")
def compression_default(genome_sjdb, compression_reads):
	return _tx_multiset(run_star(genome_sjdb, compression_reads, "--quantMode TranscriptomeSAM"))


@pytest.mark.parametrize("level", [0, 1, 9, -2], ids=["0", "1", "9", "minus2"])
def test_transcriptome_bam_compression(genome_sjdb, compression_reads, compression_default, level):
	"""--quantTranscriptomeBAMcompression (parametersDefault: -2 no BAM
	output, 0 no compression, default 1): -2 writes no
	Aligned.toTranscriptome.out.bam (the genomic output is still written);
	0 writes BGZF blocks whose deflate data are stored blocks (BTYPE 0); 1
	and 9 write compressed blocks. The decoded records equal those of a run
	without the option, with the primary flag masked (it is chosen at random,
	CHANGES.md 2.5.1a), and every read has a record on the transcript it
	was drawn from at its offset."""
	run = run_star(genome_sjdb, compression_reads,
		"--quantMode TranscriptomeSAM --quantTranscriptomeBAMcompression %d" % level)
	if level == -2:
		assert not run.exists("Aligned.toTranscriptome.out.bam")
		assert {r.qname for r in run.alignments()} == set(compression_reads.by_name)
		return
	types = hc.bgzf_block_types(run.path("Aligned.toTranscriptome.out.bam"))
	assert types and all((t == 0) == (level == 0) for t in types)
	got = _tx_multiset(run)
	assert got == compression_default
	placed = {(x[0], x[1], x[2]) for x in got}
	for r in compression_reads:
		assert (r.name, r.tx, r.extra["tx_offset"] + 1) in placed, r.name
