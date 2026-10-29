"""Area D: STARsolo feature assignment.

One run per --soloStrand with --soloFeatures Gene GeneFull GeneFull_Ex50pAS
GeneFull_ExonOverIntron SJ Velocyto on reads whose relation to each gene is
known (exonic, spliced, intronic, exon/intron spanning, antisense, gene
overlaps, nested genes). Every read category gets its own cell, so each
cell's counts show how that category was assigned.
"""

from collections import OrderedDict

import pytest

import helpers_d as hd
import star_helpers as sh
from star_helpers import ReadSet, run_star

FEATURES = ("Gene", "GeneFull", "GeneFull_Ex50pAS", "GeneFull_ExonOverIntron")
STRANDS = ("Forward", "Reverse", "Unstranded")

# Each read's relation to each gene whose locus it touches:
# (gene, orientation 's' sense / 'a' antisense, relation, velocyto class).
# relation: 'exon' = inside the exons of a transcript and concordant with its
# junctions; 'exon50' = more than half of the bases exonic but not
# concordant; 'intron' = inside a transcript, at most half exonic; 'partial'
# = overlaps the gene locus but extends past the transcript end.
# velocyto class (velocyto.py rules, docs/STARsolo.md "Velocyto"): 'S'
# exonic in every model (spliced), 'U' intronic or spanning an exon/intron
# boundary (unspliced), 'A' exonic in one model and intronic in another
# (ambiguous).
CATEGORIES = OrderedDict([
	# G07 (-) single exon chr1:30001-31000; G03 (+) exon 1 chr1:11001-11400
	("exon_s", [(("chr1", 30401, "-"), [("G07", "s", "exon", "S")]),
		(("chr1", 11101, "+"), [("G03", "s", "exon", "S")])]),
	("exon_a", [(("chr1", 30401, "+"), [("G07", "a", "exon", "S")]),
		(("chr1", 11101, "-"), [("G03", "a", "exon", "S")])]),
	# across the annotated junctions of T01a (1301-2300) and T03 (11401-11800)
	("splice_s", [(("T01a", 250, True), [("G01", "s", "exon", "S")]),
		(("T03", 350, True), [("G03", "s", "exon", "S")])]),
	("splice_a", [(("T01a", 250, False), [("G01", "a", "exon", "S")]),
		(("T03", 350, False), [("G03", "a", "exon", "S")])]),
	# inside G09's 5 kb intron (chr2:12301-17300), away from GN1/GN2
	("intron_s", [(("chr2", 14101, "+"), [("G09", "s", "intron", "U")]),
		(("chr2", 15101, "+"), [("G09", "s", "intron", "U")])]),
	("intron_a", [(("chr2", 14501, "-"), [("G09", "a", "intron", "U")])]),
	# 50 bases of G03 exon 1 and 40 of its intron
	("span_s", [(("chr1", 11351, "+"), [("G03", "s", "exon50", "U")])]),
	("span_a", [(("chr1", 11351, "-"), [("G03", "a", "exon50", "U")])]),
	# 87 bases of G03 exon 1 and 3 of its intron: too little intron for
	# velocyto.py (MIN_FLANK 5; STAR uses 6, Transcriptome_classifyAlign.cpp),
	# so velocyto counts it as spliced, while Gene does not count it
	("exon_edge", [(("chr1", 11314, "+"), [("G03", "s", "exon50", "S")])]),
	# G01 exon 2 (2301-2450): exonic in T01a, intronic in T01b
	("exon2_T01a", [(("chr1", 2331, "+"), [("G01", "s", "exon", "A")])]),
	("intergenic", [(("chr3", 18001, "+"), [])]),
	# exons of G05 (+) and G06 (-) overlap at chr1:21201-21400
	("overlap", [(("chr1", 21251, "+"), [("G05", "s", "exon", "S"), ("G06", "a", "exon", "S")])]),
	# G07 locus ends at 30001; the read covers 29951-30040 on its strand
	("gene_end", [(("chr1", 29951, "-"), [("G07", "s", "partial", None)])]),
	# unannotated junction 30401-31800 (T10skip) between two G10 exons
	("novel_sj", [(("T10skip", 350, True), [("G10", "s", "exon50", None)])]),
	# GN1 (+) nests in G09's intron: inside GN1 exon 1 (13001-13200)
	("nested_exon", [(("chr2", 13021, "+"), [("GN1", "s", "exon", "S"), ("G09", "s", "intron", "U")])]),
	# 60 bases of GN1 exon 1 and 30 of its intron
	("nested_span", [(("chr2", 13141, "+"), [("GN1", "s", "exon50", "U"), ("G09", "s", "intron", "U")])]),
	# inside GN2 (-) exon 16001-16400 on the + strand: antisense to GN2,
	# sense in G09's intron
	("nested_as", [(("chr2", 16101, "+"), [("GN2", "a", "exon", "S"), ("G09", "s", "intron", "U")])]),
])

# Cells left out of some comparisons, with the reason.
# gene_end: parametersDefault describes GeneFull_ExonOverIntron and
# GeneFull_Ex50pAS as counting "all reads overlapping genes' exons and
# introns", but STAR 2.7.11b counts only reads inside the gene/transcript
# for both (reported separately).
# gene_end, novel_sj: velocyto.py's treatment of reads extending past a
# transcript or with unannotated junctions is not documented for STARsolo.
SKIP = {"GeneFull_ExonOverIntron": {"gene_end"}, "GeneFull_Ex50pAS": {"gene_end"},
	"Velocyto": {"gene_end", "novel_sj"}}


def _ok(orient, strand):
	"""Whether a read in this orientation to a gene counts for it
	(parametersDefault soloStrand: Forward = read strand same as the RNA)."""
	return strand == "Unstranded" or (orient == "s") == (strand == "Forward")


def expected_genes(relations, strand):
	"""{feature: set of genes} for one read, from parametersDefault
	soloFeatures:
	- Gene: genes with a transcript the read is concordant with;
	- GeneFull: genes whose locus the read overlaps;
	- GeneFull_ExonOverIntron: Gene genes if any, else genes containing it;
	- GeneFull_Ex50pAS: sense concordant > (antisense concordant: not
	  counted) > sense >50% exonic > (antisense >50%: not counted) > sense
	  inside a transcript; antisense-only reads are not counted.
	"""
	sense = [(g, rel) for g, o, rel, _ in relations if _ok(o, strand)]
	anti = [(g, rel) for g, o, rel, _ in relations if not _ok(o, strand)]
	out = {}
	out["Gene"] = {g for g, rel in sense if rel == "exon"}
	out["GeneFull"] = {g for g, rel in sense}
	out["GeneFull_ExonOverIntron"] = out["Gene"] or {g for g, rel in sense if rel != "partial"}
	ex50 = set()
	for level in ("exon", "exon50", "intron"):
		s = {g for g, rel in sense if rel == level}
		a = {g for g, rel in anti if rel == level}
		if s:
			ex50 = s
			break
		if a and level != "intron":
			break
	out["GeneFull_Ex50pAS"] = ex50
	return out


def expected_velocyto(relations, strand):
	"""(gene, class) or None: reads touching one gene (on the counted
	strand, inside a transcript) get that gene's class; multi-gene reads are
	not counted."""
	rel = [(g, v) for g, o, r, v in relations if _ok(o, strand) and r != "partial"]
	genes = {g for g, v in rel}
	if len(genes) != 1:
		return None
	return rel[0]


def _build_reads(sim, cells):
	pool = hd.UmiPool(sim)
	pairs = []
	for cat, cb in cells.items():
		for spec, relations in CATEGORIES[cat]:
			name = sim.name("cell")
			if spec[0] in sim.g.transcripts:
				cdna = sim.transcript_read(spec[0], spec[1], 90, sense=spec[2], name=name)
			else:
				cdna = sim.genomic_read(spec[0], spec[1], spec[2], 90, name=name)
			p = hd.tenx_pair(sim, cdna, cb, pool.take())
			p.extra["category"] = cat
			p.extra["relations"] = relations
			pairs.append(p)
	return ReadSet(pairs, kind="solo")


@pytest.fixture(scope="module")
def feature_genome(genome_factory, synth_genome):
	return hd.area_d_genome(genome_factory, synth_genome)


@pytest.fixture(scope="module")
def feature_reads(synth_genome, tenx_whitelist):
	wl = tenx_whitelist[0]
	cells = OrderedDict((cat, wl[i]) for i, cat in enumerate(CATEGORIES))
	sim = hd.module_sim(synth_genome, __name__, "feature_reads")
	return _build_reads(sim, cells), cells


@pytest.fixture(scope="module")
def feature_runs(feature_genome, feature_reads, tenx_whitelist):
	reads, _ = feature_reads
	out = {}
	for strand in STRANDS:
		out[strand] = run_star(feature_genome, reads, hd.solo_args(tenx_whitelist[1], "--soloFeatures",
			"Gene", "GeneFull", "GeneFull_Ex50pAS", "GeneFull_ExonOverIntron", "SJ", "Velocyto",
			"--soloStrand", strand), name="features_" + strand)
	return out


def _expected_matrix(reads, feature, strand):
	counts, stats = {}, {"noNoFeature": 0, "MultiFeature": 0, "counted": 0}
	for p in reads:
		genes = expected_genes(p.extra["relations"], strand)[feature]
		if not genes:
			stats["noNoFeature"] += 1
		elif len(genes) > 1:
			stats["MultiFeature"] += 1
		else:
			key = (genes.pop(), p.extra["cb"])
			counts[key] = counts.get(key, 0) + 1
			stats["counted"] += 1
	return counts, stats


@pytest.mark.parametrize("strand", STRANDS)
@pytest.mark.parametrize("feature", FEATURES)
def test_gene_features_per_strand(feature_runs, feature_reads, feature_genome, strand, feature):
	"""Each gene-level feature's raw matrix holds, per cell, one UMI for
	every read that the documented rule (expected_genes) assigns to exactly
	one gene on the counted strand. Reads with no gene are noNoFeature and
	reads with several genes MultiFeature in Features.stats; with
	--soloMultiMappers Unique the latter are not counted. yesUMIs and
	yesCellBarcodes follow from the matrix. features.tsv lists the GTF genes
	(sorted by ID, which is GTF order here) with field 3 'Gene Expression'."""
	reads, cells = feature_reads
	run = feature_runs[strand]
	expected, stats = _expected_matrix(reads, feature, strand)
	skip = {cells[c] for c in SKIP.get(feature, ())}
	got = run.solo(feature)
	assert {k: v for k, v in got.counts.items() if k[1] not in skip} == \
		{k: v for k, v in expected.items() if k[1] not in skip}
	ids = [g.gene_id for g in feature_genome.synth.annotated_genes()] + [g[0] for g in hd.EXTRA_GENES]
	assert got.feature_ids() == ids
	assert {f[2] for f in got.features} == {"Gene Expression"}
	fs = run.solo_features_stats(feature)
	assert fs["noUnmapped"] == 0
	assert fs["MultiFeature"] == stats["MultiFeature"]
	assert fs["subMultiFeatureMultiGenomic"] == 0
	assert fs["noNoFeature"] + fs["MultiFeature"] + fs["yesWLmatch"] == len(reads)
	assert fs["yesWLmatch"] == fs["yessubWLmatch_UniqueFeature"] == fs["yessubWLmatchExact"]
	assert fs["yesUMIs"] == sum(got.counts.values())
	assert fs["yesCellBarcodes"] == len({cb for _, cb in got.counts})
	if not skip:
		assert fs["noNoFeature"] == stats["noNoFeature"]
		assert fs["yesWLmatch"] == stats["counted"]


@pytest.mark.parametrize("strand", STRANDS)
def test_velocyto_spliced_unspliced_ambiguous(feature_runs, feature_reads, strand):
	"""Velocyto writes spliced.mtx, unspliced.mtx and ambiguous.mtx with
	the same entries: every (gene, cell) with a count in any of them, zero
	where that class has none. Exonic and spliced reads are spliced,
	intronic and exon/intron-spanning reads unspliced, a read exonic in
	T01a but intronic in T01b ambiguous; reads touching two genes on the
	counted strand are not counted (expected_velocyto)."""
	reads, cells = feature_reads
	run = feature_runs[strand]
	skip = {cells[c] for c in SKIP["Velocyto"]}
	expected = {cls: {} for cls in "SUA"}
	for p in reads:
		ev = expected_velocyto(p.extra["relations"], strand)
		if ev is None or p.extra["cb"] in skip:
			continue
		key = (ev[0], p.extra["cb"])
		for cls in "SUA":
			expected[cls].setdefault(key, 0)
		expected[ev[1]][key] += 1
	for cls, mtx in zip("SUA", ("spliced.mtx", "unspliced.mtx", "ambiguous.mtx")):
		got = run.solo("Velocyto", mtx=mtx)
		assert {k: v for k, v in got.counts.items() if k[1] not in skip} == expected[cls], mtx


@pytest.mark.parametrize("strand", STRANDS)
def test_sj_feature_counts_junctions(feature_runs, feature_reads, strand):
	"""The SJ feature counts, per cell, the UMIs of reads crossing each
	junction in SJ.out.tab, whatever --soloStrand says (parametersDefault:
	"splice junctions: reported in SJ.out.tab"). Matrix rows are SJ.out.tab
	lines and SJ.out.tab holds exactly the junctions of the simulated
	spliced reads; SJ/raw/features.tsv is a link to SJ.out.tab. Reads
	without a junction are noNoFeature."""
	reads, _ = feature_reads
	run = feature_runs[strand]
	sj = run.sj()
	truth_sj = {}
	for p in reads:
		for s, e in p.r1.introns:
			truth_sj.setdefault((p.r1.contig, s, e), []).append(p.extra["cb"])
	assert sorted((r.contig, r.start, r.end) for r in sj) == sorted(truth_sj)
	row = {(r.contig, r.start, r.end): i + 1 for i, r in enumerate(sj)}
	expected = {}
	for key, cbs in truth_sj.items():
		for cb in cbs:
			k = (row[key], cb)
			expected[k] = expected.get(k, 0) + 1
	raw = run.solo_dir("SJ") / "raw"
	barcodes = (raw / "barcodes.tsv").read_text().split()
	nr, nc, entries = sh.parse_mtx(raw / "matrix.mtx")
	assert (nr, nc) == (len(sj), len(barcodes))
	assert {(r, barcodes[c - 1]): v for (r, c), v in entries.items()} == expected
	assert (raw / "features.tsv").read_text() == run.path("SJ.out.tab").read_text()
	fs = run.solo_features_stats("SJ")
	n_spliced = sum(len(cbs) for cbs in truth_sj.values())
	assert fs["noNoFeature"] == len(reads) - n_spliced
	assert fs["yesWLmatch"] == fs["yesUMIs"] == n_spliced
	assert fs["yesCellBarcodes"] == len({cb for cbs in truth_sj.values() for cb in cbs})
