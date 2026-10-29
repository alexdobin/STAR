"""Area D: STARsolo output options (--soloOutFileNames,
--soloOutFormatFeaturesGeneField3, --soloCellReadStats) and a barcode inside
the cDNA mate (--soloBarcodeMate)."""

import pytest

import helpers_d as hd
import star_helpers as sh
from star_helpers import ReadSet, SimRead, run_star

STAT_NAMES = ["cbMatch", "cbPerfect", "cbMMunique", "cbMMmultiple", "genomeU", "genomeM", "featureU", "featureM",
	"exonic", "intronic", "exonicAS", "intronicAS", "mito", "countedU", "countedM",
	"nUMIunique", "nGenesUnique", "nUMImulti", "nGenesMulti"]


@pytest.fixture(scope="module")
def out_reads(synth_genome, tenx_whitelist):
	"""Cell wl[0]: 3 exonic G01 reads, 1 intergenic, 1 unmapped (random
	cDNA); cell wl[1]: 2 exonic G08 reads and 1 G03 read whose CB has one
	substitution; 1 G07 read whose CB matches no whitelist CB. Every UMI is
	new."""
	wl = tenx_whitelist[0]
	sim = hd.module_sim(synth_genome, __name__, "out_reads")
	pool = hd.UmiPool(sim)
	pairs = []

	def add(kind, cdna, cb, cb_read=None, gene=None):
		p = hd.tenx_pair(sim, cdna, cb, pool.take(), cb_read=cb_read, gene=gene)
		p.extra["kind"] = kind
		pairs.append(p)

	for i in range(3):
		add("exon", hd.exonic_read(sim, "T01a", from_end=150 + 10 * i), wl[0], gene="G01")
	add("intergenic", hd.genomic_read(sim, "chr3", 18101, "+"), wl[0])
	add("unmapped", sim.random_read(90, name=sim.name("cell")), wl[0])
	for i in range(2):
		add("exon", hd.exonic_read(sim, "T08", from_end=150 + 10 * i), wl[1], gene="G08")
	add("mm", hd.exonic_read(sim, "T03", from_end=200), wl[1], cb_read=hd.substitute(wl[1], 7), gene="G03")
	far = sim.random_seq(16)
	while min(sh.hamming(far, w) for w in wl) < 2:
		far = sim.random_seq(16)
	add("far", hd.exonic_read(sim, "T07", from_end=200), None, cb_read=far, gene="G07")
	return ReadSet(pairs, kind="solo")


def test_cell_read_stats(genome_sjdb, out_reads, tenx_whitelist):
	"""--soloCellReadStats Standard writes Gene/CellReads.stats: a header,
	a CBnotInPasslist line and one line per whitelist CB with reads, each
	counting that CB's reads per category (parametersDefault "Output reads
	statistics for each CB"): cbPerfect / cbMMunique exact or one-mismatch
	CB, genomeU unique genome alignment (the unmapped read has none),
	featureU one gene, exonic inside a transcript, countedU counted in the
	matrix, and the cell's unique UMIs and genes. The read without a
	whitelist CB is on the CBnotInPasslist line, mapped and genic but not
	counted."""
	wl = tenx_whitelist[0]
	run = run_star(genome_sjdb, out_reads, hd.solo_args(tenx_whitelist[1], "--soloCellReadStats", "Standard"))
	lines = [l.split("\t") for l in (run.solo_dir("Gene") / "CellReads.stats").read_text().splitlines()]
	assert lines[0] == ["CB"] + STAT_NAMES
	rows = {l[0]: dict(zip(STAT_NAMES, (int(x) for x in l[1:]))) for l in lines[1:]}
	assert set(rows) == {"CBnotInPasslist", wl[0], wl[1]}
	want = {}
	for cb in (wl[0], wl[1]):
		reads = [p for p in out_reads if p.extra["cb"] == cb]
		genic = [p for p in reads if p.extra["gene"]]
		want[cb] = dict.fromkeys(STAT_NAMES, 0)
		want[cb].update({"cbMatch": len(reads), "cbPerfect": sum(p.extra["kind"] != "mm" for p in reads),
			"cbMMunique": sum(p.extra["kind"] == "mm" for p in reads),
			"genomeU": sum(p.extra["kind"] != "unmapped" for p in reads),
			"featureU": len(genic), "exonic": len(genic), "countedU": len(genic), "nUMIunique": len(genic),
			"nGenesUnique": len({p.extra["gene"] for p in genic})})
		assert rows[cb] == want[cb], cb
	no_cb = rows["CBnotInPasslist"]
	assert (no_cb["cbPerfect"], no_cb["genomeU"], no_cb["featureU"], no_cb["exonic"], no_cb["countedU"]) == (0, 1, 1, 1, 0)


def test_out_file_names_and_field3(genome_sjdb, out_reads, tenx_whitelist):
	"""--soloOutFileNames MySolo/ genes.tsv cells.tsv counts.mtx renames the
	output directory and the three matrix files (Barcodes.stats moves with
	the directory); --soloOutFormatFeaturesGeneField3 - drops the third
	features column (parametersDefault)."""
	run = run_star(genome_sjdb, out_reads, hd.solo_args(tenx_whitelist[1], "--soloOutFileNames", "MySolo/", "genes.tsv",
		"cells.tsv", "counts.mtx", "--soloOutFormatFeaturesGeneField3", "-"))
	raw = run.outdir / "MySolo" / "Gene" / "raw"
	assert not run.exists("Solo.out")
	assert sorted(p.name for p in raw.iterdir()) == ["cells.tsv", "counts.mtx", "genes.tsv"]
	assert (run.outdir / "MySolo" / "Barcodes.stats").exists()
	genes = [l.split("\t") for l in (raw / "genes.tsv").read_text().splitlines()]
	assert genes == [[g.gene_id, g.name] for g in run.genome.synth.annotated_genes()]
	cells = (raw / "cells.tsv").read_text().split()
	assert cells == sorted(tenx_whitelist[0])
	_, _, entries = sh.parse_mtx(raw / "counts.mtx")
	got = {(genes[r - 1][0], cells[c - 1]): v for (r, c), v in entries.items()}
	assert got == sh.expected_solo_counts([p for p in out_reads if p.extra["gene"] and p.extra["cb"]])


def test_field3_custom_string(genome_sjdb, out_reads, tenx_whitelist):
	"""--soloOutFormatFeaturesGeneField3 'Custom Type' is written as the
	third features.tsv column of every gene."""
	run = run_star(genome_sjdb, out_reads, hd.solo_args(tenx_whitelist[1], "--soloOutFormatFeaturesGeneField3",
		"Custom Type"))
	assert {f[2] for f in run.solo("Gene").features} == {"Custom Type"}


def test_barcode_inside_cdna_mate(genome_sjdb, out_reads, tenx_whitelist):
	"""docs/STARsolo.md "Barcode and cDNA on the same mate": with
	--soloBarcodeMate 1 --clip5pNbases 28 single-end reads made of CB+UMI
	followed by the cDNA give the same matrix as separate barcode reads,
	and the alignments soft-clip the 28 barcode bases."""
	pairs = [p for p in out_reads if p.extra["kind"] in ("exon", "mm")]
	se = ReadSet([SimRead(p.name, p.r2.seq + p.r1.seq, p.r2.qual + p.r1.qual) for p in pairs])
	run = run_star(genome_sjdb, se, hd.solo_args(tenx_whitelist[1], "--soloBarcodeMate", "1", "--clip5pNbases", "28"))
	assert run.solo("Gene").counts == sh.expected_solo_counts(pairs)
	by_name = {p.name: p for p in pairs}
	for rec in run.records():
		r1 = by_name[rec.qname].r1
		ops = sh.cigar_ops(r1.cigar)
		ops = [(28, "S")] + ops if r1.strand == "+" else ops + [(28, "S")]
		assert (rec.rname, rec.pos, rec.cigar, rec.strand) == (r1.contig, r1.pos, sh.ops_to_cigar(ops), r1.strand)
