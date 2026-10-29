"""Area D: STARsolo SAM attributes (CR CY UR UY CB UB GX GN sS sQ sM) in a
coordinate-sorted BAM (docs/STARsolo.md "BAM tags")."""

import pytest

import helpers_d as hd
import star_helpers as sh
from star_helpers import ReadSet, run_star

TAGS = ["NH", "HI", "nM", "AS", "CR", "CY", "UR", "UY", "CB", "UB", "GX", "GN", "sS", "sQ", "sM"]


@pytest.fixture(scope="module")
def tag_reads(synth_genome, tenx_whitelist):
	"""Reads with known barcode outcomes (all cDNA reads exonic and sense
	unless noted):
	exact: CB in the whitelist; umi3 x3 and its 1-mismatch neighbour x1 in
	the same gene and cell (1MM_All collapses them onto the 3-read UMI);
	mm: CB with one substitution, corrected to its whitelist CB, which also
	has an exact read (1MM_multi needs one);
	intergenic: whitelist CB, no gene; nomatch: CB >= 2 mismatches from the
	whitelist; badumi: UMI with an N."""
	wl = tenx_whitelist[0]
	sim = hd.module_sim(synth_genome, __name__, "tag_reads")
	pool = hd.UmiPool(sim)
	pairs = []

	def add(kind, cdna, cb, umi, cb_read=None, umi_read=None, gene=None):
		p = hd.tenx_pair(sim, cdna, cb, umi, cb_read=cb_read, umi_read=umi_read, gene=gene)
		p.extra["kind"] = kind
		pairs.append(p)

	u = pool.take()
	u1 = hd.substitute(u, 4)
	for umi_read in (u, u, u, u1):
		add("exact", hd.exonic_read(sim, "T08", from_end=200 + 7 * len(pairs)), wl[3], u, umi_read=umi_read, gene="G08")
	add("exact", hd.exonic_read(sim, "T01a", from_end=180), wl[3], pool.take(), gene="G01")
	add("exact", hd.exonic_read(sim, "T02", from_end=200), wl[8], pool.take(), gene="G02")
	add("mm", hd.exonic_read(sim, "T03", from_end=200), wl[8], pool.take(), cb_read=hd.substitute(wl[8], 9), gene="G03")
	add("intergenic", hd.genomic_read(sim, "chr3", 18101, "+"), wl[3], pool.take())
	nomatch = sim.random_seq(16)
	while min(sh.hamming(nomatch, w) for w in wl) < 2:
		nomatch = sim.random_seq(16)
	add("nomatch", hd.exonic_read(sim, "T07", from_end=300), None, pool.take(), cb_read=nomatch, gene="G07")
	v = pool.take()
	add("badumi", hd.exonic_read(sim, "T07", from_end=250), wl[3], v, umi_read=v[:6] + "N" + v[7:], gene="G07")
	return ReadSet(pairs, kind="solo")


@pytest.fixture(scope="module")
def tag_run(genome_sjdb, tag_reads, tenx_whitelist):
	return run_star(genome_sjdb, tag_reads, hd.solo_args(tenx_whitelist[1], "--outSAMtype", "BAM",
		"SortedByCoordinate", "--outSAMattributes", *TAGS), name="solo_bam_tags")


def test_solo_bam_tags(tag_run, tag_reads):
	"""Per record: CR/UR the sequenced (raw) CB/UMI, CY/UY their
	qualities, sS/sQ the whole barcode read and its qualities; CB/UB the
	corrected CB and UMI ('-' when the read has no valid CB or UMI; the
	1-mismatch UMI gets the 3-read UMI under 1MM_All); GX/GN the gene ID and
	name of the alignment when it lies in one gene, '-' otherwise, whatever
	the barcode; sM 0 for an exact
	whitelist match, 1 for a single one-mismatch match, negative when the CB
	did not match."""
	genes = {g.gene_id: g.name for g in tag_run.genome.synth.annotated_genes()}
	recs = tag_run.records()
	assert sorted(r.qname for r in recs) == sorted(tag_reads.by_name)
	for rec in recs:
		p = tag_reads.by_name[rec.qname]
		kind, bc, bq = p.extra["kind"], p.r2.seq, p.r2.qual
		assert sh.record_truth(rec) == p.r1.truth()
		assert (rec.tag("CR"), rec.tag("CY"), rec.tag("UR"), rec.tag("UY")) == (bc[:16], bq[:16], bc[16:], bq[16:])
		assert (rec.tag("sS"), rec.tag("sQ")) == (bc, bq)
		valid = kind in ("exact", "mm", "intergenic")
		assert rec.tag("CB") == (p.extra["cb"] if valid else "-")
		assert rec.tag("UB") == (p.extra["umi"] if valid else "-")
		gene = p.extra["gene"]
		assert (rec.tag("GX"), rec.tag("GN")) == ((gene, genes[gene]) if gene else ("-", "-"))
		if kind in ("exact", "intergenic"):
			assert rec.tag("sM") == 0
		if kind == "mm":
			assert rec.tag("sM") == 1
		if kind == "nomatch":
			assert rec.tag("sM") < 0


def test_solo_bam_matrix_agrees(tag_run, tag_reads):
	"""The Gene matrix from the same run: cell wl[3] has G08 (one UMI after
	collapsing) and G01; cell wl[8] has G02 and the corrected G03; the read without a
	whitelist CB and the read with an N in its UMI are not counted."""
	expected = {}
	for p in tag_reads:
		if p.extra["kind"] in ("exact", "mm"):
			key = (p.extra["gene"], p.extra["cb"])
			expected.setdefault(key, set()).add(p.extra["umi"])
	assert tag_run.solo("Gene").counts == {k: len(v) for k, v in expected.items()}


def test_cb_ub_require_sorted_bam(genome_sjdb, tag_reads, tenx_whitelist):
	"""docs/STARsolo.md: CB/UB "require sorted BAM output"; with unsorted
	BAM STAR exits with a parameter error."""
	run = run_star(genome_sjdb, tag_reads, hd.solo_args(tenx_whitelist[1], "--outSAMtype", "BAM", "Unsorted",
		"--outSAMattributes", "NH", "CB", "UB"), check=False)
	assert run.returncode != 0
	assert "can only be output in the sorted BAM file" in run.stderr + run.log_out()
