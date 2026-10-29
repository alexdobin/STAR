"""Area D: STARsolo UMI collapsing (--soloUMIdedup) and multi-gene UMI
filtering (--soloUMIfiltering), compared with reference implementations of
the documented rules in helpers_d."""

from collections import OrderedDict

import pytest

import helpers_d as hd
from star_helpers import ReadSet, run_star

ALL_DEDUP = list(hd.DEDUP)

# (cell index, transcript, {UMI recipe: reads}); recipes: "u" a fresh UMI,
# "u+3" the UMI of "u" with base 3 changed, "u+3+9" that UMI with base 9
# changed too (2 mismatches from "u").
DEDUP_GROUPS = [
	(0, "T01a", OrderedDict([("u", 5), ("u+5", 1)])),           # 5 vs 1 read, 1MM
	(0, "T03", OrderedDict([("u", 3), ("u+5", 2)])),            # 3 vs 2: UMI-tools collapses, 1MM_Directional not
	(1, "T01a", OrderedDict([("u", 2), ("u+5", 2)])),           # equal counts
	(1, "T07", OrderedDict([("u", 4), ("u+3", 2), ("u+3+9", 1)])),  # chain: 1MM_CR is not transitive
	(2, "T08", OrderedDict([("u", 1), ("u+2+7", 1)])),          # 2 mismatches apart
	(2, "T10", OrderedDict([("u", 3)])),
]

# One cell (index 4): UMI m1 with 3 reads in G01 and 1 in G03, m2 with 2
# reads in G07 and 2 in G08, m3 with 1 read in G01 only.
FILTER_GROUPS = [("m1", "T01a", 3), ("m1", "T03", 1), ("m2", "T07", 2), ("m2", "T08", 2), ("m3", "T01a", 1)]


def _umi(recipe, bases):
	parts = recipe.split("+")
	u = bases[parts[0]]
	for p in parts[1:]:
		u = hd.substitute(u, int(p))
	return u


def _exonic(sim, tx_id):
	t = sim.g.transcripts[tx_id]
	offs = sim.valid_offsets(tx_id, 90, 25, region=(t.length - 300, t.length - 90))
	return sim.transcript_read(tx_id, sim.rng.choice(offs), 90, sense=True, name=sim.name("cell"))


@pytest.fixture(scope="module")
def umi_reads(synth_genome, tenx_whitelist):
	wl = tenx_whitelist[0]
	sim = hd.module_sim(synth_genome, __name__, "umi_reads")
	pool = hd.UmiPool(sim)
	pairs = []
	for cell, tx, recipes in DEDUP_GROUPS:
		bases = {"u": pool.take()}
		gene = sim.g.transcripts[tx].gene_id
		for recipe, n in recipes.items():
			u = _umi(recipe, bases)
			assert len(set(u)) > 1
			for _ in range(n):
				pairs.append(hd.tenx_pair(sim, _exonic(sim, tx), wl[cell], u, gene=gene))
	bases = {k: pool.take() for k in ("m1", "m2", "m3")}
	for key, tx, n in FILTER_GROUPS:
		for _ in range(n):
			pairs.append(hd.tenx_pair(sim, _exonic(sim, tx), wl[4], bases[key], gene=sim.g.transcripts[tx].gene_id))
	# shuffled so that reads of one UMI are not adjacent in the input, as
	# STAR's collapsing must not rely on input order
	sim.rng.shuffle(pairs)
	return ReadSet(pairs, kind="solo")


def _expected(reads, method, cells=None):
	"""{(gene, cell): count} from the reference collapsing of each
	(gene, cell)'s UMI read counts."""
	out = {}
	for key, counts in hd.umi_read_counts(reads).items():
		if cells is None or key[1] in cells:
			out[key] = hd.DEDUP[method](counts)
	return out


@pytest.fixture(scope="module")
def all_dedup_run(genome_sjdb, umi_reads, tenx_whitelist):
	return run_star(genome_sjdb, umi_reads, hd.solo_args(tenx_whitelist[1], "--soloUMIdedup", *ALL_DEDUP),
		name="umi_all_dedup")


def test_reference_dedup_values():
	"""The reference rules give the hand-worked counts for the groups
	(1MM_All, UMItools, Directional, Exact, NoDedup, 1MM_CR):
	5+1 -> 1 1 1 2 6 1; 3+2 -> 1 1 2 2 5 1 (3 >= 2*2-1 but 3 < 2*2);
	2+2 -> 1 2 2 2 4 1; chain 4-2-1 -> 1 1 1 3 7 2 (CR moves the 1-read UMI
	to its 2-read neighbour, which itself moves to the 4-read one; the
	correction is not applied twice); 2 mismatches -> 2 2 2 2 2 2."""
	u, v, w = "ACGTACGTACGT", "ACGTAAGTACGT", "ACGTAAGTATGT"
	cases = [({u: 5, v: 1}, [1, 1, 1, 2, 6, 1]), ({u: 3, v: 2}, [1, 1, 2, 2, 5, 1]),
		({u: 2, v: 2}, [1, 2, 2, 2, 4, 1]), ({u: 4, v: 2, w: 1}, [1, 1, 1, 3, 7, 2]),
		({u: 1, w: 1}, [2, 2, 2, 2, 2, 2])]
	for counts, want in cases:
		assert [hd.DEDUP[m](counts) for m in ALL_DEDUP] == want


@pytest.mark.parametrize("method", ALL_DEDUP)
def test_umi_dedup_methods_in_one_run(all_dedup_run, umi_reads, method):
	"""With several --soloUMIdedup types, raw/ holds one
	umiDedup-<type>.mtx per type (and no matrix.mtx), each equal to that
	type's documented collapsing (parametersDefault soloUMIdedup) applied to
	every (gene, cell): 1MM_All counts 1-mismatch connected components,
	1MM_Directional_UMItools merges b into a when count(a) >= 2*count(b)-1,
	1MM_Directional when count(a) >= 2*count(b), Exact counts distinct UMIs,
	NoDedup counts reads, 1MM_CR applies CellRanger's one-pass correction."""
	raw = all_dedup_run.solo_dir("Gene") / "raw"
	assert not (raw / "matrix.mtx").exists()
	got = all_dedup_run.solo("Gene", mtx="umiDedup-%s.mtx" % method)
	assert got.counts == _expected(umi_reads, method)


def test_umi_dedup_first_type_drives_stats(all_dedup_run, umi_reads):
	"""Features.stats yesUMIs and the filtered matrix use the first
	--soloUMIdedup type (1MM_All here); yessubWLmatch_UniqueFeature counts
	every read. The default CellRanger2.2 3000 0.99 10 filter
	(helpers_d.knee_filter) calls every cell here, so filtered/ holds the
	whole 1MM_All matrix."""
	expected = _expected(umi_reads, "1MM_All")
	fs = all_dedup_run.solo_features_stats("Gene")
	assert fs["yesUMIs"] == sum(expected.values())
	assert fs["yessubWLmatch_UniqueFeature"] == len(umi_reads)
	totals = {}
	for (g, cb), n in expected.items():
		totals[cb] = totals.get(cb, 0) + n
	called = hd.knee_filter(totals, 3000, 0.99, 10)
	assert called == set(totals)
	filt = all_dedup_run.solo("Gene", matrix="filtered", mtx="umiDedup-1MM_All.mtx")
	assert filt.barcodes == sorted(called)
	assert filt.counts == expected


@pytest.mark.parametrize("method", ["1MM_CR", "Exact"])
def test_single_umi_dedup_type_writes_matrix_mtx(genome_sjdb, umi_reads, tenx_whitelist, method):
	"""With one --soloUMIdedup type the counts go to raw/matrix.mtx."""
	run = run_star(genome_sjdb, umi_reads, hd.solo_args(tenx_whitelist[1], "--soloUMIdedup", method))
	assert run.solo("Gene").counts == _expected(umi_reads, method)
	assert run.solo_features_stats("Gene")["yesUMIs"] == sum(_expected(umi_reads, method).values())


FILTER_EXPECTED = {
	# parametersDefault soloUMIfiltering:
	# "-": basic filtering only, m1 counts in G01 and G03, m2 in G07 and G08
	"-": {"G01": 2, "G03": 1, "G07": 1, "G08": 1},
	# MultiGeneUMI: the lower-count copy of m1 (1 read in G03) is removed;
	# m2 has no lower-count copy (2 and 2 reads) and is kept in both genes
	"MultiGeneUMI": {"G01": 2, "G07": 1, "G08": 1},
	# (MultiGeneUMI_All removes nothing in STAR 2.7.11b, contradicting
	# parametersDefault; reported separately)
	# MultiGeneUMI_CR (CellRanger >= 3): a multi-gene UMI goes to the gene
	# with the most reads, and is dropped on a tie (m2)
	"MultiGeneUMI_CR": {"G01": 2},
}


@pytest.mark.parametrize("filtering", list(FILTER_EXPECTED))
def test_umi_filtering_multigene(genome_sjdb, umi_reads, tenx_whitelist, filtering):
	"""--soloUMIfiltering on one cell whose UMIs m1 (3 reads G01, 1 read
	G03) and m2 (2 reads G07, 2 reads G08) map to two genes. Cells whose
	UMIs map to one gene keep their dedup counts (1MM_All, or 1MM_CR which
	MultiGeneUMI_CR requires)."""
	wl = tenx_whitelist[0]
	dedup = "1MM_CR" if filtering == "MultiGeneUMI_CR" else "1MM_All"
	run = run_star(genome_sjdb, umi_reads, hd.solo_args(tenx_whitelist[1], "--soloUMIdedup", dedup,
		"--soloUMIfiltering", filtering))
	got = run.solo("Gene").counts
	assert {g: v for (g, cb), v in got.items() if cb == wl[4]} == FILTER_EXPECTED[filtering]
	assert {k: v for k, v in got.items() if k[1] != wl[4]} == _expected(umi_reads, dedup, set(wl[:4]))


def test_multigene_umi_cr_requires_1mm_cr(genome_sjdb, umi_reads, tenx_whitelist):
	"""parametersDefault: MultiGeneUMI_CR "Only works with --soloUMIdedup
	1MM_CR"; with 1MM_All STAR exits with a parameter error."""
	run = run_star(genome_sjdb, umi_reads, hd.solo_args(tenx_whitelist[1], "--soloUMIfiltering",
		"MultiGeneUMI_CR"), check=False)
	assert run.returncode != 0
	assert "MultiGeneUMI_CR only works with --soloUMIdedup 1MM_CR" in run.stderr + run.log_out()
