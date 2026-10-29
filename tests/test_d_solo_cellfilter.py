"""Area D: STARsolo cell calling (--soloCellFilter), filtered/ matrices,
Summary.csv and UMIperCellSorted.txt."""

import pytest

import helpers_d as hd
import star_helpers as sh
from star_helpers import run_star

# cell index: {gene: UMIs}; totals 40 30 20 12 5 3 2 1
CELLS = {
	0: {"G01": 20, "G03": 10, "G07": 10}, 1: {"G01": 10, "G08": 20}, 2: {"G09": 20},
	3: {"G10": 6, "G11": 6}, 4: {"G01": 5}, 5: {"G02": 3}, 6: {"G03": 2}, 7: {"G04": 1},
}
# EmptyDrops_CR set: 12 low-count cells drawing from a common gene mix and
# one below-knee cell (index 20) with a profile of its own
AMBIENT = ["G01", "G02", "G03", "G08", "G09"]


def _cells(wl, ambient=False):
	cells = {wl[i]: d for i, d in CELLS.items()}
	if ambient:
		for i in range(8, 20):
			cells[wl[i]] = {AMBIENT[i % 5]: 1 + (i % 3), AMBIENT[(i + 2) % 5]: 1}
		cells[wl[20]] = {"G04": 4, "G11": 4}
	return cells


@pytest.fixture(scope="module")
def cf_reads(synth_genome, tenx_whitelist):
	sim = hd.module_sim(synth_genome, __name__, "cf_reads")
	return sim.tenx_reads(tenx_whitelist[0], _cells(tenx_whitelist[0]), reads_per_umi=(1, 2))


@pytest.fixture(scope="module")
def ambient_reads(synth_genome, tenx_whitelist):
	sim = hd.module_sim(synth_genome, __name__, "ambient_reads")
	return sim.tenx_reads(tenx_whitelist[0], _cells(tenx_whitelist[0], ambient=True))


def _run(genome, reads, wl_path, cell_filter, threads=None, extra=()):
	return run_star(genome, reads, hd.solo_args(wl_path, "--soloCellFilter", *cell_filter.split()) + list(extra),
		threads=threads)


def _totals(reads):
	totals = {}
	for (g, cb), n in sh.expected_solo_counts(reads).items():
		totals[cb] = totals.get(cb, 0) + n
	return totals


def _q30(quals):
	"""Fraction of bases with Phred >= 30 (quality character >= '?')."""
	s = "".join(quals)
	return sum(c >= "?" for c in s) / float(len(s))


@pytest.mark.parametrize("params", ["4 0.75 6", "4 0.75 3", "3000 0.99 10"])
def test_cellranger22_knee_filtering(genome_sjdb, cf_reads, tenx_whitelist, params):
	"""--soloCellFilter CellRanger2.2 <nExpected> <maxPercentile>
	<maxMinRatio> keeps the cells helpers_d.knee_filter calls (CellRanger
	2.2: robust max = total of the cell ranked round(nExpected*(1-maxPerc)),
	cutoff = round(robust max / maxMinRatio)): 4 0.75 6 -> cutoff 5, five
	cells; 4 0.75 3 -> cutoff 10, four cells; the default 3000 0.99 10 ->
	every cell. filtered/ holds those barcodes in whitelist order, their raw
	counts, and the same features.tsv."""
	run = _run(genome_sjdb, cf_reads, tenx_whitelist[1], "CellRanger2.2 " + params)
	totals = _totals(cf_reads)
	n_exp, perc, ratio = params.split()
	called = hd.knee_filter(totals, int(n_exp), float(perc), float(ratio))
	assert len(called) == {"4 0.75 6": 5, "4 0.75 3": 4, "3000 0.99 10": 8}[params]
	raw, filt = run.solo("Gene"), run.solo("Gene", "filtered")
	assert filt.barcodes == sorted(called)
	assert filt.counts == {k: v for k, v in raw.counts.items() if k[1] in called}
	assert filt.features == raw.features
	assert filt.shape == (len(raw.features), len(called))
	assert run.solo_summary("Gene")["Estimated Number of Cells"] == len(called)


def test_filtered_summary_statistics(genome_sjdb, cf_reads, tenx_whitelist):
	"""Summary.csv with five called cells (CellRanger2.2 4 0.75 6), every
	value from the simulated reads: reads, saturation 1 - UMIs/reads, Q30
	fractions of the barcode (CB+UMI) and cDNA qualities, mapping
	fractions, and for the called cells their unique reads (sum, integer
	mean, median), UMIs, genes per cell and genes detected (medians of an
	odd number of cells). UMIperCellSorted.txt lists every cell's UMI total
	in decreasing order."""
	run = _run(genome_sjdb, cf_reads, tenx_whitelist[1], "CellRanger2.2 4 0.75 6")
	counts = sh.expected_solo_counts(cf_reads)
	totals = _totals(cf_reads)
	called = hd.knee_filter(totals, 4, 0.75, 6)
	reads_per_cell = {}
	for p in cf_reads:
		reads_per_cell[p.extra["cb"]] = reads_per_cell.get(p.extra["cb"], 0) + 1
	genes_per_cell = {cb: sum(1 for (g, c) in counts if c == cb) for cb in called}
	n, in_cells = len(called), sorted(reads_per_cell[cb] for cb in called)
	s = run.solo_summary("Gene")
	assert s["Number of Reads"] == len(cf_reads)
	assert s["Reads With Valid Barcodes"] == 1
	assert s["Sequencing Saturation"] == pytest.approx(1 - sum(totals.values()) / float(len(cf_reads)), abs=1e-6)
	assert s["Q30 Bases in CB+UMI"] == pytest.approx(_q30(p.r2.qual for p in cf_reads), abs=1e-6)
	assert s["Q30 Bases in RNA read"] == pytest.approx(_q30(p.r1.qual for p in cf_reads), abs=1e-6)
	assert s["Reads Mapped to Genome: Unique+Multiple"] == s["Reads Mapped to Genome: Unique"] == 1
	assert s["Reads Mapped to Gene: Unique+Multiple Gene"] == "NoMulti"
	assert s["Reads Mapped to Gene: Unique Gene"] == 1
	assert s["Estimated Number of Cells"] == n
	assert s["Unique Reads in Cells Mapped to Gene"] == sum(in_cells)
	assert s["Fraction of Unique Reads in Cells"] == pytest.approx(sum(in_cells) / float(len(cf_reads)), abs=1e-6)
	assert s["Mean Reads per Cell"] == sum(in_cells) // n
	assert s["Median Reads per Cell"] == in_cells[n // 2]
	umis = sorted(totals[cb] for cb in called)
	assert s["UMIs in Cells"] == sum(umis)
	assert s["Mean UMI per Cell"] == sum(umis) // n
	assert s["Median UMI per Cell"] == umis[n // 2]
	assert s["Mean Gene per Cell"] == sum(genes_per_cell.values()) // n
	assert s["Median Gene per Cell"] == sorted(genes_per_cell.values())[n // 2]
	assert s["Total Gene Detected"] == len({g for (g, c) in counts if c in called})
	got = [int(x) for x in (run.solo_dir("Gene") / "UMIperCellSorted.txt").read_text().split()]
	assert got == sorted(totals.values(), reverse=True)


def test_no_cell_filter(genome_sjdb, cf_reads, tenx_whitelist):
	"""--soloCellFilter None writes no filtered/ directory, no
	UMIperCellSorted.txt and no called-cell lines in Summary.csv; raw/ is
	unchanged."""
	run = _run(genome_sjdb, cf_reads, tenx_whitelist[1], "None")
	d = run.solo_dir("Gene")
	assert not (d / "filtered").exists()
	assert not (d / "UMIperCellSorted.txt").exists()
	assert "Estimated Number of Cells" not in run.solo_summary("Gene")
	assert run.solo("Gene").counts == sh.expected_solo_counts(cf_reads)


def test_top_cells_more_than_detected(genome_sjdb, cf_reads, tenx_whitelist):
	"""--soloCellFilter TopCells 100 with 8 detected cells keeps all of
	them (the top 100 by UMI count). (STAR 2.7.11b keeps N+1 cells for N
	below the number of detected cells, contrary to "followed by the exact
	number of cells"; reported separately.)"""
	run = _run(genome_sjdb, cf_reads, tenx_whitelist[1], "TopCells 100")
	filt = run.solo("Gene", "filtered")
	assert filt.barcodes == sorted(_totals(cf_reads))
	assert filt.counts == sh.expected_solo_counts(cf_reads)


def test_emptydrops_without_ambient_cells_is_knee_filtering(genome_sjdb, cf_reads, tenx_whitelist):
	"""EmptyDrops_CR takes ambient profiles from the cells ranked indMin to
	indMax (45000-90000 by default); with 8 cells there are none, so the
	called cells are those of the knee step with the same first three
	parameters (CellRanger2.2 4 0.75 6)."""
	run = _run(genome_sjdb, cf_reads, tenx_whitelist[1], "EmptyDrops_CR 4 0.75 6 45000 90000 500 0.01 20000 0.01 10000")
	called = hd.knee_filter(_totals(cf_reads), 4, 0.75, 6)
	assert run.solo("Gene", "filtered").barcodes == sorted(called)


def test_emptydrops_calls_extra_cells(genome_sjdb, ambient_reads, tenx_whitelist):
	"""EmptyDrops_CR 4 0.75 3 2 8 2 ...: the knee step calls the four cells
	with >= 10 UMIs; EmptyDrops then tests the cells below the knee with
	>= umiMin=2 UMIs against the ambient profile of the cells ranked 2-8 and
	adds the significant ones (docs/STARsolo.md: it "calls extra cells
	compared to the knee filtering"). Checked: the knee cells are kept,
	extra cells are candidates, filtered counts equal raw counts, the
	result is the same at 1 and 4 threads. The extra set {4, 6, 8, 20}
	pins the current output."""
	wl = tenx_whitelist[0]
	cf = "EmptyDrops_CR 4 0.75 3 2 8 2 0.01 20000 0.01 10000"
	one = _run(genome_sjdb, ambient_reads, tenx_whitelist[1], cf, threads=1)
	# 4 threads on input chunks of 8698 bytes per mate, so several threads map reads
	four = _run(genome_sjdb, ambient_reads, tenx_whitelist[1], cf, threads=4,
		extra=["--limitIObufferSize", "220000", "50000000"])
	totals = _totals(ambient_reads)
	knee = hd.knee_filter(totals, 4, 0.75, 3)
	filt = one.solo("Gene", "filtered")
	called = set(filt.barcodes)
	assert knee <= called
	assert all(totals[cb] >= 2 for cb in called - knee)
	assert filt.counts == {k: v for k, v in one.solo("Gene").counts.items() if k[1] in called}
	assert four.solo("Gene", "filtered").counts == filt.counts
	assert sorted(wl.index(cb) for cb in called - knee) == [4, 6, 8, 20]
	assert one.solo_summary("Gene")["Estimated Number of Cells"] == len(called)
