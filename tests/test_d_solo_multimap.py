"""Area D: STARsolo multi-gene reads (--soloMultiMappers Uniform, PropUnique,
Rescue, EM) and the UniqueAndMult-*.mtx matrices.

The area-D annotation puts gene GR1 (+) on the dup2 copy at chr1:26001 and
GR2 (-) on the reverse-complemented copy at chr2:8001, each with 200 unique
bases, so a sense read from dup2 aligns to two loci in two genes while reads
from the unique ends belong to one gene.
"""

from collections import OrderedDict

import pytest

import helpers_d as hd
import star_helpers as sh
from star_helpers import ReadSet, run_star

METHODS = ("Uniform", "PropUnique", "Rescue", "EM")

# cell index: (unique GR1 UMIs, unique GR2 UMIs, reads per multi-gene UMI)
CELLS = OrderedDict([
	(2, (3, 1, [1, 1])),
	(6, (0, 0, [1])),
	(10, (2, 0, [1])),
	(14, (0, 2, [2])),
	# EM shrinks GR2's share by a factor 0.8 per step and stops with it at
	# ~0.04, after passing 0.09: tests the stopping rule and the zeroing of
	# shares below 0.01
	(18, (1, 0, [1] * 4)),
])


@pytest.fixture(scope="module")
def mm_genome(genome_factory, synth_genome):
	return hd.area_d_genome(genome_factory, synth_genome)


@pytest.fixture(scope="module")
def mm_reads(synth_genome, tenx_whitelist):
	wl = tenx_whitelist[0]
	sim = hd.module_sim(synth_genome, __name__, "mm_reads")
	pool = hd.UmiPool(sim)
	pairs = []

	def add(contig, pos, strand, cb, umi, genes):
		cdna = sim.genomic_read(contig, pos, strand, 90, name=sim.name("cell"))
		p = hd.tenx_pair(sim, cdna, cb, umi)
		p.extra["genes"] = genes
		pairs.append(p)

	for ci, (n1, n2, multi) in CELLS.items():
		cb = wl[ci]
		for _ in range(n1):
			add("chr1", sim.rng.randint(25801, 25910), "+", cb, pool.take(), ("GR1",))
		for _ in range(n2):
			add("chr2", sim.rng.randint(9501, 9610), "-", cb, pool.take(), ("GR2",))
		for n in multi:
			umi = pool.take()
			for _ in range(n):
				add("chr1", sim.rng.randint(26101, 27300), "+", cb, umi, ("GR1", "GR2"))
	return ReadSet(pairs, kind="solo")


def _solo_mm_args(wl_path):
	return hd.solo_args(wl_path, "--soloMultiMappers", "Unique", *METHODS)


@pytest.fixture(scope="module")
def mm_run(mm_genome, mm_reads, tenx_whitelist):
	return run_star(mm_genome, mm_reads, _solo_mm_args(tenx_whitelist[1]), name="multimap")


def _per_cell(reads):
	unique, multi = {}, {}
	for p in reads:
		cb = p.extra["cb"]
		if len(p.extra["genes"]) == 1:
			unique.setdefault(cb, {}).setdefault(p.extra["genes"][0], set()).add(p.extra["umi"])
		else:
			multi.setdefault(cb, {})[p.extra["umi"]] = p.extra["genes"]
	unique = {cb: {g: len(u) for g, u in d.items()} for cb, d in unique.items()}
	multi = {cb: list(d.values()) for cb, d in multi.items()}
	return unique, multi


def test_multigene_reads_alignment_truth(mm_run, mm_reads):
	"""The dup2 reads align to their two copies (NH 2), the reads from the
	unique ends to one locus, which is what makes them multi-gene and
	unique-gene reads."""
	nh = {}
	for rec in mm_run.records():
		nh[rec.qname] = rec.tag("NH")
	for p in mm_reads:
		assert nh[p.name] == len(p.extra["genes"])


def test_unique_matrix_ignores_multigene_reads(mm_run, mm_reads):
	"""raw/matrix.mtx holds unique-gene UMIs only ("Unique: count only
	reads that map to unique genes"). Features.stats counts every dup2
	read as MultiFeature and subMultiFeatureMultiGenomic; yesWLmatch
	counts unique and multi-gene reads, yessubWLmatch_UniqueFeature and
	yesUMIs the unique ones. Summary.csv reports the multi-gene fraction."""
	unique, multi = _per_cell(mm_reads)
	expected = {(g, cb): n for cb, d in unique.items() for g, n in d.items()}
	assert mm_run.solo("Gene").counts == expected
	n_multi = sum(len(p.extra["genes"]) > 1 for p in mm_reads)
	fs = mm_run.solo_features_stats("Gene")
	assert fs["MultiFeature"] == fs["subMultiFeatureMultiGenomic"] == n_multi
	assert fs["yesWLmatch"] == len(mm_reads)
	assert fs["yessubWLmatch_UniqueFeature"] == len(mm_reads) - n_multi
	assert fs["yesUMIs"] == sum(expected.values())
	summary = mm_run.solo_summary("Gene")
	assert summary["Reads Mapped to Gene: Unique+Multiple Gene"] == pytest.approx(1.0)
	assert summary["Reads Mapped to Gene: Unique Gene"] == pytest.approx(
		(len(mm_reads) - n_multi) / float(len(mm_reads)), abs=1e-6)
	assert summary["Reads Mapped to Genome: Unique"] == pytest.approx(
		(len(mm_reads) - n_multi) / float(len(mm_reads)), abs=1e-6)


@pytest.mark.parametrize("method", METHODS)
def test_unique_and_mult_matrices(mm_run, mm_reads, method):
	"""UniqueAndMult-<method>.mtx = unique counts + the cell's multi-gene
	UMIs distributed as docs/STARsolo.md "Multi-gene reads" describes
	(helpers_d.distribute_multigene). Hand-worked for the cell with GR1 3,
	GR2 1 unique UMIs and 2 multi-gene UMIs: Uniform 4/2, PropUnique
	4.5/1.5, Rescue (weights 3+1 : 1+1) 4.333/1.667, EM converging to the
	maximum-likelihood 4.5/1.5. Genes that get a zero share still have an
	entry."""
	unique, multi = _per_cell(mm_reads)
	expected = {}
	for cb in set(unique) | set(multi):
		u = unique.get(cb, {})
		add = hd.distribute_multigene(u, multi.get(cb, []), method)
		for g in set(u) | set(add):
			expected[(g, cb)] = u.get(g, 0) + add.get(g, 0.0)
	got = mm_run.solo("Gene", mtx="UniqueAndMult-%s.mtx" % method).counts
	assert sorted(got) == sorted(expected)
	for k, v in expected.items():
		assert got[k] == pytest.approx(v, rel=2e-5, abs=1e-6), (method, k)


def test_reference_distribution_hand_values():
	"""distribute_multigene on the documented example cell."""
	u, m = {"A": 3, "B": 1}, [("A", "B"), ("A", "B")]
	assert hd.distribute_multigene(u, m, "Uniform") == {"A": 1.0, "B": 1.0}
	assert hd.distribute_multigene(u, m, "PropUnique") == {"A": 1.5, "B": 0.5}
	r = hd.distribute_multigene(u, m, "Rescue")
	assert (r["A"], r["B"]) == pytest.approx((4.0 / 3, 2.0 / 3))
	em = hd.distribute_multigene(u, m, "EM")
	assert abs(em["A"] - 1.5) < 0.02 and abs(em["A"] + em["B"] - 2) < 1e-9


# STAR reads at most limitIObufferSize[0] / 2 - 2 * 651 - 2 * 50000 bytes of
# each mate per input chunk (Parameters.cpp, IncludeDefine.h): 8698 here, so
# the read set below is split into several chunks mapped by different
# threads. With the default 30 MB buffer every test read set fits in one
# chunk and only one thread gets reads.
SMALL_CHUNKS = ["--limitIObufferSize", "220000", "50000000"]


@pytest.fixture(scope="module")
def chunk_reads(synth_genome, tenx_whitelist, mm_reads):
	"""(single-gene 10x reads, those plus mm_reads): 16 cells with 4 of
	G01 G02 G03 G07 G08 G10 each, 1-3 UMIs per gene and 1-3 reads per UMI,
	and 4 reads whose CB has one substitution; the combined set is in random
	order."""
	wl = tenx_whitelist[0]
	sim = hd.module_sim(synth_genome, __name__, "chunk_reads")
	sim._n = 50000  # read names must not repeat those of mm_reads
	genes = ("G01", "G02", "G03", "G07", "G08", "G10")
	cells = OrderedDict((wl[i], OrderedDict((genes[(i + j) % 6], 1 + (i + j) % 3) for j in range(4))) for i in range(16))
	base = sim.tenx_reads(wl, cells, reads_per_umi=(1, 3), cb_errors={wl[0]: 2, wl[5]: 2})
	pairs = base.records + mm_reads.records
	sim.rng.shuffle(pairs)
	return base, ReadSet(pairs, kind="solo")


def test_solo_outputs_do_not_depend_on_input_chunks(mm_genome, chunk_reads, mm_reads, tenx_whitelist):
	"""Several input chunks mapped by 4 threads (SMALL_CHUNKS) give the
	Solo.out of one chunk mapped by 1 thread, so the merging of per-thread
	barcode counts, read records and statistics is checked (and EM, which
	sums floating-point shares, does not depend on the order reads arrive
	in): every file under Solo.out (Barcodes.stats; for Gene, GeneFull, SJ and Velocyto
	Features.stats, Summary.csv, raw and filtered matrices including the
	UniqueAndMult-*.mtx of every --soloMultiMappers method) is byte-identical, the Gene matrix equals the
	truth (10x reads plus the unique-gene part of mm_reads; 1-mismatch CBs
	are corrected only if their whitelist CB's exact reads from all chunks
	are counted), and the coordinate-sorted BAM records with CB/UB/GX/GN are
	the same."""
	base, reads = chunk_reads
	args = hd.solo_args(tenx_whitelist[1], "--soloFeatures", "Gene", "GeneFull", "SJ", "Velocyto",
		"--soloMultiMappers", "Unique", "Uniform", "PropUnique", "Rescue", "EM", "--outSAMtype", "BAM", "SortedByCoordinate",
		"--outSAMattributes", "NH", "HI", "AS", "nM", "CB", "UB", "GX", "GN")
	one = run_star(mm_genome, reads, args, threads=1)
	four = run_star(mm_genome, reads, args + SMALL_CHUNKS, threads=4)
	files = sorted(p.relative_to(one.outdir) for p in (one.outdir / "Solo.out").rglob("*") if p.is_file())
	assert files == sorted(p.relative_to(four.outdir) for p in (four.outdir / "Solo.out").rglob("*") if p.is_file())
	assert len(files) > 20
	for f in files:
		assert (one.outdir / f).read_bytes() == (four.outdir / f).read_bytes(), str(f)
	expected = sh.expected_solo_counts(base)
	unique, _ = _per_cell(mm_reads)
	expected.update({(g, cb): n for cb, d in unique.items() for g, n in d.items()})
	assert one.solo("Gene").counts == expected
	assert sh.record_multiset(four.records(), bam_types=True) == sh.record_multiset(one.records(), bam_types=True)
	assert one.log_final()["Number of input reads"] == len(reads)
