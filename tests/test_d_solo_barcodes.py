"""Area D: STARsolo cell barcode and UMI checks (--soloCBmatchWLtype,
--soloCBwhitelist None, --soloBarcodeReadLength, barcode geometry, gzipped
input), Barcodes.stats and the barcode part of Features.stats/Summary.csv."""

from collections import OrderedDict

import pytest

import helpers_d as hd
import star_helpers as sh
from star_helpers import ReadSet, SimPair, SimRead, run_star

MODES = ("Exact", "1MM", "1MM_multi", "1MM_multi_pseudocounts", "1MM_multi_Nbase_pseudocounts")
HIGH, LOW, TOP = "F", "#", "I"  # Phred 37, 2, 40


def _far_from(b, others, d=3):
	return all(sh.hamming(b, o) >= d for o in others)


@pytest.fixture(scope="module")
def cb_design(synth_genome, tmp_path_factory):
	"""Whitelist built around W1:
	W2 = W1 with bases 3 and 10 changed (2 mismatches), W4 = W1 with base 14
	changed (1 mismatch), W7 = W1 with bases 6 and 12 changed; W3, W5, W6
	and ten fillers are >= 3 mismatches from everything else. W1 and W2 get
	2 exactly matching reads each, W5 one, W3, W4, W6, W7 none.

	Reads (sequenced CB; all cDNA reads lie in G07's exon, every UMI is new):
	X  = W1 with base 3 set to W2's: 1MM from W1 (base 3) and W2 (base 10);
	     quality Phred 2 at base 3 and 40 at base 10
	Y  = the same sequence with Phred 37 everywhere
	Q  = W1 with base 6 set to W7's: 1MM from W1 and from W7
	Y5 = W5 with one change, 1MM from W5 only
	Y6 = W6 with one change, 1MM from W6 only
	Z1 = W5 with N at base 2 (one whitelist match); Z2 = W1 with N at base
	14 (matches W1 and W4); Z3 = W5 with N at bases 2 and 9
	FAR = >= 2 mismatches from every whitelist barcode
	NU = W5 with an N in the UMI; HU = W5 with a homopolymer UMI
	"""
	sim = hd.module_sim(synth_genome, __name__, "cb_design", const_qual=HIGH)
	w1 = sim.random_seq(16)
	w2 = hd.substitute(hd.substitute(w1, 3), 10)
	w4 = hd.substitute(w1, 14)
	w7 = hd.substitute(hd.substitute(w1, 6), 12)
	wl = OrderedDict([("W1", w1), ("W2", w2), ("W4", w4), ("W7", w7)])
	for name in ["W3", "W5", "W6"] + ["F%d" % i for i in range(10)]:
		b = sim.random_seq(16)
		while not _far_from(b, wl.values()):
			b = sim.random_seq(16)
		wl[name] = b
	path = hd.write_whitelist(tmp_path_factory.mktemp("wl_d") / "wl_design.txt", list(wl.values()))
	x = w1[:3] + w2[3] + w1[4:]
	q = w1[:6] + w7[6] + w1[7:]
	far = sim.random_seq(16)
	while not _far_from(far, wl.values(), 2):
		far = sim.random_seq(16)
	w5, w6 = wl["W5"], wl["W6"]
	qual_x = HIGH * 3 + LOW + HIGH * 6 + TOP + HIGH * 5
	specs = [("E1", w1, None), ("E1", w1, None), ("E2", w2, None), ("E2", w2, None), ("E5", w5, None),
		("X", x, qual_x), ("Y", x, None), ("Q", q, None), ("Y5", hd.substitute(w5, 8), None),
		("Y6", hd.substitute(w6, 11), None), ("Z1", w5[:2] + "N" + w5[3:], None),
		("Z2", w1[:14] + "N" + w1[15:], None), ("Z3", w5[:2] + "N" + w5[3:9] + "N" + w5[10:], None),
		("FAR", far, None), ("NU", w5, None), ("HU", w5, None)]
	pool = hd.UmiPool(sim)
	pairs = []
	for kind, cb_read, qual in specs:
		umi = pool.take()
		umi_read = {"NU": umi[:5] + "N" + umi[6:], "HU": "A" * 12}.get(kind, umi)
		cdna = sim.genomic_read("chr1", 30101 + 7 * len(pairs), "-", 90, name=sim.name("cell"))
		bq = (qual or HIGH * 16) + HIGH * 12
		p = hd.tenx_pair(sim, cdna, None, umi, cb_read=cb_read, umi_read=umi_read, bc_qual=bq, gene="G07")
		p.extra["kind"] = kind
		pairs.append(p)
	return wl, path, ReadSet(pairs, kind="solo")


# Whitelist barcode each read is counted for, per --soloCBmatchWLtype
# (parametersDefault soloCBmatchWLtype; None = not counted).
# 1MM_multi: posterior p(Wi) ~ (exact reads of Wi) * 10^(-Q/10) at the
# mismatch, only whitelist CBs with exact reads take part: X -> W1
# (0.9992), Y -> 0.5 (not called), Q -> W1 (W7 has no exact reads).
# Pseudocounts add 1 to every whitelist CB: Q -> W1 3 : W7 1 (0.75, not
# called), Y6 counted although W6 has no exact reads. Without Nbase, an N
# that matches two whitelist CBs (Z2) is rejected; with it the posterior is
# W1 3 : W4 1 (not called).
ASSIGN = {
	"E1": "W1", "E2": "W2", "E5": "W5",
	"X": {"Exact": None, "1MM": None, "default": "W1"},
	"Y": None,
	"Q": {"1MM_multi": "W1", "default": None},
	"Y5": {"Exact": None, "default": "W5"},
	"Y6": {"1MM_multi_pseudocounts": "W6", "1MM_multi_Nbase_pseudocounts": "W6", "default": None},
	"Z1": {"Exact": None, "default": "W5"},
	"Z2": None, "Z3": None, "FAR": None, "NU": None, "HU": None,
}

# Barcodes.stats per mode: categories that are not zero.
BARCODE_STATS = {
	"Exact": {"yesWLmatchExact": 5, "noNoWLmatch": 8, "noNinCB": 1, "noNinUMI": 1, "noUMIhomopolymer": 1},
	"1MM": {"yesWLmatchExact": 5, "yesOneWLmatchWithMM": 3, "noTooManyWLmatches": 4, "noNoWLmatch": 1,
		"noNinCB": 1, "noNinUMI": 1, "noUMIhomopolymer": 1},
	"1MM_multi": {"yesWLmatchExact": 5, "yesOneWLmatchWithMM": 3, "yesMultWLmatchWithMM": 3,
		"noTooManyWLmatches": 1, "noNoWLmatch": 1, "noNinCB": 1, "noNinUMI": 1, "noUMIhomopolymer": 1},
	"1MM_multi_pseudocounts": {"yesWLmatchExact": 5, "yesOneWLmatchWithMM": 3, "yesMultWLmatchWithMM": 3,
		"noTooManyWLmatches": 1, "noNoWLmatch": 1, "noNinCB": 1, "noNinUMI": 1, "noUMIhomopolymer": 1},
	"1MM_multi_Nbase_pseudocounts": {"yesWLmatchExact": 5, "yesOneWLmatchWithMM": 3, "yesMultWLmatchWithMM": 4,
		"noNoWLmatch": 1, "noNinCB": 1, "noNinUMI": 1, "noUMIhomopolymer": 1},
}

# Features.stats: reads whose CB matched with one mismatch to a whitelist CB
# without exact reads (noMMtoWLwithoutExact) or whose posterior was too low
# (noTooManyWLmatches).
FEATURE_REJECTS = {
	"Exact": (0, 0), "1MM": (1, 0), "1MM_multi": (1, 1),
	"1MM_multi_pseudocounts": (0, 2), "1MM_multi_Nbase_pseudocounts": (0, 3),
}


def _assigned(kind, mode):
	a = ASSIGN[kind]
	if isinstance(a, dict):
		return a.get(mode, a["default"])
	return a


@pytest.mark.parametrize("mode", MODES)
def test_cb_whitelist_matching(genome_sjdb, cb_design, mode):
	"""Per --soloCBmatchWLtype, each read is counted for the whitelist CB in
	ASSIGN (or not at all): raw Gene counts per cell, Barcodes.stats,
	Features.stats (yesWLmatch, yessubWLmatchExact, noMMtoWLwithoutExact,
	noTooManyWLmatches) and Summary.csv 'Reads With Valid Barcodes' =
	1 - (invalid barcode reads + Features.stats rejects) / reads."""
	wl, path, reads = cb_design
	run = run_star(genome_sjdb, reads, hd.solo_args(path, "--soloCBmatchWLtype", mode), name="cb_" + mode)
	expected = {}
	for p in reads:
		w = _assigned(p.extra["kind"], mode)
		if w:
			expected[("G07", wl[w])] = expected.get(("G07", wl[w]), 0) + 1
	gene = run.solo("Gene")
	assert gene.counts == expected
	assert gene.barcodes == sorted(wl.values())
	bstats = run.solo_barcodes_stats()
	want = dict.fromkeys(bstats, 0)
	want.update(BARCODE_STATS[mode])
	assert dict(bstats) == want
	no_exact, too_many = FEATURE_REJECTS[mode]
	fs = run.solo_features_stats("Gene")
	assert fs["yesWLmatch"] == sum(expected.values())
	assert fs["yessubWLmatchExact"] == 5
	assert (fs["noMMtoWLwithoutExact"], fs["noTooManyWLmatches"]) == (no_exact, too_many)
	invalid = sum(v for k, v in BARCODE_STATS[mode].items() if k.startswith("no")) + no_exact + too_many
	summary = run.solo_summary("Gene")
	assert summary["Number of Reads"] == len(reads)
	assert summary["Reads With Valid Barcodes"] == pytest.approx(1 - invalid / float(len(reads)), abs=1e-6)
	assert summary["Reads Mapped to Gene: Unique Gene"] == pytest.approx(sum(expected.values()) / float(len(reads)),
		abs=1e-6)


def test_no_whitelist_counts_every_barcode(genome_sjdb, cb_design):
	"""--soloCBwhitelist None accepts every CB as sequenced (no
	correction): each read with an N-free CB and a valid UMI is counted for
	its own CB, and barcodes.tsv lists those CBs sorted."""
	_, _, reads = cb_design
	run = run_star(genome_sjdb, reads, hd.solo_args("None"))
	expected = {}
	for p in reads:
		if "N" in p.extra["cb_read"] or p.extra["kind"] in ("NU", "HU"):
			continue
		expected[("G07", p.extra["cb_read"])] = expected.get(("G07", p.extra["cb_read"]), 0) + 1
	gene = run.solo("Gene")
	assert gene.counts == expected
	assert gene.barcodes == sorted({cb for _, cb in expected})


@pytest.fixture(scope="module")
def simple_reads(synth_genome, tenx_whitelist):
	sim = hd.module_sim(synth_genome, __name__, "simple_reads")
	wl = tenx_whitelist[0]
	cells = {wl[1]: {"G01": 2, "G08": 1}, wl[5]: {"G03": 3}, wl[9]: {"G11": 1, "G02": 2}}
	return sim.tenx_reads(wl, cells, reads_per_umi=(1, 2), cb_errors={wl[5]: 1})


def _with_barcode_reads(reads, transform, quals=None):
	pairs = []
	for p in reads:
		bc = transform(p.r2.seq, p)
		r2 = SimRead(p.name, bc, (quals or (lambda q, p: "F" * len(bc)))(p.r2.qual, p), kind="barcode")
		pairs.append(SimPair(p.name, p.r1, r2, None, "solo", extra=p.extra))
	return ReadSet(pairs, kind="solo")


def test_barcode_read_length_check(genome_sjdb, simple_reads, tenx_whitelist):
	"""--soloBarcodeReadLength 1 (default) requires the barcode read to be
	CB+UMI = 28 bases: 30-base barcode reads stop STAR with an input error.
	With --soloBarcodeReadLength 0 (no check) or 30 the extra bases are
	ignored and the matrix equals the truth."""
	longer = _with_barcode_reads(simple_reads, lambda s, p: s + "GG")
	bad = run_star(genome_sjdb, longer, hd.solo_args(tenx_whitelist[1]), check=False)
	assert bad.returncode != 0
	assert "not equal to expected 28" in bad.stderr + bad.log_out()
	for n in ("0", "30"):
		run = run_star(genome_sjdb, longer, hd.solo_args(tenx_whitelist[1], "--soloBarcodeReadLength", n))
		assert run.solo("Gene").counts == sh.expected_solo_counts(simple_reads)


def test_barcode_read_must_be_last_file(genome_sjdb, simple_reads, tenx_whitelist, tmp_path):
	"""docs/STARsolo.md: the 1st file has to be the cDNA read and the 2nd
	the barcode read. Given in the other order, the 90-base cDNA reads are
	taken as barcode reads and fail the 28-base length check."""
	files = simple_reads.write(tmp_path, prefix="tenx")
	run = run_star(genome_sjdb, [str(files[1]), str(files[0])], hd.solo_args(tenx_whitelist[1]), check=False)
	assert run.returncode != 0
	assert "not equal to expected 28" in run.stderr + run.log_out()


def test_gzipped_cdna_and_barcode_files(genome_sjdb, simple_reads, tenx_whitelist):
	"""cDNA and barcode FASTQ files both gzipped and read through
	--readFilesCommand gzip -cd give the truth matrix (and the same one as
	plain files); Barcodes.stats counts the one corrected CB."""
	run = run_star(genome_sjdb, simple_reads, hd.solo_args(tenx_whitelist[1]), gz=True)
	assert "--readFilesCommand" in run.argv
	assert run.solo("Gene").counts == sh.expected_solo_counts(simple_reads)
	bstats = run.solo_barcodes_stats()
	assert bstats["yesOneWLmatchWithMM"] == 1
	assert bstats["yesWLmatchExact"] == len(simple_reads) - 1


def test_barcode_geometry_umi_before_cb(genome_sjdb, simple_reads, tenx_whitelist):
	"""--soloUMIstart/--soloCBstart place UMI and CB anywhere in the
	barcode read: with 4 leading bases, the UMI at 5-16 and the CB at 17-32
	(--soloBarcodeReadLength 32) the matrix equals the truth. CR/UR in the
	BAM are the sequenced CB and UMI."""
	moved = _with_barcode_reads(simple_reads, lambda s, p: "ACGT" + s[16:] + s[:16])
	run = run_star(genome_sjdb, moved, ["--soloType", "CB_UMI_Simple", "--soloCBwhitelist", str(tenx_whitelist[1]),
		"--soloUMIstart", "5", "--soloUMIlen", "12", "--soloCBstart", "17", "--soloCBlen", "16",
		"--soloBarcodeReadLength", "32", "--outSAMtype", "BAM", "Unsorted", "--outSAMattributes", "NH", "CR", "UR"])
	assert run.solo("Gene").counts == sh.expected_solo_counts(simple_reads)
	for rec in run.records():
		p = simple_reads.by_name[rec.qname]
		assert (rec.tag("CR"), rec.tag("UR")) == (p.extra["cb_read"], p.extra["umi"])
