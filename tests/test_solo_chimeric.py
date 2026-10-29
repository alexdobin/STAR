"""Area D: STARsolo CB_UMI_Simple (10x 3' v3 geometry) and chimeric output."""

import pytest

import star_helpers as sh
from star_helpers import ReadSet, run_star


def _solo_args(whitelist_path, *extra):
	return ["--soloType", "CB_UMI_Simple", "--soloCBwhitelist", str(whitelist_path), "--soloUMIlen", "12"] + list(extra)


@pytest.fixture(scope="module")
def solo_cells(tenx_whitelist):
	wl = tenx_whitelist[0]
	return {
		wl[0]: {"G01": 3, "G02": 2, "G07": 4},
		wl[3]: {"G03": 1, "G08": 5, "G11": 2},
		wl[7]: {"G10": 2, "G01": 1},
		wl[12]: {"G04": 3},
		wl[20]: {"G09": 1, "G02": 1, "G08": 1},
	}


def test_solo_gene_matrix_matches_truth(genome_sjdb, tenx_whitelist, solo_cells, sim):
	"""CB_UMI_Simple with a 16 bp CB + 12 bp UMI barcode read and 1-3 reads
	per UMI: the raw Gene matrix holds, for every (gene, cell), the number of
	distinct UMIs simulated there (exact CBs plus CBs with one substitution,
	which 1MM_multi corrects). features.tsv lists the annotated genes as
	(id, name, Gene Expression), barcodes.tsv the sorted whitelist. Every read
	is exonic and sense, so GeneFull equals Gene. Features.stats: every read
	matches the whitelist and one gene; yesUMIs is the UMI total;
	Barcodes.stats splits exact and one-mismatch CB matches."""
	wl, wl_path = tenx_whitelist
	reads = sim.tenx_reads(wl, solo_cells, reads_per_umi=(1, 3), cb_errors={wl[0]: 2, wl[3]: 1})
	run = run_star(genome_sjdb, reads, _solo_args(wl_path, "--soloFeatures", "Gene", "GeneFull"))
	expected = sh.expected_solo_counts(reads)
	gene = run.solo("Gene")
	assert gene.counts == expected
	assert gene.features == [(g.gene_id, g.name, "Gene Expression") for g in sim.g.annotated_genes()]
	assert gene.barcodes == sorted(wl)
	assert gene.shape == (len(gene.features), len(wl))
	assert run.solo("GeneFull").counts == expected

	n_err = sum(p.extra["cb_read"] != p.extra["cb"] for p in reads)
	assert n_err == 3
	fstats = run.solo_features_stats("Gene")
	assert fstats["yesWLmatch"] == len(reads)
	assert fstats["yessubWLmatch_UniqueFeature"] == len(reads)
	assert fstats["yesUMIs"] == sum(expected.values())
	assert fstats["yesCellBarcodes"] == len(solo_cells)
	bstats = run.solo_barcodes_stats()
	assert bstats["yesWLmatchExact"] == len(reads) - n_err
	assert bstats["yesOneWLmatchWithMM"] == n_err
	assert run.solo_summary("Gene")["Number of Reads"] == len(reads)


def test_solo_exact_cb_matching_drops_mismatched_barcodes(genome_sjdb, tenx_whitelist, solo_cells, sim):
	"""With --soloCBmatchWLtype Exact the reads whose CB carries one
	substitution are not counted: the matrix equals the truth computed
	without them, and Barcodes.stats reports them as noNoWLmatch."""
	wl, wl_path = tenx_whitelist
	reads = sim.tenx_reads(wl, solo_cells, reads_per_umi=1, cb_errors={wl[7]: 3})
	run = run_star(genome_sjdb, reads, _solo_args(wl_path, "--soloCBmatchWLtype", "Exact"))
	exact = ReadSet([p for p in reads if p.extra["cb_read"] == p.extra["cb"]], kind="solo")
	assert len(exact) == len(reads) - 3
	assert run.solo("Gene").counts == sh.expected_solo_counts(exact)
	assert run.solo_barcodes_stats()["noNoWLmatch"] == 3


def test_solo_matrix_same_with_four_threads(genome_sjdb, tenx_whitelist, solo_cells, sim):
	"""The raw Gene matrix and Features.stats do not depend on --runThreadN:
	4 threads on input chunks of 8698 bytes per mate (--limitIObufferSize
	220000 50000000), so more than one thread maps reads."""
	wl, wl_path = tenx_whitelist
	reads = sim.tenx_reads(wl, solo_cells, reads_per_umi=(1, 4), cb_errors={wl[12]: 2})
	one = run_star(genome_sjdb, reads, _solo_args(wl_path), threads=1)
	four = run_star(genome_sjdb, reads, _solo_args(wl_path, "--limitIObufferSize", "220000", "50000000"), threads=4)
	assert one.solo().counts == four.solo().counts == sh.expected_solo_counts(reads)
	assert one.solo_features_stats() == four.solo_features_stats()


def _chimeric_reads(sim, n_per_combo=3):
	"""50+50 base reads joining unique, gene-free loci on two contigs for all
	four strand combinations (ReadSimulator.clean_chimeric_read)."""
	return ReadSet([sim.clean_chimeric_read(strands) for strands in (("+", "+"), ("+", "-"), ("-", "+"), ("-", "-"))
		for _ in range(n_per_combo)])


def test_chimeric_junction_columns(genome_plain, sim):
	"""With --chimSegmentMin 20, each 50+50 chimeric read gets one line in
	Chimeric.out.junction whose first 14 columns equal the values worked out
	from the two pieces (manual 6.4): donor breakpoint = first base after
	piece A on its strand, acceptor breakpoint = base before piece B on its
	strand, junction type 0, repeats 0 (pieces chosen so the breakpoint cannot
	move), and SAM-like starts/CIGARs on the +
	strand. Log.final.out counts every read as chimeric."""
	reads = _chimeric_reads(sim)
	run = run_star(genome_plain, reads, "--chimSegmentMin 20")
	rows = run.chimeric()
	assert sorted(r["read_name"] for r in rows) == sorted(reads.by_name)
	for row in rows:
		truth = reads.by_name[row["read_name"]]
		assert [row[k] for k in sh.CHIMERIC_FIELDS] == list(truth.extra["junction"].values())
	assert run.log_final()["Number of chimeric reads"] == len(reads)


def test_chimeric_within_bam(genome_plain, sim):
	"""--chimOutType WithinBAM writes each chimeric read as two records in the
	BAM, one of them supplementary (0x800, hard-clipped), at the two pieces'
	positions and strands; each record's SA tag names the other piece."""
	reads = _chimeric_reads(sim, n_per_combo=2)
	run = run_star(genome_plain, reads, "--chimSegmentMin 20 --chimOutType WithinBAM --outSAMtype BAM Unsorted")
	by_name = run.alignments().by_name()
	assert sorted(by_name) == sorted(reads.by_name)
	for name, recs in by_name.items():
		(ca, pa, sa, _), (cb, pb, sb, _) = reads.by_name[name].extra["pieces"]
		assert len(recs) == 2
		assert sum(r.is_supplementary for r in recs) == 1
		assert "H" in next(r for r in recs if r.is_supplementary).cigar
		assert sorted((r.rname, r.pos, r.strand) for r in recs) == sorted([(ca, pa, sa), (cb, pb, sb)])
		for r in recs:
			other = next(x for x in recs if x is not r)
			sa_fields = r.tag("SA").rstrip(";").split(",")
			assert (sa_fields[0], int(sa_fields[1]), sa_fields[2], sa_fields[3]) == (
				other.rname, other.pos, other.strand, other.cigar)
