"""Area C: splice junctions inserted at the mapping step (--sjdbFileChrStartEnd,
--sjdbGTFfile with its tag options, --sjdbOverhang, --sjdbInsertSave,
--sjdbScore) and 2-pass mapping (--twopassMode Basic, --twopass1readsN).

Reads with a 4-base overhang are the probe for "junction inserted": below
--alignSJoverhangMin 5 they cannot be spliced across an unannotated
junction, but --alignSJDBoverhangMin 3 lets them cross an inserted one
(parametersDefault).
"""

import re

import pytest

import helpers_c as hc
import star_helpers as sh
from star_helpers import ReadSet, record_multiset, record_truth, run_star


def _short(sim, contig, intron, strand="+", left=96, name=None):
	"""Read crossing `intron` with a 4-base overhang (left=96) or a 96-base
	right block (left=4)."""
	return hc.junction_read(sim, contig, intron, left, 100 - left, strand, name=name)


def _assert_truth(run, reads):
	recs = run.alignments().by_name()
	for r in reads:
		assert [record_truth(x) for x in recs.get(r.name, [])] == [r.truth()], r.name


def _unspliced(run, names):
	recs = run.alignments().by_name()
	for n in names:
		assert recs[n] and all("N" not in x.cigar for x in recs[n]), n


def test_sjdb_file_chr_start_end_on_the_fly(genome_plain, sim, tmp_path):
	"""--sjdbFileChrStartEnd at the mapping step against an index built
	without junctions (manual 2.2.4, 3.3.1, 9.1). File 1 lists T10skip
	(given twice), T02 intron 1 with strand '-', and the T11 GT/AT intron
	with strand '.'; file 2 is the SJ.out.tab of a first run over T03 (GC/AG)
	and T12nc (non-canonical), whose numeric strand column is accepted.
	Both files are read ("multiple files ... concatenated") and the
	duplicate is collapsed, so 5 junctions are inserted with the default
	--sjdbOverhang 100 (_STARgenome/sjdbInfo.txt header "5<TAB>100"). One
	4-base-overhang read per inserted junction gets its true spliced
	alignment and an annotated SJ.out.tab row (unique 1, overhang 4: the
	overhang filter does not apply to annotated junctions); the same kind of
	read over the T05 intron, not listed, stays unspliced."""
	first = [hc.junction_read(sim, "chr1", (11401, 11800), 50, 50)] + \
		[hc.junction_read(sim, "chr3", (15201, 15800), l, 100 - l, s) for l, s in ((50, "+"), (40, "-"), (60, "+"))]
	run1 = run_star(genome_plain, ReadSet(first))
	assert [(r.contig, r.start, r.end, r.strand, r.annotated) for r in run1.sj()] == [
		("chr1", 11401, 11800, 1, 0), ("chr3", 15201, 15800, 0, 0)]
	listed = tmp_path / "junctions.tab"
	listed.write_text("chr2\t30401\t31800\t+\nchr2\t30401\t31800\t+\nchr1\t6401\t7000\t-\nchr3\t12501\t13000\t.\n")
	probes = [_short(sim, "chr2", (30401, 31800), "+", name="p_t10skip"),
		_short(sim, "chr1", (6401, 7000), "-", left=4, name="p_t02"),
		_short(sim, "chr3", (12501, 13000), "-", name="p_t11"),
		_short(sim, "chr1", (11401, 11800), "+", left=4, name="p_t03"),
		_short(sim, "chr3", (15201, 15800), "+", name="p_t12nc")]
	control = _short(sim, "chr1", (20601, 21000), "+", name="p_t05")
	run = run_star(genome_plain, ReadSet(probes + [control]),
		["--sjdbFileChrStartEnd", str(listed), str(run1.path("SJ.out.tab"))])
	assert run.path("_STARgenome/sjdbInfo.txt").read_text().splitlines()[0] == "5\t100"
	_assert_truth(run, probes)
	_unspliced(run, ["p_t05"])
	expected = hc.expected_sj(hc.truth_alignments(probes), sim.g,
		{(r.contig,) + tuple(r.introns[0]) for r in probes})
	assert [(r.annotated, r.unique, r.max_overhang) for r in expected] == [(1, 1, 4)] * 5
	assert run.sj() == expected


def test_sjdb_overhang_and_insert_save_all(genome_plain, genome_sjdb, sim):
	"""--sjdbGTFfile with --sjdbOverhang 75 at the mapping step inserts the
	16 GTF introns (_STARgenome/sjdbInfo.txt header "16<TAB>75";
	sjdbList.out.tab lists each intron with its transcript's strand) and
	aligns the reads as the index built with the GTF and --sjdbOverhang 99
	does, 4-base-overhang reads included. --sjdbInsertSave All writes a
	complete genome into _STARgenome (manual 3.3.1): mapping against it with
	no sjdb options gives the same alignments and SJ.out.tab. The default
	Basic saves only the junction and transcript files, not Genome, SA and
	SAindex."""
	short = [_short(sim, "chr1", (1301, 2300), "+"), _short(sim, "chr1", (7201, 8000), "-", left=4),
		_short(sim, "chr2", (2301, 3000), "+", left=4), _short(sim, "chr3", (12501, 13000), "-")]
	reads = sim.transcript_reads(60) + ReadSet(short)
	gtf = ["--sjdbGTFfile", str(sim.g.gtf_path), "--sjdbOverhang", "75"]
	fly = run_star(genome_plain, reads, gtf + ["--sjdbInsertSave", "All"])
	saved = fly.path("_STARgenome")
	assert saved.joinpath("sjdbInfo.txt").read_text().splitlines()[0] == "16\t75"
	listed = {tuple(l.split("\t")) for l in saved.joinpath("sjdbList.out.tab").read_text().splitlines()}
	expected_list = {(j.contig, str(j.start), str(j.end), j.strand) for j in sim.g.junctions()}
	assert len(expected_list) == 16 and listed == expected_list
	_assert_truth(fly, short)
	pre = run_star(genome_sjdb, reads)
	assert record_multiset(fly.records()) == record_multiset(pre.records())
	assert fly.sj() == pre.sj()
	for f in ("Genome", "SA", "SAindex", "genomeParameters.txt", "chrName.txt"):
		assert saved.joinpath(f).exists(), f
	reuse = run_star(saved, reads)
	assert record_multiset(reuse.records()) == record_multiset(fly.records())
	assert reuse.sj() == fly.sj()
	basic = run_star(genome_plain, ReadSet(short), gtf)
	names = {p.name for p in basic.path("_STARgenome").iterdir()}
	assert "sjdbList.out.tab" in names and not names & {"Genome", "SA", "SAindex"}
	_assert_truth(basic, short)


@pytest.mark.parametrize("case", [("gtf", 15), ("gtf", 16), ("twopass", 16), ("twopass", 17)],
	ids=["gtf15", "gtf16", "twopass16", "twopass17"])
def test_limit_sjdb_insert_nsj(genome_plain, genome_sjdb, sim, case):
	"""--limitSjdbInsertNsj is the most junctions inserted on the fly,
	including those found in the 1st pass of a 2-pass run
	(parametersDefault). The GTF has 16 introns: inserting it at the mapping
	step works with a limit of 16 (the 4-base-overhang read splices) and
	stops with a LIMIT error at 15. A 2-pass run on the GTF index whose 1st
	pass finds the unannotated T10skip junction inserts 17: it works with 17
	and stops at 16."""
	mode, limit = case
	if mode == "gtf":
		probe = _short(sim, "chr1", (1301, 2300), "+")
		reads = [probe]
		genome, args, needed = genome_plain, ["--sjdbGTFfile", str(sim.g.gtf_path)], 16
	else:
		probe = _short(sim, "chr2", (30401, 31800), "+")
		reads = [probe] + [sim.transcript_read("T10skip", o) for o in (330, 360)]
		genome, args, needed = genome_sjdb, ["--twopassMode", "Basic"], 17
	run = run_star(genome, ReadSet(reads), args + ["--limitSjdbInsertNsj", str(limit)], check=False)
	if limit >= needed:
		assert run.returncode == 0, run.stderr
		_assert_truth(run, reads)
	else:
		assert run.returncode != 0
		assert "=%d is larger than the limitSjdbInsertNsj=%d" % (needed, limit) in run.stderr


def _custom_gtf(synth, path):
	"""The synthetic GTF rewritten for the tag options: contigs without the
	'chr' prefix, exon lines typed 'myexon', extra attributes tx_alias
	"TA_<transcript>" and gene_alias "GA_<gene>". Two decoy lines typed
	'exon' describe a fake transcript DECOYT with intron chr1:44101-44600."""
	lines = []
	for line in synth.gtf_text().splitlines():
		f = line.split("\t")
		if f[2] != "exon":
			continue
		gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
		tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
		f[0] = f[0][len("chr"):]
		f[2] = "myexon"
		f[8] = 'gene_alias "GA_%s"; tx_alias "TA_%s"; %s' % (gid, tid, f[8])
		lines.append(f)
	decoy = 'gene_id "DECOY"; transcript_id "DECOYT"; gene_alias "GA_DECOY"; tx_alias "TA_DECOY";'
	lines.append(["1", "x", "exon", "44001", "44100", ".", "+", ".", decoy])
	lines.append(["1", "x", "exon", "44601", "44800", ".", "+", ".", decoy])
	return hc.write_gtf(path, lines)


def test_sjdb_gtf_chr_prefix_feature_and_tags(genome_plain, sim, tmp_path):
	"""A GTF whose contigs lack 'chr', whose exons are typed 'myexon' and
	whose parent IDs sit in tx_alias/gene_alias is read at the mapping step
	with --sjdbGTFchrPrefix chr --sjdbGTFfeatureExon myexon
	--sjdbGTFtagExonParentTranscript tx_alias --sjdbGTFtagExonParentGene
	gene_alias (parametersDefault, manual 2.2.3): the 4-base-overhang reads
	over GTF introns splice, the read over the decoy 'exon' intron does not,
	ReadsPerGene.out.tab rows are GA_<gene> in GTF order without the decoy
	gene, and Aligned.toTranscriptome.out.bam references are TA_<transcript>.
	Gene counts by hand: 2 sense G01 reads, 1 sense G02 read, 2 antisense
	G07 reads, 1 intergenic decoy read. Without the options STAR finds no
	usable exon line and exits with an error."""
	gtf = _custom_gtf(sim.g, tmp_path / "custom.gtf")
	g01 = [_short(sim, "chr1", (1301, 2300), "+", name="g01a"), _short(sim, "chr1", (1301, 3450), "+", left=4, name="g01b")]
	g02 = [_short(sim, "chr1", (7201, 8000), "-", left=4, name="g02")]
	g07 = [sim.transcript_read("T07", o, sense=False) for o in (100, 600)]
	decoy = _short(sim, "chr1", (44101, 44600), "+", name="decoy")
	reads = ReadSet(g01 + g02 + g07 + [decoy])
	args = ["--sjdbGTFfile", str(gtf), "--sjdbGTFchrPrefix", "chr", "--sjdbGTFfeatureExon", "myexon",
		"--sjdbGTFtagExonParentTranscript", "tx_alias", "--sjdbGTFtagExonParentGene", "gene_alias",
		"--quantMode", "GeneCounts", "TranscriptomeSAM"]
	run = run_star(genome_plain, reads, args)
	_assert_truth(run, g01 + g02 + g07)
	_unspliced(run, ["decoy"])
	genes = ["GA_" + g.gene_id for g in sim.g.annotated_genes()]
	counts = run.gene_counts()
	assert list(counts) == ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"] + genes
	expected = {"N_noFeature": (1, 3, 4), "GA_G01": (2, 2, 0), "GA_G02": (1, 1, 0), "GA_G07": (2, 0, 2)}
	for k, v in counts.items():
		assert v == expected.get(k, (0, 0, 0)), k
	tx = run.transcriptome()
	assert dict(tx.references) == {"TA_" + t.tx_id: t.length for t in sim.g.transcripts.values() if t.annotated}
	by_name = tx.by_name()
	# g02 covers 7197-7200 and 8001-8096; on the - strand T02 its first base is 8096
	assert [(r.rname, r.pos, r.cigar) for r in by_name["g02"]] == [("TA_T02", hc.tx_offset_of(sim.g, "T02", 8096) + 1, "100M")]
	plain = run_star(genome_plain, reads, ["--sjdbGTFfile", str(gtf)], check=False)
	assert plain.returncode != 0
	assert "no valid exon lines in the GTF file" in plain.stderr


@pytest.mark.parametrize("score", [None, 0, 1, 3, -3], ids=["default2", "0", "1", "3", "minus3"])
def test_sjdb_score_breaks_tie_between_copies(genome_plain, sim, tmp_path, score):
	"""An intron inside the duplicated dup2 segment is listed in
	--sjdbFileChrStartEnd for the chr1 copy only. Reads crossing it align
	equally well to both copies except that the chr1 alignment crosses a
	database junction and gains --sjdbScore (default 2, parametersDefault).
	With --outFilterMultimapScoreRange 1: a gain of 2 or 3 makes chr1 the
	only alignment; 0 or 1 keeps both (NH 2); -3 leaves only chr2. The read
	with a 4-base overhang at the chr1 copy can only use the database
	junction (no novel junction with an overhang < 5), so it is unique at
	chr1 whenever sjdbScore >= 0. SJ.out.tab: chr1 row annotated 1, chr2
	row annotated 0, counts following the alignments."""
	a, b = hc.segment_intron(sim.g, "dup2")
	(c1, s1, e1), (c2, s2, e2) = hc.segment_intron_loci(sim.g, "dup2", a, b)
	listed = tmp_path / "dup2.tab"
	listed.write_text("%s\t%d\t%d\t+\n" % (c1, s1, e1))
	dup = [hc.segment_junction_read(sim, "dup2", a, b, 50, 50, copy=i, strand="+-"[i]) for i in (0, 1)]
	short4 = hc.junction_read(sim, c1, (s1, e1), 96, 4, "+", name="short4")
	args = ["--sjdbFileChrStartEnd", str(listed)]
	if score is not None:
		args += ["--sjdbScore", str(score)]
	run = run_star(genome_plain, ReadSet(dup + [short4]), args)
	gain = 2 if score is None else score
	recs = run.alignments().by_name()
	sj = {(r.contig, r.start, r.end): r for r in run.sj()}
	for r in dup:
		loci = r.extra["loci"]
		got = sorted(record_truth(x) for x in recs[r.name])
		if gain > 1:
			assert got == [loci[0]]
		elif gain < -1:
			assert got == [loci[1]]
		else:
			assert got == sorted(loci)
			assert all(x.tag("NH") == 2 for x in recs[r.name])
			if gain == 1:
				assert [x.rname for x in recs[r.name] if not x.is_secondary] == [c1]
	if gain >= 0:
		assert [record_truth(x) for x in recs["short4"]] == [short4.truth()]
	if gain > 1:
		assert sj[(c1, s1, e1)][3:] == (1, 1, 1, 3, 0, 50)
		assert (c2, s2, e2) not in sj
	elif gain >= 0:
		assert sj[(c1, s1, e1)][3:] == (1, 1, 1, 1, 2, 50)
		assert sj[(c2, s2, e2)][3:] == (2, 2, 0, 0, 2, 50)
	else:
		assert sj[(c2, s2, e2)][3:] == (2, 2, 0, 2, 0, 50)


@pytest.mark.parametrize("n1", [8, 9])
def test_twopass1readsN(genome_plain, sim, n1):
	"""--twopass1readsN N maps only the first N reads in the 1st pass
	(manual 9.2). Input order: 3 reads over the unannotated T10skip intron
	with a 4-base overhang, 5 unspliced reads, then 8 reads over T10skip
	with long overhangs. With N=8 the 1st pass sees no long read, finds no
	junction (_STARpass1/SJ.out.tab empty, _STARpass1/Log.final.out 8 input
	reads), so in the 2nd pass the short reads stay unspliced and the
	junction is unannotated with 8 unique reads. With N=9 the 1st pass maps
	one long read, its junction is inserted, the short reads splice, and
	the junction is reported annotated with 11 unique reads. The final
	Log.final.out counts the 16 reads of the 2nd pass only."""
	short = [_short(sim, "chr2", (30401, 31800), "+", name="short1"), _short(sim, "chr2", (30401, 31800), "-", left=4, name="short2"),
		_short(sim, "chr2", (30401, 31800), "-", name="short3")]
	filler = list(sim.transcript_reads(5, spliced=False))
	long_reads = [sim.transcript_read("T10skip", o, sense=bool(i % 2)) for i, o in enumerate(range(325, 365, 5))]
	reads = ReadSet(short + filler + long_reads)
	run = run_star(genome_plain, reads, "--twopassMode Basic --twopass1readsN %d" % n1)
	assert run.path("_STARpass1/Log.final.out").exists()
	pass1 = sh.parse_log_final(run.path("_STARpass1/Log.final.out"))
	assert pass1["Number of input reads"] == n1
	pass1_sj = sh.parse_sj_out(run.path("_STARpass1/SJ.out.tab"))
	assert run.log_final()["Number of input reads"] == 16
	_assert_truth(run, filler + long_reads)
	sj = {(r.contig, r.start, r.end): r for r in run.sj()}
	if n1 == 8:
		assert pass1_sj == []
		_unspliced(run, [r.name for r in short])
		assert sj[("chr2", 30401, 31800)][3:8] == (1, 1, 0, 8, 0)
	else:
		assert [(r.contig, r.start, r.end, r.unique) for r in pass1_sj] == [("chr2", 30401, 31800, 1)]
		_assert_truth(run, short)
		assert sj[("chr2", 30401, 31800)][3:8] == (1, 1, 1, 11, 0)


def test_twopass_inserts_filtered_junctions_and_counts_once(genome_sjdb, sim):
	"""--twopassMode Basic with --quantMode GeneCounts on the GTF index.
	Pass 1 junctions are inserted only if they pass the SJ.out.tab filters:
	T10skip (3 long reads, GT/AG) is inserted, so its two 4-base-overhang
	reads splice in pass 2 and the junction is reported annotated with 5
	unique reads; T12nc (2 long reads, non-canonical, below
	--outSJfilterCountUniqueMin 3) is not, so its 4-base-overhang read stays
	unspliced and the junction is absent from SJ.out.tab. Quantification
	runs in the 2nd pass only (twoPassRunPass1 turns it off for pass 1), and
	Log.final.out describes the 2nd pass only: 12 input reads, 12 unique,
	splices 7 = 5 GT/AG (annotated, pass-1 junction) + 2 non-canonical.
	ReadsPerGene by hand: G10 5 sense reads, G01 4 sense reads, the 3 T12nc
	reads (G12 is not in the GTF) no feature."""
	t10_long = [sim.transcript_read("T10skip", o) for o in (330, 345, 360)]
	t10_short = [_short(sim, "chr2", (30401, 31800), "+", name="t10s1"), _short(sim, "chr2", (30401, 31800), "+", left=4, name="t10s2")]
	t12_long = [hc.junction_read(sim, "chr3", (15201, 15800), 50, 50, "+"), hc.junction_read(sim, "chr3", (15201, 15800), 40, 60, "+")]
	t12_short = _short(sim, "chr3", (15201, 15800), "+", name="t12s")
	g01 = [sim.transcript_read("T01a", o) for o in (10, 60, 110, 160)]
	reads = ReadSet(t10_long + t10_short + t12_long + [t12_short] + g01)
	run = run_star(genome_sjdb, reads, "--twopassMode Basic --quantMode GeneCounts")
	_assert_truth(run, t10_long + t10_short + t12_long + g01)
	_unspliced(run, ["t12s"])
	sj = {(r.contig, r.start, r.end): r for r in run.sj()}
	assert sj[("chr2", 30401, 31800)][3:8] == (1, 1, 1, 5, 0)
	assert ("chr3", 15201, 15800) not in sj
	assert len(sj) == 1
	log = run.log_final()
	assert log["Number of input reads"] == 12
	assert log["Uniquely mapped reads number"] == 12
	assert log["Number of splices: Total"] == 7
	assert log["Number of splices: Annotated (sjdb)"] == 5
	assert log["Number of splices: GT/AG"] == 5
	assert log["Number of splices: Non-canonical"] == 2
	counts = run.gene_counts()
	expected = {"N_noFeature": (3, 3, 12), "G01": (4, 4, 0), "G10": (5, 5, 0)}
	for k, v in counts.items():
		assert v == expected.get(k, (0, 0, 0)), k
