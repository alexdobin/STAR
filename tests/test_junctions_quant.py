"""Area C: SJ.out.tab and its filters, sjdb insertion on the fly, 2-pass
mode, --quantMode GeneCounts and TranscriptomeSAM, Log.final.out."""

from collections import Counter

import star_helpers as sh
from star_helpers import ReadSet, record_multiset, record_truth, run_star


def _expected_sj(reads, synth, annotated):
	"""SJRecord rows expected from unique, single-locus simulated reads."""
	counts = Counter()
	overhang = {}
	for r in reads:
		blocks = sh.cigar_blocks(r.pos, r.cigar)
		for i, (s, e) in enumerate(r.introns):
			key = (r.contig, s, e)
			counts[key] += 1
			left = blocks[i][1] - blocks[i][0] + 1
			right = blocks[i + 1][1] - blocks[i + 1][0] + 1
			overhang[key] = max(overhang.get(key, 0), min(left, right))
	rows = []
	for (contig, s, e), n in sorted(counts.items(), key=lambda kv: (list(synth.contigs).index(kv[0][0]), kv[0][1:])):
		motif = sh.motif_code(synth.seq(contig, s, s + 1), synth.seq(contig, e - 1, e))
		rows.append(sh.SJRecord(contig, s, e, sh.MOTIF_STRAND[motif], motif, int(annotated), n, 0, overhang[(contig, s, e)]))
	return rows


def test_sj_out_matches_truth_without_annotation(genome_plain, sim):
	"""With no annotation, SJ.out.tab lists exactly the introns crossed by the
	single-junction reads: strand and motif from the genome, annotated 0,
	unique count = reads crossing, multi 0, max overhang = the largest
	min(left block, right block). All pass the default filters (overhang >=
	25 > 12, one unique read is enough for canonical motifs)."""
	reads = ReadSet([r for r in sim.transcript_reads(400, spliced=True) if len(r.introns) == 1])
	run = run_star(genome_plain, reads)
	assert run.sj() == _expected_sj(reads, sim.g, annotated=False)


def test_sj_out_annotated_and_filtered(genome_sjdb, sim):
	"""With the GTF inserted, a crossed annotated junction (GC/AG, 1 read) is
	reported with annotated 1 and the unannotated T10skip junction (3 reads)
	with annotated 0; the non-canonical T12nc junction (2 reads) is dropped by
	the default outSJfilterCountUniqueMin of 3 for non-canonical motifs.
	--outSJfilterCountUniqueMin 4 4 4 4 with --outSJfilterCountTotalMin
	4 4 4 4 also drops T10skip but keeps the annotated junction (the filters
	do not apply to annotated junctions)."""
	t10skip = [sim.transcript_read("T10skip", o) for o in (330, 350, 360)]  # novel 30401-31800, 3 reads
	t12 = [sim.transcript_read("T12nc", o) for o in (130, 150)]  # novel non-canonical, 2 reads
	t03 = [sim.transcript_read("T03", 350)]  # annotated GC/AG, 1 read
	reads = ReadSet(t10skip + t12 + t03)
	sj = {(r.contig, r.start, r.end): r for r in run_star(genome_sjdb, reads).sj()}
	assert sj[("chr2", 30401, 31800)][3:8] == (1, 1, 0, 3, 0)
	assert sj[("chr1", 11401, 11800)][3:8] == (1, 3, 1, 1, 0)
	assert ("chr3", 15201, 15800) not in sj  # non-canonical needs 3 unique reads by default
	assert len(sj) == 2
	filt = "--outSJfilterCountUniqueMin 4 4 4 4 --outSJfilterCountTotalMin 4 4 4 4"
	sj = {(r.contig, r.start, r.end) for r in run_star(genome_sjdb, reads, filt).sj()}
	assert sj == {("chr1", 11401, 11800)}


def test_sjdb_inserted_on_the_fly_matches_prebuilt_index(genome_plain, genome_sjdb, sim):
	"""Mapping to the plain index with --sjdbGTFfile/--sjdbOverhang 99 inserts
	the annotated junctions on the fly and gives the same records and
	SJ.out.tab as the index built with the GTF, including reads with a
	4-base overhang (below alignSJoverhangMin 5, so only alignable through the
	annotation), which get their true spliced alignment."""
	short = [sim.transcript_read("T01a", 300 - 96), sim.transcript_read("T02", 500 - 4),
		sim.transcript_read("T05", 600 - 96, sense=False)]
	reads = sim.transcript_reads(150) + ReadSet(short)
	gtf = ["--sjdbGTFfile", str(sim.g.gtf_path), "--sjdbOverhang", "99"]
	fly = run_star(genome_plain, reads, gtf)
	pre = run_star(genome_sjdb, reads)
	assert record_multiset(fly.records()) == record_multiset(pre.records())
	assert fly.sj() == pre.sj()
	recs = pre.alignments().by_name()
	for r in short:
		assert [record_truth(x) for x in recs[r.name]] == [r.truth()]


def test_twopass_basic_aligns_short_overhangs_to_novel_junction(genome_plain, sim):
	"""The T10 exon-skipping junction chr2:30401-31800 is not annotated.
	Eight reads cross it with long overhangs; three cross it with 4 bases on
	one side, below alignSJoverhangMin 5, so 1-pass mapping cannot splice
	them. --twopassMode Basic inserts the 1st-pass junction, after which the
	short reads align spliced (alignSJDBoverhangMin 3) with their true CIGAR,
	and SJ.out.tab reports the junction as annotated (manual 5.5) with all
	11 unique reads."""
	long_reads = [sim.transcript_read("T10skip", o) for o in range(325, 365, 5)]
	short_reads = [sim.transcript_read("T10skip", 400 - 96), sim.transcript_read("T10skip", 400 - 4),
		sim.transcript_read("T10skip", 400 - 96, sense=False)]
	reads = ReadSet(long_reads + short_reads)
	one = run_star(genome_plain, reads).alignments().by_name()
	for r in short_reads:
		assert all("N" not in rec.cigar for rec in one.get(r.name, []))
	two = run_star(genome_plain, reads, "--twopassMode Basic")
	recs = two.alignments().by_name()
	for r in long_reads + short_reads:
		assert [record_truth(x) for x in recs[r.name]] == [r.truth()]
	sj = {(x.contig, x.start, x.end): x for x in two.sj()}
	assert sj[("chr2", 30401, 31800)][3:8] == (1, 1, 1, 11, 0)


def test_gene_counts_hand_worked(genome_sjdb, sim):
	"""ReadsPerGene.out.tab for a hand-built read set. Columns: unstranded,
	read strand = gene strand, read strand = opposite of gene.
	5 sense reads in G01 (+): G01 (5, 5, 0); col 3 sees them as no feature.
	4 antisense reads in G02 (-): G02 (4, 0, 4); col 2 no feature.
	3 + strand reads inside the G05 (+) / G06 (-) exon overlap 21201-21400:
	ambiguous unstranded, G05 in col 2, G06 in col 3.
	2 intergenic reads: no feature in all columns. 3 dup3 reads: multimapping.
	2 random reads: unmapped."""
	g01 = [sim.transcript_read("T01a", o) for o in (5, 40, 80, 120, 160)]
	g02 = [sim.transcript_read("T02", o, sense=False) for o in (10, 50, 90, 130)]
	ovl = [sim.genomic_read("chr1", p, "+") for p in (21205, 21250, 21290)]
	inter = [sim.genomic_read("chr1", p, "+") for p in (44500, 46500)]
	multi = [sim.segment_read("dup3", o) for o in (10, 100, 200)]
	rnd = [sim.random_read() for _ in range(2)]
	reads = ReadSet(g01 + g02 + ovl + inter + multi + rnd)
	counts = run_star(genome_sjdb, reads, "--quantMode GeneCounts").gene_counts()
	expected = {"N_unmapped": (2, 2, 2), "N_multimapping": (3, 3, 3), "N_noFeature": (2, 6, 7),
		"N_ambiguous": (3, 0, 0), "G01": (5, 5, 0), "G02": (4, 0, 4), "G05": (0, 3, 0), "G06": (0, 0, 3)}
	genes = [g.gene_id for g in sim.g.annotated_genes()]
	assert list(counts) == ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"] + genes
	for k, v in counts.items():
		assert v == expected.get(k, (0, 0, 0)), k


def test_transcriptome_bam_coordinates(genome_sjdb, sim):
	"""--quantMode TranscriptomeSAM writes Aligned.toTranscriptome.out.bam with
	alignments in transcript coordinates: RNAME the transcript, POS =
	transcript offset + 1 (5' to 3'), CIGAR 100M, FLAG 0 for sense and 16 for
	antisense reads, for + and - strand and spliced transcripts. A read in
	exon 1 of G01, shared by T01a and T01b, gets one record per transcript
	with NH 2."""
	reads = ReadSet([sim.transcript_read("T03", 350), sim.transcript_read("T03", 100, sense=False),
		sim.transcript_read("T02", 350), sim.transcript_read("T02", 20, sense=False),
		sim.transcript_read("T08", 280), sim.transcript_read("T01a", 30, name="shared")])
	run = run_star(genome_sjdb, reads, "--quantMode TranscriptomeSAM")
	tx = run.transcriptome().by_name()
	assert sorted(tx) == sorted(reads.by_name)
	for name, recs in tx.items():
		truth = reads.by_name[name]
		if name == "shared":
			assert sorted((r.rname, r.pos) for r in recs) == [("T01a", 31), ("T01b", 31)]
			assert all(r.tag("NH") == 2 and r.flag in (0, 256) for r in recs)
			continue
		assert len(recs) == 1
		rec = recs[0]
		assert (rec.rname, rec.pos, rec.cigar) == (truth.tx, truth.extra["tx_offset"] + 1, "100M")
		assert rec.flag == (0 if truth.extra["sense"] else 16)
		assert rec.seq == (sim.g.transcript_seq(truth.tx)[truth.extra["tx_offset"]:][:100])


def test_log_final_counts_hand_worked(genome_plain, sim):
	"""Log.final.out for 100 unique reads (20 spliced with one junction each),
	15 dup3 reads, 10 rep12 reads and 25 random reads: 150 input reads of
	length 100; 100 unique (66.67%); 15 multi (10.00%); 10 too many loci
	(6.67%); the 25 random reads are split between the three unmapped
	categories; 20 splices, all GT/AG-type here; mismatch rate 0."""
	uniq = sim.transcript_reads(80, tx_ids=["T07"]) + ReadSet(
		[r for r in sim.transcript_reads(60, tx_ids=["T01a", "T02", "T05", "T06"], spliced=True)
			if len(r.introns) == 1][:20])
	assert len(uniq) == 100
	other = ReadSet([sim.segment_read("dup3", sim.rng.randint(0, 500)) for _ in range(15)]
		+ [sim.segment_read("rep12", sim.rng.randint(0, 200)) for _ in range(10)]
		+ [sim.random_read() for _ in range(25)])
	log = run_star(genome_plain, uniq + other).log_final()
	assert log["Number of input reads"] == 150
	assert log["Average input read length"] == 100
	assert log["Uniquely mapped reads number"] == 100
	assert log["Uniquely mapped reads %"] == 66.67
	assert log["Number of reads mapped to multiple loci"] == 15
	assert log["% of reads mapped to multiple loci"] == 10.0
	assert log["Number of reads mapped to too many loci"] == 10
	assert log["% of reads mapped to too many loci"] == 6.67
	unmapped = sum(log["Number of reads unmapped: %s" % k] for k in ("too many mismatches", "too short", "other"))
	assert unmapped == 25
	assert log["Number of splices: Total"] == 20
	assert log["Number of splices: GT/AG"] == 20
	assert log["Number of splices: Annotated (sjdb)"] == 0
	assert log["Mismatch rate per base, %"] == 0.0
	assert log["Average mapped length"] == 100.0
	assert log["Number of chimeric reads"] == 0
