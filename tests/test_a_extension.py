"""Area A, extension and stitching limits: --alignEndsType, soft-clipping at
reference ends, junction overhang minimums, --alignIntronMin/Max,
--alignMatesGapMax and --alignSplicedMateMapLmin(OverLmate). Expected
alignments come from the simulation truth or from the local-extension rule
in helpers_a (extendAlign.cpp: +1 per match, -1 per mismatch, keep the best
strictly positive prefix)."""

import pytest

import helpers_a as ha
import star_helpers as sh
from star_helpers import ReadSet, SimPair, record_truth, run_star

ATTRS = ["--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD"]

# mismatch offsets (genome orientation) near the read ends
END_PATTERNS = [(99,), (98,), (97,), (0,), (1,), (2,), (0, 99), (1, 98), (0, 2), (97, 99), (0, 1, 2), (98, 99)]

# (extend left end of a '+' read, extend right end of a '+' read) per mode
# for a single-end read (mate 1); parametersDefault alignEndsType.
SE_FORCED = {
	"Local": (False, False),
	"EndToEnd": (True, True),
	"Extend5pOfRead1": (True, False),
	"Extend5pOfReads12": (True, False),
}


def _expected_end_alignment(synth, read, force_left, force_right):
	"""(pos, cigar, AS, nM) of an unspliced read with mismatches only near
	its ends, from helpers_a.end_clips."""
	ref = synth.seq(read.contig, read.pos, read.pos + 99)
	fwd = read.sam_seq()
	off = [i for i in range(100) if fwd[i] != ref[i]]
	left, right = ha.end_clips(fwd, ref, [i for i in off if i < 50], [i for i in off if i >= 50], force_left, force_right)
	nmm = sum(1 for i in off if left <= i < 100 - right)
	m = 100 - left - right
	cig = sh.ops_to_cigar([(left, "S"), (m, "M"), (right, "S")])
	return read.pos + left, cig, m - 2 * nmm + ha.length_score(m), nmm


@pytest.mark.parametrize("mode", list(SE_FORCED))
def test_end_mismatches_by_align_ends_type(genome_plain, sim, mode):
	"""Single-end reads with mismatches in their first or last three bases.
	Local soft-clips an end unless extending over it gains score; EndToEnd
	extends both ends; Extend5pOfRead1 and Extend5pOfReads12 extend the 5'
	end of read 1 (the left end of a '+' read, the right end of a '-' read)
	and treat the 3' end locally. Example (Local): a mismatch at the last base
	gives 99M1S, AS 99 - 2 = 97; at the second-last base 98M2S (extending
	scores -1 + 1 = 0, not > 0); at the third-last base 100M with nM 1, AS 98
	- 2 = 96. NM/MD follow the CIGAR."""
	reads = []
	for strand in "+-":
		for mm in END_PATTERNS:
			reads.append(sim.genomic_read(*sim.unique_locus(120), strand=strand, mismatches=mm))
	reads = ReadSet(reads)
	aln = run_star(genome_plain, reads, ATTRS + ["--alignEndsType", mode]).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	clipped = 0
	for rec in aln:
		r = reads.by_name[rec.qname]
		fl, fr = SE_FORCED[mode]
		if r.strand == "-":
			fl, fr = fr, fl
		pos, cig, score, nmm = _expected_end_alignment(sim.g, r, fl, fr)
		assert (rec.rname, rec.pos, rec.cigar, rec.strand) == (r.contig, pos, cig, r.strand), r.mismatches
		assert (rec.tag("AS"), rec.tag("nM")) == (score, nmm)
		assert (rec.tag("NM"), rec.tag("MD")) == sh.expected_nm_md(sim.g, r.contig, pos, cig, rec.seq)
		clipped += "S" in cig
	assert clipped == {"Local": 20, "EndToEnd": 0, "Extend5pOfRead1": 12, "Extend5pOfReads12": 12}[mode]


@pytest.mark.parametrize("mode,ext5", [
	("Local", (False, False)),
	("EndToEnd", (True, True)),
	("Extend5pOfRead1", (True, False)),
	("Extend5pOfReads12", (True, True)),
])
def test_mate_5p_extension_by_align_ends_type(genome_plain, sim, mode, ext5):
	"""Pairs (insert 300) with a mismatch at both end bases of each mate:
	Local clips them (1S98M1S); EndToEnd keeps every base (100M);
	Extend5pOfRead1 extends only mate 1's 5' end; Extend5pOfReads12 extends
	the 5' ends of both mates. ext5 = (mate 1, mate 2) 5' end extended. The
	3' ends stay local except under EndToEnd."""
	pairs = ReadSet([sim.genomic_pair(*sim.unique_locus(400), insert=300, strand="+-"[i % 2],
		mismatches1=(0, 99), mismatches2=(0, 99)) for i in range(8)])
	aln = run_star(genome_plain, pairs, ATTRS + ["--alignEndsType", mode]).alignments()
	by_name = aln.by_name()
	assert sorted(by_name) == sorted(pairs.by_name)
	for name, recs in by_name.items():
		p = pairs.by_name[name]
		kept = 0
		for rec in recs:
			truth = p.r1 if rec.mate == 1 else p.r2
			five = ext5[rec.mate - 1]
			three = mode == "EndToEnd"
			left, right = (not five, not three) if truth.strand == "+" else (not three, not five)
			cig = sh.ops_to_cigar([(int(left), "S"), (100 - left - right, "M"), (int(right), "S")])
			assert record_truth(rec) == (truth.contig, truth.pos + left, cig, truth.strand)
			kept += 2 - left - right
		assert all(r.tag("nM") == kept for r in recs)


def test_soft_clip_at_reference_ends(genome_plain, sim):
	"""Reads with 10 bases past the start of chr2 or the end of chr1/chr3
	align with those bases soft-clipped under the default
	--alignSoftClipAtReferenceEnds Yes (10S90M / 90M10S at the reference
	end). With No they are not output and are counted as unmapped
	(parametersDefault: "prohibit"). Reads that end exactly on the last base
	or start on the first base need no clipping and align in both modes."""
	g = sim.g
	reads = []
	extra = "ACGTTGCAAC"
	for i, (c, side) in enumerate([("chr2", "start"), ("chr1", "end"), ("chr3", "end"), ("chr3", "start")]):
		n = len(g.contigs[c])
		for strand in "+-":
			name = "%s_%s_%s" % (c, side, "p" if strand == "+" else "m")
			if side == "start":
				reads.append(ha.oriented(name + "_over", extra + g.seq(c, 1, 90), c, 1, "10S90M", strand))
				reads.append(ha.oriented(name + "_flush", g.seq(c, 1, 100), c, 1, "100M", strand))
			else:
				reads.append(ha.oriented(name + "_over", g.seq(c, n - 89, n) + extra, c, n - 89, "90M10S", strand))
				reads.append(ha.oriented(name + "_flush", g.seq(c, n - 99, n), c, n - 99, "100M", strand))
	reads = ReadSet(reads)
	aln = run_star(genome_plain, reads).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		assert record_truth(rec) == reads.by_name[rec.qname].truth()
	run = run_star(genome_plain, reads, "--alignSoftClipAtReferenceEnds No")
	aln = run.alignments()
	flush = sorted(n for n in reads.by_name if n.endswith("_flush"))
	assert sorted(aln.by_name()) == flush
	for rec in aln:
		assert record_truth(rec) == reads.by_name[rec.qname].truth()
	assert sum(ha.unmapped_counts(run.log_final())) == len(reads) - len(flush)


def _overhang_reads(sim, txs, overhangs):
	"""Reads over the first intron of each transcript with `k` bases on one
	side: '<tx>_R<k>' has k bases right of the intron, '<tx>_L<k>' k bases
	left of it; senses alternate."""
	reads = []
	for tx in txs:
		for k in overhangs:
			reads.append(ha.tx_junction_read(sim, tx, 0, 100 - k, sense=True, name="%s_R%d" % (tx, k)))
			reads.append(ha.tx_junction_read(sim, tx, 0, k, sense=False, name="%s_L%d" % (tx, k)))
	return ReadSet(reads)


def _short_side_clipped(synth, read):
	"""(pos, cigar) of a junction read whose short side is not spliced: the
	long block plus the local extension of the short side into the intron."""
	(s1, e1), (s2, e2) = sh.cigar_blocks(read.pos, read.cigar)
	score, pos, cig = ha.unspliced_fallback(synth.contigs[read.contig], read.sam_seq(), s1, e1 - s1 + 1, s2)[0]
	return pos, cig


def test_sj_overhang_min_unannotated(genome_plain, sim):
	"""Without annotation a junction needs blocks of at least
	--alignSJoverhangMin bases on both sides (stitchWindowAligns.cpp; the
	junctions here cannot slide, so the repeat allowance is 0). 6- and 7-base
	overhangs are spliced with the default 5; with --alignSJoverhangMin 7 the
	6-base side is no longer spliced: the read aligns with its long block,
	extended into the intron by the local rule, the rest soft-clipped."""
	reads = _overhang_reads(sim, ["T09", "T10"], (6, 7))
	for args in ([], ["--alignSJoverhangMin", "7"]):
		aln = run_star(genome_plain, reads, args).alignments()
		assert sorted(aln.by_name()) == sorted(reads.by_name)
		for rec in aln:
			r = reads.by_name[rec.qname]
			if args and rec.qname.endswith("6"):
				pos, cig = _short_side_clipped(sim.g, r)
				assert (rec.rname, rec.pos, rec.cigar, rec.strand) == (r.contig, pos, cig, r.strand)
			else:
				assert record_truth(rec) == r.truth()


def test_sjdb_overhang_min(genome_sjdb, sim):
	"""Annotated junctions are spliced with overhangs down to
	--alignSJDBoverhangMin (default 3) through the inserted junction
	sequences; with --alignSJDBoverhangMin 5, 3- and 4-base overhangs are not
	spliced (long block plus local extension, rest soft-clipped) while 5-base
	overhangs still are."""
	reads = _overhang_reads(sim, ["T01a", "T05", "T09", "T10"], (3, 4, 5))
	for args in ([], ["--alignSJDBoverhangMin", "5"]):
		aln = run_star(genome_sjdb, reads, args).alignments()
		assert sorted(aln.by_name()) == sorted(reads.by_name)
		for rec in aln:
			r = reads.by_name[rec.qname]
			if args and int(rec.qname[-1]) < 5:
				pos, cig = _short_side_clipped(sim.g, r)
				assert (rec.rname, rec.pos, rec.cigar, rec.strand) == (r.contig, pos, cig, r.strand)
			else:
				assert record_truth(rec) == r.truth()


def _gap_reads(sim, n_each=4, gap=22):
	"""Reads with a `gap`-base genomic gap after 50 bases, half with GT...AG
	at the gap ends and half non-canonical; gap position cannot slide and
	the junction scan (helpers_a.junction_placement) keeps it in place."""
	g = sim.g
	reads = []
	kinds = {"gtag": 0, "nonc": 0}
	while min(kinds.values()) < n_each:
		c, p = sim.unique_locus(150)
		s = p + 50
		if not sim.deletion_is_unambiguous(c, s, gap):
			continue
		code = sh.motif_code(g.seq(c, s, s + 1), g.seq(c, s + gap - 2, s + gap - 1))
		kind = "gtag" if code == 1 else ("nonc" if code == 0 else None)
		if kind is None or kinds[kind] >= n_each:
			continue
		r = sim.genomic_read(c, p, "+-"[len(reads) % 2], deletion=(50, gap), name="%s%d" % (kind, kinds[kind]))
		if ha.junction_placement(g.contigs[c], r.sam_seq(), 50, s - 1, s + gap)[0] != 0:
			continue
		kinds[kind] += 1
		reads.append(r)
	return ReadSet(reads)


def test_intron_min_separates_deletions_from_introns(genome_plain, sim):
	"""A 22-base genomic gap is an intron when 22 >= --alignIntronMin
	(parametersDefault) and a deletion otherwise. --alignIntronMin 22:
	50M22N50M, AS = 100 - 2 (length term) + motif penalty (0 for GT/AG, -8
	non-canonical) = 98 or 90. --alignIntronMin 23 with --scoreDelBase 0:
	50M22D50M, AS = 100 - 2 (scoreDelOpen) - 2 = 96. --alignIntronMin 23 with
	the default deletion penalties: 100 - 2 - 44 - 2 = 52 for the deletion
	and 49 for 50M50S, both below int(0.66 * 99): unmapped as too short."""
	reads = _gap_reads(sim)
	aln = run_star(genome_plain, reads, ATTRS + ["--alignIntronMin", "22"]).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		r = reads.by_name[rec.qname]
		assert record_truth(rec) == (r.contig, r.pos, "50M22N50M", r.strand)
		assert rec.tag("AS") == (98 if rec.qname.startswith("gtag") else 90)
	aln = run_star(genome_plain, reads, ATTRS + ["--alignIntronMin", "23", "--scoreDelBase", "0"]).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		r = reads.by_name[rec.qname]
		assert record_truth(rec) == r.truth()
		assert rec.tag("AS") == 96 and rec.tag("NM") == 22
	run = run_star(genome_plain, reads, ATTRS + ["--alignIntronMin", "23"])
	assert len(run.alignments()) == 0
	assert ha.unmapped_counts(run.log_final()) == (0, len(reads), 0, 0)


def test_intron_max(genome_plain, sim):
	"""--alignIntronMax 1000: the 1000-base intron of T01a is allowed (a gap
	is rejected only when larger than the maximum), the 5000-base intron of
	T09 is not. T09 reads are then aligned by their longer block extended
	into the intron (helpers_a.unspliced_fallback); a 50/50 read keeps 49-ish
	points, below int(0.66 * 99), and is unmapped as too short."""
	reads = []
	for tx in ("T01a", "T09"):
		for left in (70, 50, 30):
			for sense in (True, False):
				reads.append(ha.tx_junction_read(sim, tx, 0, left, sense=sense, name="%s_%d_%d" % (tx, left, sense)))
	reads = ReadSet(reads)
	run = run_star(genome_plain, reads, ATTRS + ["--alignIntronMax", "1000"])
	got = {r.qname: r for r in run.alignments()}
	short = 0
	for r in reads:
		if r.tx == "T01a":
			assert record_truth(got[r.name]) == r.truth()
			continue
		(s1, e1), (s2, e2) = sh.cigar_blocks(r.pos, r.cigar)
		best = ha.unspliced_fallback(sim.g.contigs[r.contig], r.sam_seq(), s1, e1 - s1 + 1, s2)
		score, pos, cig = best[0]
		assert score != 65
		if score >= 66:
			assert best[1][0] < score
			assert (got[r.name].pos, got[r.name].cigar, got[r.name].tag("AS")) == (pos, cig, score)
		else:
			assert r.name not in got
			short += 1
	assert short == 2
	assert ha.unmapped_counts(run.log_final()) == (0, 2, 0, 0)


def test_mates_gap_max(genome_plain, sim):
	"""--alignMatesGapMax 100: pairs with 50 or 100 bases between the mates
	map as proper pairs at their true positions; with 101 or 150 bases the
	mates cannot be stitched (stitchAlignToTranscript.cpp rejects a gap
	larger than the maximum) and one mate alone (score 98) is below
	int(0.66 * 200): unmapped as too short. Without the option all pairs map."""
	pairs = []
	for gap in (50, 100, 101, 150):
		for strand in "+-":
			pairs.append(sim.genomic_pair(*sim.unique_locus(400), insert=200 + gap, strand=strand,
				name="gap%d%s" % (gap, "p" if strand == "+" else "m")))
	pairs = ReadSet(pairs)
	for args, ok in (([], (50, 100, 101, 150)), (["--alignMatesGapMax", "100"], (50, 100))):
		run = run_star(genome_plain, pairs, args)
		by_name = run.alignments().by_name()
		want = sorted(n for n in pairs.by_name if int(n[3:-1]) in ok)
		assert sorted(by_name) == want
		for name in want:
			p = pairs.by_name[name]
			recs = {r.mate: r for r in by_name[name]}
			assert record_truth(recs[1]) == p.r1.truth() and record_truth(recs[2]) == p.r2.truth()
			assert recs[1].is_proper_pair
		assert ha.unmapped_counts(run.log_final()) == (0, len(pairs) - len(want), 0, 0)


def test_spliced_mate_map_min(genome_plain, sim):
	"""Reads of 25 Ns, then 45 + 30 bases across an intron: the spliced
	alignment maps 75 bases. stitchWindowAligns.cpp drops a spliced
	alignment whose mate maps fewer than --alignSplicedMateMapLmin or
	(uint)(--alignSplicedMateMapLminOverLmate * mate length) bases. Limits of
	75 (absolute) and 0.75 keep it (25S45M<N>30M); 76 and 0.76 drop it, and
	the unspliced remainder (45 bases) is too short, so the read is unmapped."""
	reads = []
	for i, tx in enumerate(["T05", "T09", "T10", "T02", "T03", "T08"]):
		r = ha.tx_junction_read(sim, tx, 0, 70, sense=True)
		fwd = "N" * 25 + r.sam_seq()[25:]
		(s1, e1), (s2, e2) = sh.cigar_blocks(r.pos, r.cigar)
		cig = "25S45M%dN30M" % (s2 - e1 - 1)
		reads.append(ha.oriented("sm%d" % i, fwd, r.contig, r.pos + 25, cig, "+-"[i % 2]))
	reads = ReadSet(reads)
	for args, mapped in (("--alignSplicedMateMapLmin 75", True), ("--alignSplicedMateMapLmin 76", False),
			("--alignSplicedMateMapLminOverLmate 0.75", True), ("--alignSplicedMateMapLminOverLmate 0.76", False)):
		run = run_star(genome_plain, reads, args)
		aln = run.alignments()
		if mapped:
			assert sorted(aln.by_name()) == sorted(reads.by_name)
			for rec in aln:
				assert record_truth(rec) == reads.by_name[rec.qname].truth()
		else:
			assert len(aln) == 0
			assert ha.unmapped_counts(run.log_final()) == (0, len(reads), 0, 0)


@pytest.mark.parametrize("dist,spliced_5kb", [(9, False), (30, True)])
def test_window_bins_bound_intron_length(genome_plain, sim, dist, spliced_5kb):
	"""With --alignIntronMax 0 the longest intron is about
	2^winBinNbits * winAnchorDistNbins (parametersDefault): with
	--winBinNbits 8, 9 bins (2304 bases) still join the 400- and 1000-base
	introns of T03 and T01a but not the 5000-base T09 intron, whose reads
	then align with their longer block extended into the intron
	(helpers_a.unspliced_fallback); 30 bins (7680 bases) join all three."""
	reads = ReadSet([ha.tx_junction_read(sim, tx, 0, left, sense=bool(i), name="%s_%d_%d" % (tx, left, i))
		for tx in ("T01a", "T03", "T09") for left in (70, 30) for i in range(2)])
	got = {r.qname: r for r in run_star(genome_plain, reads, ATTRS + ["--winBinNbits", "8", "--winAnchorDistNbins",
		str(dist)]).alignments()}
	assert sorted(got) == sorted(reads.by_name)
	for r in reads:
		rec = got[r.name]
		if r.tx != "T09" or spliced_5kb:
			assert record_truth(rec) == r.truth()
		else:
			(s1, e1), (s2, e2) = sh.cigar_blocks(r.pos, r.cigar)
			score, pos, cig = ha.unspliced_fallback(sim.g.contigs[r.contig], r.sam_seq(), s1, e1 - s1 + 1, s2)[0]
			assert (rec.pos, rec.cigar, rec.tag("AS")) == (pos, cig, score)


def _occurs(synth, s):
	"""True if s occurs on either strand of the genome."""
	rc = sh.revcomp(s)
	return any(s in seq or rc in seq for seq in synth.contigs.values())


@pytest.mark.parametrize("seed_map_min", [5, 6])
def test_seed_map_min(genome_plain, sim, seed_map_min):
	"""--seedMapMin is the "min length of seeds to be mapped"
	(parametersDefault): a k-base overhang past a junction gets its own seed
	when k > seedMapMin and none when k < seedMapMin. The search from the
	other read end does not find it either: the k overhang bases plus one
	more occur elsewhere in the genome (checked), so the longest match from
	that end is not at the true locus. --alignSJoverhangMin 3 lets every
	overhang here be spliced once seeded. Overhangs of seedMapMin - 1 bases
	are therefore not spliced (long block plus local extension, rest
	soft-clipped) and overhangs of seedMapMin + 1 bases are. k = seedMapMin
	is left out: STAR does not search a piece of exactly seedMapMin bases,
	which the documented minimum would allow."""
	reads = _overhang_reads(sim, ["T09", "T10", "T05"], (seed_map_min - 1, seed_map_min + 1))
	args = ["--alignSJoverhangMin", "3"] + ([] if seed_map_min == 5 else ["--seedMapMin", str(seed_map_min)])
	aln = run_star(genome_plain, reads, args).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		r = reads.by_name[rec.qname]
		k = int(rec.qname.split("_")[1][1:])
		if k > seed_map_min:
			assert record_truth(rec) == r.truth()
			continue
		fwd = r.sam_seq()
		assert _occurs(sim.g, fwd[-(k + 1):] if "_R" in rec.qname else fwd[:k + 1])
		pos, cig = _short_side_clipped(sim.g, r)
		assert (rec.rname, rec.pos, rec.cigar, rec.strand) == (r.contig, pos, cig, r.strand)
