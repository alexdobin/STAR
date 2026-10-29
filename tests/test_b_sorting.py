"""Area B: coordinate sorting (order, ties, unmapped reads, sorting bins,
threads and RAM), BAM compression, --outStd, SAM header options, and the
SAM, unsorted-BAM and sorted-BAM writers giving the same records."""

import gzip
import subprocess

import pytest

import helpers_b as hb
import star_helpers as sh
from star_helpers import ReadSet, record_multiset, record_truth, run_star

BOTH_BAMS = ["--outSAMtype", "BAM", "Unsorted", "SortedByCoordinate"]
SINGLE_MATE_OK = ["--outFilterScoreMinOverLread", "0.3", "--outFilterMatchNminOverLread", "0.3"]


def _sort_key(aln):
	index = {name: i for i, (name, _) in enumerate(aln.references)}
	return [(index[r.rname], r.pos) if not r.is_unmapped else (len(index), 0) for r in aln]


def _lines(aln):
	return [r.to_line() for r in aln]


def _tie_reads(sim):
	"""Reads whose order in the input is known: three reads at the same
	position (named so that name order differs from input order), reads one
	base either side, unique and 3-copy reads, and unmapped reads spread
	through the file."""
	c, p = sim.unique_locus(150)
	reads = [sim.genomic_read(c, p + 1, "+", name="right")]
	reads += [sim.genomic_read(c, p, s, name=n) for n, s in (("tieC", "+"), ("tieA", "-"), ("tieB", "+"))]
	reads.append(sim.genomic_read(c, p - 1, "-", name="left"))
	for i in range(30):
		reads.append(sim.genomic_read(*sim.unique_locus(), strand=sim.rng.choice("+-")))
		if i % 10 == 0:
			reads.append(sim.random_read(name="unm%d" % (9 - i // 10)))
	reads += [sim.segment_read("dup3", sim.rng.randint(0, 500)) for _ in range(5)]
	return ReadSet(reads)


def test_sorted_bam_order_ties_and_unmapped(genome_plain, sim):
	"""Aligned.sortedByCoord.out.bam holds the same records as
	Aligned.out.bam (--outSAMtype BAM Unsorted SortedByCoordinate), ordered by
	reference (the @SQ order) and POS, with unmapped records (--outSAMunmapped
	Within) after all mapped ones and the header line @HD VN:1.4
	SO:coordinate (SAM specification 1.3); the unsorted file has no SO.
	Records at the same position come in input order, whatever their names,
	and unmapped records in input order (pins the current output: ties are
	broken by read index, BAMbinSortByCoordinate.cpp and
	BAMbinSortUnmapped.cpp). Paired-end: the same checks, every mate pair
	keeping RNEXT/PNEXT/TLEN."""
	reads = _tie_reads(sim)
	run = run_star(genome_plain, reads, BOTH_BAMS + ["--outSAMunmapped", "Within"])
	uns, srt = run.alignments("Aligned.out.bam"), run.alignments("Aligned.sortedByCoord.out.bam")
	assert record_multiset(srt) == record_multiset(uns)
	assert len(srt) == len(reads) + 10
	keys = _sort_key(srt)
	assert keys == sorted(keys)
	assert srt.header[0] == "@HD\tVN:1.4\tSO:coordinate" and uns.header[0] == "@HD\tVN:1.4"
	input_index = {r.name: i for i, r in enumerate(reads)}
	for a, b, ka, kb in zip(srt.records, srt.records[1:], keys, keys[1:]):
		if ka == kb:
			assert input_index[a.qname] <= input_index[b.qname]
	names = [r.qname for r in srt]
	start = names.index("left")
	assert names[start:start + 5] == ["left", "tieC", "tieA", "tieB", "right"]
	assert [r.qname for r in srt if r.is_unmapped] == ["unm9", "unm8", "unm7"]

	pairs = sim.transcript_pairs(30) + ReadSet([sim.genomic_pair(*sim.unique_locus(150), 100) for _ in range(4)])
	run = run_star(genome_plain, pairs, BOTH_BAMS)
	uns, srt = run.alignments("Aligned.out.bam"), run.alignments("Aligned.sortedByCoord.out.bam")
	assert record_multiset(srt) == record_multiset(uns)
	keys = _sort_key(srt)
	assert keys == sorted(keys)
	for rec in srt:
		truth = getattr(pairs.by_name[rec.qname], "r%d" % rec.mate)
		assert record_truth(rec) == truth.truth()


def _many_reads(sim):
	"""About 600 records: transcript and unique genomic reads, 3-copy reads,
	random reads, and 30 pairs of reads at the same position, one of each
	pair in the first half of the input and the other at the end."""
	tx = list(sim.transcript_reads(400))
	ties = []
	for i in range(30):
		c, p = sim.unique_locus()
		ties.append((sim.genomic_read(c, p, "+", name="tieA%02d" % i), sim.genomic_read(c, p, "-", name="tieB%02d" % i)))
	return ReadSet(tx[:200] + [a for a, _ in ties] + tx[200:]
		+ [sim.genomic_read(*sim.unique_locus(), strand=sim.rng.choice("+-")) for _ in range(30)]
		+ [sim.segment_read("dup3", sim.rng.randint(0, 500), strand=sim.rng.choice("+-")) for _ in range(20)]
		+ [sim.random_read() for _ in range(20)] + [b for _, b in ties])


def test_sorted_bam_invariant_to_sorting_options(genome_sjdb, sim):
	"""The decoded sorted BAM (header without the command-line lines, and
	every record in order) is the same with 1 and 4 mapping threads (also
	with --limitIObufferSize 110000 50000000, which splits the input into
	chunks of about 35 reads mapped by different threads, so reads at the
	same position come from different threads), --outBAMsortingThreadN 1 and
	4, --outBAMsortingBinsN 3 and 100,
	--limitBAMsortRAM one byte above the 'Max memory needed for sorting' of
	the reference run (bins then wait for each other), and a 50 kB output
	buffer (--limitIObufferSize 30000000 50000), which makes STAR fix the bin
	boundaries after the first ~200 records and spill bins to disk many
	times. None of these options is documented to change the output."""
	reads = _many_reads(sim)
	base = ["--outSAMtype", "BAM", "SortedByCoordinate", "--outSAMunmapped", "Within"]
	ref = run_star(genome_sjdb, reads, base, threads=1)
	want = (hb.stable_header(ref.alignments().header), _lines(ref.alignments()))
	assert len(want[1]) > 550
	keys = _sort_key(ref.alignments())
	assert keys == sorted(keys)
	variants = [
		([], 4),
		(["--limitIObufferSize", "110000", "50000000"], 4),
		(["--outBAMsortingThreadN", "1"], 4),
		(["--outBAMsortingThreadN", "4"], 1),
		(["--outBAMsortingBinsN", "3"], None),
		(["--outBAMsortingBinsN", "100"], None),
		(["--limitBAMsortRAM", str(hb.max_sort_memory(ref) + 1)], None),
		(["--limitIObufferSize", "30000000", "50000", "--outBAMsortingBinsN", "50"], None),
	]
	for extra, threads in variants:
		run = run_star(genome_sjdb, reads, base + extra, threads=threads)
		aln = run.alignments()
		assert (hb.stable_header(aln.header), _lines(aln)) == want, extra
		if "50000" in extra:
			n_binned = int(run.log_out().split("BAM sorting: ")[1].split()[0])
			assert 0 < n_binned < len(want[1]) - 40


def test_bam_compression_levels(genome_sjdb, sim):
	"""--outBAMcompression 0, 1 (default) and 9 (parametersDefault: 0 = no
	compression, 10 = maximum) give the same decoded records in both BAM
	files. Level 0 stores the data uncompressed, so the file is at least as
	large as its decompressed content; level 9 gives a smaller file than
	level 1, and level 1 than level 0."""
	reads = _many_reads(sim)
	sizes, decoded = {}, {}
	for level in ("0", "1", "9"):
		run = run_star(genome_sjdb, reads, BOTH_BAMS + ["--outBAMcompression", level])
		for name in ("Aligned.out.bam", "Aligned.sortedByCoord.out.bam"):
			raw = run.path(name).read_bytes()
			sizes[level, name] = len(raw)
			if level == "0":
				assert len(raw) >= len(gzip.decompress(raw))
			aln = sh.read_bam(raw)
			lines = _lines(aln) if "sorted" in name else record_multiset(aln)
			decoded.setdefault(name, []).append((hb.stable_header(aln.header), lines))
	for name, versions in decoded.items():
		assert versions[0] == versions[1] == versions[2]
		assert sizes["9", name] < sizes["1", name] < sizes["0", name]


def test_out_std(genome_sjdb, sim):
	"""--outStd SAM / BAM_Unsorted / BAM_SortedByCoordinate (parametersDefault)
	send the alignments to stdout instead of Aligned.out.sam /
	Aligned.out.bam / Aligned.sortedByCoord.out.bam: decoded, stdout equals
	the file a normal run writes (header without the command-line lines; the
	same records, in the same order for the sorted BAM). The output file is
	not created and the log lines go to Log.std.out."""
	reads = sim.transcript_reads(60) + ReadSet([sim.random_read() for _ in range(5)])
	common = ["--outSAMunmapped", "Within"]
	sam = run_star(genome_sjdb, reads, common).alignments()
	bams = run_star(genome_sjdb, reads, common + BOTH_BAMS)
	cases = [
		("SAM", [], sam, "Aligned.out.sam"),
		("BAM_Unsorted", ["--outSAMtype", "BAM", "Unsorted"], bams.alignments("Aligned.out.bam"), "Aligned.out.bam"),
		("BAM_SortedByCoordinate", ["--outSAMtype", "BAM", "SortedByCoordinate"],
			bams.alignments("Aligned.sortedByCoord.out.bam"), "Aligned.sortedByCoord.out.bam"),
	]
	for mode, fmt, ref, fname in cases:
		run = run_star(genome_sjdb, reads, common + fmt + ["--outStd", mode])
		got = sh.parse_sam(run.stdout) if mode == "SAM" else sh.read_bam(run.stdout_bytes)
		assert hb.stable_header(got.header) == hb.stable_header(ref.header)
		if mode == "BAM_SortedByCoordinate":
			assert _lines(got) == _lines(ref)
		else:
			assert record_multiset(got) == record_multiset(ref)
		assert len(got) == len(reads)
		assert not run.exists(fname)
		assert "finished successfully" in run.path("Log.std.out").read_text()


@pytest.mark.parametrize("fmt", list(sh.OUTPUT_FORMATS))
def test_header_lines(genome_sjdb, sim, tmp_path, fmt):
	"""Header without options: @HD VN:1.4 (plus SO:coordinate in the sorted
	BAM), one @SQ per contig with its length, @PG ID:STAR PN:STAR VN:<STAR
	--version> CL:<command line>, @CO user command line. --outSAMheaderHD
	replaces the @HD line, --outSAMheaderPG adds an @PG line before STAR's,
	--outSAMheaderCommentFile adds its non-blank lines (parametersDefault),
	and --outSAMattrRGline adds @RG. @HD is first (SAM specification 1.3);
	the order of the other lines pins the current output."""
	star_version = subprocess.run([str(sh.find_star()), "--version"], stdout=subprocess.PIPE).stdout.decode().strip()
	so = ["SO:coordinate"] if fmt == "sorted" else []
	sq = ["@SQ\tSN:%s\tLN:%d" % (n, len(s)) for n, s in sim.g.contigs.items()]
	reads = sim.transcript_reads(10)

	def check_tail(header, n_fixed):
		pg, co = header[n_fixed], header[-1]
		assert pg.startswith("@PG\tID:STAR\tPN:STAR\tVN:%s\tCL:" % star_version)
		assert co.startswith("@CO\tuser command line: %s " % sh.find_star())
		assert "--outFileNamePrefix" in pg and "--outFileNamePrefix" in co

	h = run_star(genome_sjdb, reads, sh.OUTPUT_FORMATS[fmt]).alignments().header
	assert h[:4] == ["\t".join(["@HD", "VN:1.4"] + so)] + sq and len(h) == 6
	check_tail(h, 4)

	comments = tmp_path / "comments.txt"
	comments.write_text("@CO\tfirst comment\n\n   \n@CO\tsecond comment\n")
	args = ["--outSAMheaderHD", "@HD", "VN:1.6", "--outSAMheaderPG", "@PG", "ID:other", "PN:prep", "VN:0.1",
		"--outSAMheaderCommentFile", str(comments), "--outSAMattrRGline", "ID:g1", "SM:s1"]
	h = run_star(genome_sjdb, reads, args + sh.OUTPUT_FORMATS[fmt]).alignments().header
	assert h[:5] == ["\t".join(["@HD", "VN:1.6"] + so)] + sq + ["@PG\tID:other\tPN:prep\tVN:0.1"]
	assert h[6:-1] == ["@CO\tfirst comment", "@CO\tsecond comment", "@RG\tID:g1\tSM:s1"] and len(h) == 10
	check_tail(h, 5)
	assert "--outSAMheaderCommentFile" in h[5]


def test_tmp_dir_and_keep(genome_sjdb, sim, tmp_path):
	"""--outTmpDir puts STAR's temporary files (here the BAM sorting bins) in
	the given directory and --outTmpKeep None (default) removes it at the
	end; --outTmpKeep All keeps it, and without --outTmpDir the temporary
	directory is <prefix>_STARtmp (parametersDefault). The sorted BAM is the
	same in all three runs."""
	reads = sim.transcript_reads(50)
	base = ["--outSAMtype", "BAM", "SortedByCoordinate"]
	ref = run_star(genome_sjdb, reads, base)
	assert not ref.exists("_STARtmp")
	tmp = tmp_path / "startmp"
	gone = run_star(genome_sjdb, reads, base + ["--outTmpDir", str(tmp)])
	assert not tmp.exists()
	kept = run_star(genome_sjdb, reads, base + ["--outTmpKeep", "All"])
	assert (kept.path("_STARtmp") / "BAMsort").is_dir()
	tmp2 = tmp_path / "startmp2"
	kept2 = run_star(genome_sjdb, reads, base + ["--outTmpDir", str(tmp2), "--outTmpKeep", "All"])
	assert (tmp2 / "BAMsort").is_dir() and not kept2.exists("_STARtmp")
	want = _lines(ref.alignments())
	for run in (gone, kept, kept2):
		assert _lines(run.alignments()) == want


def _writer_reads(sim):
	se = _many_reads(sim)
	extra = []
	for i in range(10):
		c, p = sim.unique_locus(120)
		extra.append(sim.genomic_read(c, p, sim.rng.choice("+-"), mismatches=(20, 61)))
		c, p = sim.unique_locus(120)
		extra.append(sim.genomic_read(c, p, sim.rng.choice("+-"), insertion=(40, 2)))
		c, p = sim.unique_locus(120)
		if sim.deletion_is_unambiguous(c, p + 45, 2):
			extra.append(sim.genomic_read(c, p, sim.rng.choice("+-"), deletion=(45, 2)))
	se = se + ReadSet(extra)
	pe = sim.transcript_pairs(80) + ReadSet(
		[hb.half_pair("half%d" % i, sim.genomic_read(*sim.unique_locus(), strand=sim.rng.choice("+-"),
			name="half%d" % i), 1 + i % 2) for i in range(6)]
		+ [hb.half_pair("halfm%d" % i, sim.segment_read("dup3", 50 * i, name="halfm%d" % i), 1 + i % 2)
			for i in range(4)]
		+ [hb.unmappable_pair("unm%d" % i) for i in range(3)])
	return se, pe


def test_sam_and_bam_writers_agree(genome_sjdb, sim):
	"""SAM text (ReadAlign_outputTranscriptSAM.cpp) and BAM records
	(ReadAlign_alignBAM.cpp) are written by separate code; for the same
	options they hold the same multiset of records, unsorted and sorted:
	single-end reads (spliced, substitutions, indels, 3 loci, unmapped) with
	NH HI AS nM NM MD jM jI XS, and pairs (proper, one mate unmapped with 1
	or 3 loci, both unmapped) with NH HI AS nM NM MD jM jI MC, both with
	--outSAMunmapped Within. In SAM and unsorted BAM all records of a read
	are adjacent and the two mates of a pair alignment follow each other
	(manual 5.3). Each BAM record's bin field is reg2bin of its aligned
	interval (4680 = reg2bin(-1, 0) for unmapped records) and each BAM file
	ends with the BGZF end-of-file block (SAM specification 4.1.2, 4.2.1)."""
	se, pe = _writer_reads(sim)
	cases = [
		(se, ["--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD", "jM", "jI", "XS"]),
		(pe, SINGLE_MATE_OK + ["--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD", "jM", "jI", "MC"]),
	]
	for reads, args in cases:
		args = args + ["--outSAMunmapped", "Within"]
		sam = run_star(genome_sjdb, reads, args).alignments()
		bams = run_star(genome_sjdb, reads, args + BOTH_BAMS)
		uns = bams.alignments("Aligned.out.bam")
		want = record_multiset(sam)
		assert len(want) > len(reads)
		assert record_multiset(uns) == want
		srt = bams.alignments("Aligned.sortedByCoord.out.bam")
		assert record_multiset(srt) == want
		assert sum(r.is_unmapped for r in sam) >= 5
		for name in ("Aligned.out.bam", "Aligned.sortedByCoord.out.bam"):
			raw = bams.path(name).read_bytes()
			assert raw.endswith(hb.BGZF_EOF)
			for qname, pos, bin_, ref_len in hb.bam_bins(raw):
				assert bin_ == (hb.reg2bin(pos, pos + ref_len) if pos >= 0 else 4680), qname
		for aln in (sam, uns):
			assert hb.contiguous_by_name(aln.records)
			if reads.paired:
				proper = [r for r in aln if r.flag & 0x2]
				assert len(proper) == 160
				for x, y in zip(proper[::2], proper[1::2]):
					assert x.qname == y.qname and {x.mate, y.mate} == {1, 2}
					assert (x.pnext, y.pnext) == (y.pos, x.pos)
					assert aln.records.index(y) == aln.records.index(x) + 1


@pytest.mark.parametrize("fmt", ["bam", "sorted"])
def test_bam_seq_odd_length_and_lower_case(genome_plain, sim, fmt):
	"""BAM packs SEQ two bases per byte, the first base in the high nibble,
	so an odd-length read ends in the high nibble of a last, half-used byte,
	and the 16 codes "=ACMGRSVTWYHKDBN" have no lower case (SAM
	specification 4.2). Reads of 101 and 75 bases (odd), reads written in
	lower case or mixed case, and unmapped reads of 51 Ns (--outSAMunmapped
	Within): each decoded SEQ is the input sequence upper-cased
	(reverse-complemented for '-' alignments), QUAL is the input's, and the
	mapped reads are at their simulated alignment."""
	reads = [sim.genomic_read(*sim.unique_locus(110), strand="+-"[i % 2], length=101) for i in range(6)]
	for i in range(6):
		r = sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2], length=75)
		r.seq = r.seq.lower()
		reads.append(r)
	for i in range(4):
		r = sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2])
		r.seq = r.seq[:33].lower() + r.seq[33:66] + r.seq[66:].lower()
		reads.append(r)
	reads += [hb.n_read("nread%d" % i, length=51) for i in range(2)]
	reads = ReadSet(reads)
	recs = run_star(genome_plain, reads, ["--outSAMunmapped", "Within"] + sh.OUTPUT_FORMATS[fmt]).records()
	assert sorted(r.qname for r in recs) == sorted(reads.by_name)
	for rec in recs:
		truth = reads.by_name[rec.qname]
		if truth.kind == "unmappable":
			assert rec.is_unmapped and (rec.seq, rec.qual) == ("N" * 51, truth.qual)
			continue
		assert record_truth(rec) == truth.truth()
		assert rec.seq == truth.sam_seq().upper()
		assert rec.qual == truth.sam_qual()
