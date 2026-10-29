"""Area A, read input: FASTQ/FASTA parsing (CRLF line ends, multi-line
FASTA, lower-case and N bases), read names, --readMapNumber,
--outQSconversionAdd, --seedSplitMin, and the 5'/3' clipping options."""

import gzip

import pytest

import helpers_a as ha
import star_helpers as sh
from star_helpers import ReadSet, SimPair, record_truth, run_star

ATTRS = ["--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD"]


def _write(path, text, gz=False):
	data = text.encode()
	if gz:
		with gzip.open(str(path), "wb") as f:
			f.write(data)
	else:
		path.write_bytes(data)
	return str(path)


def _special_reads(sim):
	"""Unique reads plus three special ones: 'lower' (+ strand, written in
	lower case), 'twoN' (- strand, N at genome-orientation offsets 30 and
	70) and 'Nevery11' (+ strand, N at offsets 10, 21, ..., 98, so every
	N-free piece is 10 bases or shorter)."""
	reads = [sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2]) for i in range(8)]
	c, p = sim.unique_locus()
	low = sim.genomic_read(c, p, "+", name="lower")
	low.seq = low.seq.lower()
	reads.append(low)
	c, p = sim.unique_locus()
	reads.append(sim.genomic_read(c, p, "-", n_offsets=(30, 70), name="twoN"))
	c, p = sim.unique_locus()
	reads.append(sim.genomic_read(c, p, "+", n_offsets=tuple(range(10, 100, 11)), name="Nevery11"))
	return ReadSet(reads)


def _check_special(aln, reads, sim, fasta):
	"""Shared expectations of test_fastq_crlf_and_multiline_fasta. Nevery11
	has no N-free piece of seedSplitMin (12) bases, so it has no seeds and is
	unmapped. Ns score 0 (extendAlign.cpp skips them), so twoN scores 98 - 2 =
	96 with nM 0, while NM/MD count each N as a mismatch."""
	by_name = aln.by_name()
	assert sorted(by_name) == sorted(n for n in reads.by_name if n != "Nevery11")
	for name, recs in by_name.items():
		truth = reads.by_name[name]
		assert len(recs) == 1
		rec = recs[0]
		assert record_truth(rec) == truth.truth()
		assert rec.seq == truth.sam_seq()
		assert rec.qual == ("*" if fasta else truth.sam_qual())
		nm, md = sh.expected_nm_md(sim.g, truth.contig, truth.pos, truth.cigar, truth.sam_seq().upper())
		assert (rec.tag("NM"), rec.tag("MD"), rec.tag("nM")) == (nm, md, 0)
		assert rec.tag("AS") == (96 if name == "twoN" else 98)
	assert by_name["lower"][0].seq == reads.by_name["lower"].seq


def test_fastq_crlf_and_multiline_fasta(genome_plain, sim, tmp_path):
	"""FASTQ and FASTA written with CRLF line ends (CHANGES.md: control
	characters are removed from the ends of read lines, for files from
	Windows), the FASTA with the sequence wrapped at 60 bases (joined into one
	line by processChunks), give the true alignments. The FASTQ name lines
	carry a comment: without one, CRLF FASTQ fails, which is not tested
	here. Lower-case bases are read as A/C/G/T and written back as
	given; SEQ is the input sequence (reverse-complemented for '-')."""
	reads = _special_reads(sim)
	fq = "".join("@%s 1:N:0:1\r\n%s\r\n+\r\n%s\r\n" % (r.name, r.seq, r.qual) for r in reads)
	aln = run_star(genome_plain, [_write(tmp_path / "crlf.fq", fq)], ATTRS).alignments()
	_check_special(aln, reads, sim, fasta=False)
	fa = "".join(">%s\r\n%s\r\n%s\r\n" % (r.name, r.seq[:60], r.seq[60:]) for r in reads)
	aln = run_star(genome_plain, [_write(tmp_path / "wrapped.fa", fa)], ATTRS).alignments()
	_check_special(aln, reads, sim, fasta=True)


def test_seed_split_min(genome_plain, sim):
	"""qualitySplit cuts reads at N into pieces and only pieces of at least
	--seedSplitMin bases are searched for seeds. With an N every 11 bases the
	pieces are 10 bases long: --seedSplitMin 11 leaves no piece (unmapped,
	counted as 'other'), --seedSplitMin 10 maps every read at its origin with
	AS = 91 matches - 2 (length term) = 89, nM 0 and NM 9."""
	reads = []
	for i in range(6):
		c, p = sim.unique_locus()
		reads.append(sim.genomic_read(c, p, "+-"[i % 2], n_offsets=tuple(range(10, 100, 11))))
	reads = ReadSet(reads)
	run = run_star(genome_plain, reads, ATTRS + ["--seedSplitMin", "11"])
	assert len(run.alignments()) == 0
	assert ha.unmapped_counts(run.log_final()) == (0, 0, 6, 0)
	aln = run_star(genome_plain, reads, ATTRS + ["--seedSplitMin", "10"]).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		assert record_truth(rec) == reads.by_name[rec.qname].truth()
		assert (rec.tag("AS"), rec.tag("nM"), rec.tag("NM")) == (89, 0, 9)


@pytest.mark.parametrize("args,expected", [
	([], lambda n: n + "_x"),
	(["--readNameSeparator", "_"], lambda n: n),
	(["--readNameSeparator", "none"], lambda n: n + "_x/1"),
	(["--readNameSeparator", "x", "/"], lambda n: n + "_"),
	(["--outSAMreadID", "Number"], None),
])
def test_read_name_trimming(genome_plain, sim, tmp_path, args, expected):
	"""Read names 'NAME_x/1 comment': the comment after the space is always
	dropped; the rest is cut at the first character listed in
	--readNameSeparator (default '/'; 'none' cuts nothing; several characters
	may be given). --outSAMreadID Number replaces the name with the read
	number (index) in the file (parametersDefault); that it is 1-based pins
	the current output."""
	reads = [sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2], name="rd%d" % i) for i in range(6)]
	fq = "".join("@%s_x/1 1:N:0:ACGT\n%s\n+\n%s\n" % (r.name, r.seq, r.qual) for r in reads)
	aln = run_star(genome_plain, [_write(tmp_path / "names.fq", fq)], args).alignments()
	assert len(aln) == len(reads)
	for i, r in enumerate(reads):
		want = str(i + 1) if expected is None else expected(r.name)
		recs = [x for x in aln if x.qname == want]
		assert len(recs) == 1, (want, [x.qname for x in aln])
		assert record_truth(recs[0]) == r.truth()


@pytest.mark.parametrize("fasta,gz", [(True, False), (False, True)])
def test_paired_end_name_from_mate1(genome_plain, sim, tmp_path, fasta, gz):
	"""Paired FASTA, and gzipped paired FASTQ through --readFilesCommand
	gzip -cd: the QNAME of both records is mate 1's name (processChunks
	copies mate 1's read ID to both mates; mate 2's name line is ignored),
	both mates at their true alignment."""
	pairs = [sim.genomic_pair(*sim.unique_locus(400), insert=300, strand="+-"[i % 2], name="pair%d" % i)
		for i in range(8)]
	files = []
	for mate in (1, 2):
		lines = []
		for p in pairs:
			r = p.r1 if mate == 1 else p.r2
			name = p.name if mate == 1 else "other%s" % p.name[4:]
			lines.append(">%s\n%s\n" % (name, r.seq) if fasta else "@%s\n%s\n+\n%s\n" % (name, r.seq, r.qual))
		ext = (".fa" if fasta else ".fq") + (".gz" if gz else "")
		files.append(_write(tmp_path / ("m%d%s" % (mate, ext)), "".join(lines), gz))
	args = ["--readFilesCommand", "gzip", "-cd"] if gz else []
	aln = run_star(genome_plain, files, args).alignments()
	by_name = aln.by_name()
	assert sorted(by_name) == sorted(p.name for p in pairs)
	for p in pairs:
		recs = {r.mate: r for r in by_name[p.name]}
		assert record_truth(recs[1]) == p.r1.truth() and record_truth(recs[2]) == p.r2.truth()
		assert recs[1].qual == ("*" if fasta else p.r1.sam_qual())


def test_read_map_number(genome_plain, sim):
	"""--readMapNumber 7 maps only the first 7 reads of the file: the others
	have no record and Log.final.out counts 7 input reads."""
	reads = ReadSet([sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2]) for i in range(20)])
	run = run_star(genome_plain, reads, "--readMapNumber 7")
	aln = run.alignments()
	assert sorted(aln.by_name()) == sorted(r.name for r in reads.records[:7])
	for rec in aln:
		assert record_truth(rec) == reads.by_name[rec.qname].truth()
	assert run.log_final()["Number of input reads"] == 7


@pytest.mark.parametrize("add", [9, -9])
def test_quality_conversion_add(genome_plain, sim, add):
	"""--outQSconversionAdd adds to every quality character
	(parametersDefault); qualities here span the whole range, and the clamp
	of the result to 33..126 (readLoad.cpp) pins the current output."""
	quals = "".join(chr(33 + (i * 7) % 94) for i in range(100))
	reads = []
	for i in range(6):
		r = sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2])
		r.qual = quals[i:] + quals[:i]
		reads.append(r)
	reads = ReadSet(reads)
	aln = run_star(genome_plain, reads, ["--outQSconversionAdd", str(add)]).alignments()
	assert len(aln) == len(reads)
	for rec in aln:
		truth = reads.by_name[rec.qname]
		assert record_truth(rec) == truth.truth()
		want = "".join(chr(min(126, max(33, ord(q) + add))) for q in truth.sam_qual())
		assert rec.qual == want


def test_quality_length_mismatch_is_fatal(genome_plain, sim, tmp_path):
	"""A FASTQ record whose quality line is shorter than its sequence stops
	STAR with an error naming the problem (readLoad.cpp)."""
	r = sim.genomic_read(*sim.unique_locus(), name="bad")
	fq = "@%s\n%s\n+\n%s\n" % (r.name, r.seq, r.qual[:-3])
	run = run_star(genome_plain, [_write(tmp_path / "bad.fq", fq)], check=False)
	assert run.returncode != 0
	assert "quality string length is not equal to sequence length" in run.stderr + run.log_out()


def _soft_clip_truth(truth, left, right):
	"""Truth with `left`/`right` bases of an unspliced 'nM' read soft-clipped."""
	n = sh.cigar_query_len(truth.cigar)
	cig = sh.ops_to_cigar([(left, "S"), (n - left - right, "M"), (right, "S")])
	return (truth.contig, truth.pos + left, cig, truth.strand)


def test_clip_bases_per_mate(genome_plain, sim):
	"""--clip5pNbases 3 5 --clip3pNbases 2 4 clip each mate separately: mate 1
	loses 3 bases at its 5' end and 2 at its 3' end, mate 2 loses 5 and 4. A
	'+' mate keeps its 5' clip on the left of the CIGAR, a '-' mate on the
	right."""
	pairs = ReadSet([sim.genomic_pair(*sim.unique_locus(400), insert=300, strand="+-"[i % 2])
		for i in range(10)])
	aln = run_star(genome_plain, pairs, "--clip5pNbases 3 5 --clip3pNbases 2 4").alignments()
	by_name = aln.by_name()
	assert sorted(by_name) == sorted(pairs.by_name)
	for name, recs in by_name.items():
		p = pairs.by_name[name]
		for rec in recs:
			truth, c5, c3 = (p.r1, 3, 2) if rec.mate == 1 else (p.r2, 5, 4)
			left, right = (c5, c3) if truth.strand == "+" else (c3, c5)
			assert record_truth(rec) == _soft_clip_truth(truth, left, right)
			assert rec.is_proper_pair


def _adapter_reads(sim, n, insert=60, adapter=sh.ILLUMINA_ADAPTER, mm=()):
	"""Reads of `insert` genomic bases followed by the adapter (with
	substitutions at adapter offsets `mm`) and its first 6 bases, at loci
	where an unclipped adapter would not extend the alignment."""
	ad = list(adapter)
	for i in mm:
		ad[i] = ha.mismatched(ad[i])
	tail = ("".join(ad) + adapter[:6])[:100 - insert]
	reads = []
	while len(reads) < n:
		c, p = sim.unique_locus()
		strand = "+-"[len(reads) % 2]
		if ha.extension_past(sim.g.contigs[c], p, p + insert - 1, strand, tail) == 0:
			reads.append(sim.genomic_read(c, p, strand, length=insert, tail=tail))
	return ReadSet(reads)


def test_clip3p_after_adapter_bases(genome_plain, sim):
	"""--clip3pAfterAdapterNbases 5 clips 5 more 3' bases after the adapter:
	60 genomic bases + 40 adapter bases align with 55 bases, 55M45S at the
	read start for '+' reads and 45S55M five bases further right for '-'
	reads."""
	reads = _adapter_reads(sim, 10)
	args = ["--clip3pAdapterSeq", sh.ILLUMINA_ADAPTER, "--clip3pAfterAdapterNbases", "5"]
	aln = run_star(genome_plain, reads, args).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		truth = reads.by_name[rec.qname]
		if truth.strand == "+":
			assert record_truth(rec) == (truth.contig, truth.pos, "55M45S", "+")
		else:
			assert record_truth(rec) == (truth.contig, truth.pos + 5, "45S55M", "-")


def test_clip3p_adapter_mismatch_proportion(genome_plain, sim):
	"""The adapter copy in each read has 3 substitutions over its 34 bases.
	localSearch accepts an adapter position when mismatches/matches <=
	--clip3pAdapterMMp: 3/31 = 0.097 passes the default 0.1, so 40 bases are
	clipped and the 60 genomic bases align (60M40S / 40S60M; score 59 >=
	int(0.66 * 59)). With 0.05 only the trailing exact 6-base adapter prefix
	qualifies, the read keeps 94 bases, and 60 aligned bases (score 59) are
	below int(0.66 * 93) = 61: unmapped as too short."""
	reads = _adapter_reads(sim, 10, mm=(5, 15, 25))
	aln = run_star(genome_plain, reads, ["--clip3pAdapterSeq", sh.ILLUMINA_ADAPTER]).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		assert record_truth(rec) == reads.by_name[rec.qname].truth()
	run = run_star(genome_plain, reads, ["--clip3pAdapterSeq", sh.ILLUMINA_ADAPTER, "--clip3pAdapterMMp", "0.05"])
	assert len(run.alignments()) == 0
	assert ha.unmapped_counts(run.log_final()) == (0, 10, 0, 0)


AD2 = "GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT"


def test_clip3p_adapter_per_mate(genome_plain, sim):
	"""Pairs from 40-base fragments read 60 bases into a different adapter on
	each mate. --clip3pAdapterSeq AD1 AD2 clips both, so the pair keeps 40 +
	40 bases and maps as a proper pair (40M60S / 60S40M, both at the fragment
	start); its score 80 - 1 = 79 passes int(0.66 * 80) = 52. Giving AD1 for
	both mates leaves mate 2 unclipped: 140 bases need int(0.66 * 140) = 92,
	so the pairs are unmapped as too short."""
	pairs = []
	for i in range(8):
		c, p = sim.unique_locus(100)
		name = "short%d" % i
		r1 = sim.genomic_read(c, p, "+-"[i % 2], length=40, tail=(sh.ILLUMINA_ADAPTER * 2)[:60], name=name)
		r2 = sim.genomic_read(c, p, "-+"[i % 2], length=40, tail=(AD2 * 2)[:60], name=name)
		r1.mate, r2.mate = 1, 2
		pairs.append(SimPair(name, r1, r2))
	pairs = ReadSet(pairs)
	aln = run_star(genome_plain, pairs, ["--clip3pAdapterSeq", sh.ILLUMINA_ADAPTER, AD2, "--clip3pAdapterMMp", "0.1", "0.1"]).alignments()
	by_name = aln.by_name()
	assert sorted(by_name) == sorted(pairs.by_name)
	for name, recs in by_name.items():
		p = pairs.by_name[name]
		for rec in recs:
			assert record_truth(rec) == (p.r1 if rec.mate == 1 else p.r2).truth()
			assert rec.is_proper_pair and abs(rec.tlen) == 40
	run = run_star(genome_plain, pairs, ["--clip3pAdapterSeq", sh.ILLUMINA_ADAPTER, sh.ILLUMINA_ADAPTER, "--clip3pAdapterMMp", "0.1", "0.1"])
	assert len(run.alignments()) == 0
	assert ha.unmapped_counts(run.log_final()) == (0, 8, 0, 0)


def test_read_files_manifest_and_prefix(genome_plain, sim, tmp_path):
	"""--readFilesManifest lists one file per line with its read group; a
	read-group field without 'ID:' gets 'ID:' added, one starting with 'ID:'
	is copied verbatim to the @RG line (parametersDefault). --readFilesPrefix
	is prepended to every file name, from the manifest or from
	--readFilesIn. The files are gzipped and read through --readFilesCommand
	gzip -cd, one process per file. Every read aligns at its origin and,
	with RG requested, carries its file's read group."""
	a = ReadSet([sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2], name="ga%d" % i) for i in range(5)])
	b = ReadSet([sim.genomic_read(*sim.unique_locus(), strand="+-"[i % 2], name="gb%d" % i) for i in range(5)])
	fa = a.write(tmp_path, prefix="fa", gz=True)[0].name
	fb = b.write(tmp_path, prefix="fb", gz=True)[0].name
	manifest = tmp_path / "manifest.tsv"
	manifest.write_text("%s\t-\tgrpA\n%s\t-\tID:grpB\tSM:sb\n" % (fa, fb))
	prefix = ["--readFilesPrefix", str(tmp_path) + "/", "--readFilesCommand", "gzip", "-cd"]
	aln = run_star(genome_plain, None, prefix + ["--readFilesManifest", str(manifest),
		"--outSAMattributes", "NH", "HI", "RG"]).alignments()
	assert aln.header_lines("@RG") == ["@RG\tID:grpA", "@RG\tID:grpB\tSM:sb"]
	reads = a + b
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		assert record_truth(rec) == reads.by_name[rec.qname].truth()
		assert rec.tag("RG") == ("grpA" if rec.qname in a.by_name else "grpB")
	aln = run_star(genome_plain, None, prefix + ["--readFilesIn", fa + "," + fb]).alignments()
	assert sorted(aln.by_name()) == sorted(reads.by_name)
	for rec in aln:
		assert record_truth(rec) == reads.by_name[rec.qname].truth()
