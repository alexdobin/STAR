"""Checks of the harness itself (synthetic genome, simulator truth, parsers),
so a failure elsewhere points at STAR rather than at star_helpers."""

import shutil
import struct
import subprocess
import zlib

import pytest

import star_helpers as sh
from star_helpers import ReadSet, run_star


def _aligned_from_truth(synth, read):
	"""Read bases implied by the truth CIGAR (soft clips and insertions taken
	from the read), with a list of genomic positions where read != genome."""
	g = synth.contigs[read.contig]
	seq = read.sam_seq()
	out, diffs = [], []
	rp, qp = read.pos, 0
	for n, op in sh.cigar_ops(read.cigar):
		if op == "M":
			for i in range(n):
				if seq[qp + i] != g[rp - 1 + i]:
					diffs.append(rp + i)
			out.append(seq[qp:qp + n])
			rp += n
			qp += n
		elif op in "IS":
			out.append(seq[qp:qp + n])
			qp += n
		elif op in "DN":
			rp += n
	return "".join(out), diffs


def test_synthetic_genome_layout(synth_genome):
	"""The default genome is deterministic, 110 kb with SAindexNbases 7, and
	every annotated intron carries the motif its gene spec asks for."""
	again = sh.make_synthetic_genome()
	assert again.contigs == synth_genome.contigs
	assert synth_genome.length == 110000 and synth_genome.sa_index_nbases() == 7
	motifs = {(j.contig, j.start, j.end): (j.strand, j.motif) for j in synth_genome.junctions()}
	assert motifs[("chr1", 1301, 2300)] == ("+", 1)
	assert motifs[("chr1", 6401, 7000)] == ("-", 2)
	assert motifs[("chr1", 11401, 11800)] == ("+", 3)
	assert motifs[("chr2", 2301, 3000)] == ("-", 4)
	assert motifs[("chr1", 15301, 15900)] == ("+", 5)
	assert motifs[("chr3", 12501, 13000)] == ("-", 6)
	hidden = {(j.contig, j.start, j.end): j for j in synth_genome.junctions(annotated_only=False)}
	assert hidden[("chr3", 15201, 15800)].motif == 0
	assert not hidden[("chr2", 30401, 31800)].annotated
	for seg in synth_genome.segments.values():
		seqs = [synth_genome.seq(c.contig, c.start, c.end) for c in seg.copies]
		seqs = [s if c.strand == "+" else sh.revcomp(s) for s, c in zip(seqs, seg.copies)]
		if seg.name != "dupmm":
			assert len(set(seqs)) == 1
		else:
			assert sh.hamming(seqs[0], seqs[1]) == 10


def test_simulated_truth_is_consistent(sim):
	"""For spliced, indel, clipped and multimapper reads, the bases implied by
	the truth CIGAR at the truth position reproduce the read, and the read
	differs from the genome exactly at the recorded mismatch positions."""
	reads = list(sim.transcript_reads(50))
	seg = sim.segment_read("dup2", 700, strand="-")
	contig, pos = sim.unique_locus(130)
	reads += [sim.genomic_read(contig, pos, "+", mismatches=(12, 40), deletion=(60, 4)),
		sim.genomic_read(contig, pos, "-", insertion=(33, 3), mismatches=(80,)),
		sim.genomic_read(contig, pos, "-", length=70, tail="A" * 30),
		sim.transcript_read("T08", 280, mismatches=(15,), sense=False), seg]
	reads += list(sim.transcript_pairs(10).mates(2))
	for r in reads:
		aligned, diffs = _aligned_from_truth(sim.g, r)
		assert aligned == r.sam_seq()
		assert diffs == list(r.mismatches), r
		assert sh.cigar_query_len(r.cigar) == len(r.seq) == len(r.qual)
	assert len(seg.extra["loci"]) == 2
	for c, p, strand in seg.extra["loci"]:
		s = sim.g.seq(c, p, p + 99)
		assert (s if strand == "+" else sh.revcomp(s)) == seg.seq


def test_sam_line_round_trip():
	"""parse_sam_line and SamRecord.to_line are inverse for typical STAR
	lines, including B arrays and unmapped records."""
	lines = [
		"r1\t99\tchr1\t1001\t255\t50M300N50M\t=\t1141\t240\tACGT\tIIII\tNH:i:1\tHI:i:1\tAS:i:198\tjM:B:c,21\tjI:B:i,1051,1350\tMC:Z:100M",
		"r2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tNH:i:0\tHI:i:0\tAS:i:11\tnM:i:0\tuT:A:1",
	]
	for line in lines:
		rec = sh.parse_sam_line(line)
		assert rec.to_line() == line
	rec = sh.parse_sam_line(lines[0])
	assert rec.introns == [(1051, 1350)] and rec.end == 1400 and rec.mate == 1


def test_parsers_on_handwritten_files(tmp_path):
	"""Log.final.out, SJ.out.tab, ReadsPerGene.out.tab and MatrixMarket
	parsers on small handwritten files."""
	(tmp_path / "Log.final.out").write_text(
		"                          Number of input reads |\t150\n"
		"                        Uniquely mapped reads % |\t66.67%\n"
		"                          Average mapped length |\t99.50\n"
		"                                 UNIQUE READS:\n")
	assert sh.parse_log_final(tmp_path / "Log.final.out") == {
		"Number of input reads": 150, "Uniquely mapped reads %": 66.67, "Average mapped length": 99.5}
	(tmp_path / "SJ.out.tab").write_text("chr1\t1301\t2300\t1\t1\t1\t4\t0\t48\n")
	assert sh.parse_sj_out(tmp_path / "SJ.out.tab") == [sh.SJRecord("chr1", 1301, 2300, 1, 1, 1, 4, 0, 48)]
	(tmp_path / "rpg").write_text("N_unmapped\t2\t2\t2\nG01\t5\t5\t0\n")
	assert sh.parse_reads_per_gene(tmp_path / "rpg") == {"N_unmapped": (2, 2, 2), "G01": (5, 5, 0)}
	(tmp_path / "m.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n%\n3 2 2\n1 1 4\n3 2 1\n")
	assert sh.parse_mtx(tmp_path / "m.mtx") == (3, 2, {(1, 1): 4, (3, 2): 1})


def _bgzf_block(data):
	"""One BGZF member (SAM/BAM specification 4.1) holding `data`."""
	c = zlib.compressobj(6, zlib.DEFLATED, -15)
	payload = c.compress(data) + c.flush()
	header = struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, len(payload) + 25)
	return header + payload + struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))


BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def _bam_record(ref_id, pos0, mapq, bin_, flag, qname, cigar, seq, qual, nref, npos0, tlen, tags):
	"""BAM alignment record bytes (SAM/BAM specification 4.2), built
	independently of star_helpers.read_bam."""
	cig = [n << 4 | "MIDNSHP=X".index(op) for n, op in sh.cigar_ops(cigar)]
	l_seq = len(seq)
	nib = "=ACMGRSVTWYHKDBN"
	packed = bytearray()
	for i in range(0, l_seq, 2):
		packed.append(nib.index(seq[i]) << 4 | (nib.index(seq[i + 1]) if i + 1 < l_seq else 0))
	qb = bytes([255] * l_seq) if qual == "*" else bytes(ord(ch) - 33 for ch in qual)
	name = qname.encode() + b"\x00"
	body = (struct.pack("<iiBBHHHiiii", ref_id, pos0, len(name), mapq, bin_, len(cig), flag, l_seq, nref, npos0, tlen)
		+ name + struct.pack("<%dI" % len(cig), *cig) + bytes(packed) + qb + tags)
	return struct.pack("<i", len(body)) + body


def test_bam_decoder_on_handmade_bam():
	"""read_bam on a BAM written by hand from the SAM/BAM specification, so
	the decoder is checked without samtools or pysam: every tag type (A,
	c, C, s, S, i, I, f, Z, H and B arrays), an odd-length SEQ, QUAL absent
	(0xFF), CIGAR with S, I, D and N, an unmapped record whose mate is on
	another reference, the bin field, and a record split across two BGZF
	blocks. bam_types keeps the stored integer type codes, which
	record_multiset(bam_types=True) compares."""
	text = b"@HD\tVN:1.4\n@SQ\tSN:chrA\tLN:1000\n@SQ\tSN:chrB\tLN:2000\n"
	head = b"BAM\x01" + struct.pack("<i", len(text)) + text + struct.pack("<i", 2)
	for name, ln in ((b"chrA", 1000), (b"chrB", 2000)):
		head += struct.pack("<i", len(name) + 1) + name + b"\x00" + struct.pack("<i", ln)
	tags = (b"XAAx" + b"Xcc" + struct.pack("<b", -5) + b"XCC" + struct.pack("<B", 200)
		+ b"Xss" + struct.pack("<h", -300) + b"XSS" + struct.pack("<H", 60000)
		+ b"Xii" + struct.pack("<i", -70000) + b"XII" + struct.pack("<I", 3000000000)
		+ b"Xff" + struct.pack("<f", 1.5) + b"XZZhello\x00" + b"XHH1AE3\x00"
		+ b"BcBc" + struct.pack("<i2b", 2, -1, 2) + b"BSBS" + struct.pack("<i2H", 2, 1, 65535)
		+ b"BfBf" + struct.pack("<i2f", 2, 0.5, 2.0) + b"BIBI" + struct.pack("<iI", 1, 4000000000))
	seq1 = "NACGTACGTACGTACGTACGTACGTACGTAC"
	qual1 = "".join(chr(33 + (i % 41)) for i in range(len(seq1)))
	r1 = _bam_record(0, 100, 60, 4681, 99, "r1", "5S10M2I3M4D6M100N5M", seq1, qual1, 0, 300, 250, tags)
	r2 = _bam_record(-1, -1, 0, 4680, 0x45, "r2", "*", "ACGTN", "*", 1, 1499, 0, b"NHC\x00")
	r3 = _bam_record(1, 0, 3, 4681, 0x91, "r3", "4M", "TTTT", "IIII", 0, 9, 0, b"")
	data = head + r1 + r2 + r3
	split = len(head) + len(r1) + 7
	raw = _bgzf_block(data[:split]) + _bgzf_block(data[split:]) + BGZF_EOF
	bam = sh.read_bam(raw)
	assert bam.header == text.decode().splitlines()
	assert bam.references == [("chrA", 1000), ("chrB", 2000)]
	assert [r.to_line() for r in bam] == [
		"\t".join(["r1", "99", "chrA", "101", "60", "5S10M2I3M4D6M100N5M", "=", "301", "250", seq1, qual1,
			"XA:A:x", "Xc:i:-5", "XC:i:200", "Xs:i:-300", "XS:i:60000", "Xi:i:-70000", "XI:i:3000000000",
			"Xf:f:1.5", "XZ:Z:hello", "XH:H:1AE3", "Bc:B:c,-1,2", "BS:B:S,1,65535", "Bf:B:f,0.5,2",
			"BI:B:I,4000000000"]),
		"r2\t69\t*\t0\t0\t*\tchrB\t1500\t0\tACGTN\t*\tNH:i:0",
		"r3\t145\tchrB\t1\t3\t4M\tchrA\t10\t0\tTTTT\tIIII",
	]
	assert [r.bin for r in bam] == [4681, 4680, 4681]
	assert "".join(bam.records[0].bam_types.values()) == "AcCsSiIfZH" + "B:cB:SB:fB:I"
	assert bam.records[1].bam_types == {"NH": "C"} and bam.records[2].bam_types == {}
	retyped = sh.read_bam(_bgzf_block(head + r1 + r2.replace(b"NHC\x00", b"NHc\x00") + r3) + BGZF_EOF)
	assert sh.record_multiset(retyped) == sh.record_multiset(bam)
	assert sh.record_multiset(retyped, bam_types=True) != sh.record_multiset(bam, bam_types=True)
	r = bam.records[0]
	assert (r.end, r.introns, r.tag_types["Bc"], r.tags["Xf"]) == (100 + 10 + 3 + 4 + 6 + 100 + 5, [(124, 223)], "B:c", 1.5)
	assert sh.parse_sam_line(r.to_line()).bin is None


def test_solo_matrix_keys_and_out_names(tmp_path):
	"""read_solo_matrix keys rows by gene ID when IDs are unique and by
	(contig, start, end) for the SJ feature, whose features.tsv is
	SJ.out.tab (contig repeats); StarRun.solo* follow --soloOutFileNames."""
	d = tmp_path / "My" / "SJ" / "raw"
	d.mkdir(parents=True)
	(d / "f.tsv").write_text("chr1\t1301\t2300\t1\t1\t1\t4\t0\t48\nchr1\t1301\t3450\t1\t1\t1\t2\t0\t30\n")
	(d / "c.tsv").write_text("AAAA\nCCCC\n")
	(d / "m.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n%\n2 2 3\n1 1 4\n2 1 1\n2 2 2\n")
	g = tmp_path / "My" / "Gene" / "raw"
	g.mkdir(parents=True)
	(g / "f.tsv").write_text("G01\tGeneOne\tGene Expression\nG02\tGeneTwo\tGene Expression\n")
	(g / "c.tsv").write_text("AAAA\n")
	(g / "m.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n%\n2 1 1\n2 1 7\n")
	(tmp_path / "My" / "Barcodes.stats").write_text("    nNoAdapter    0\n    nNinBarcode    3\n")
	argv = ["STAR", "--soloOutFileNames", "My/", "f.tsv", "c.tsv", "m.mtx", "--soloType", "CB_UMI_Simple"]
	run = sh.StarRun(tmp_path, argv, subprocess.CompletedProcess(argv, 0, b"", b""), [], None)
	assert run.solo_out_names() == ["My/", "f.tsv", "c.tsv", "m.mtx"]
	sj = run.solo("SJ")
	assert sj.keys == [("chr1", 1301, 2300), ("chr1", 1301, 3450)]
	assert sj.counts == {(("chr1", 1301, 2300), "AAAA"): 4, (("chr1", 1301, 3450), "AAAA"): 1,
		(("chr1", 1301, 3450), "CCCC"): 2}
	assert run.solo("Gene").counts == {("G02", "AAAA"): 7}
	assert run.solo_barcodes_stats() == {"nNoAdapter": 0, "nNinBarcode": 3}
	default = sh.StarRun(tmp_path, ["STAR"], subprocess.CompletedProcess(["STAR"], 0, b"", b""), [], None)
	assert default.solo_dir("Gene") == tmp_path / "Solo.out" / "Gene"


def _bam_run(genome_sjdb, sim):
	reads = sim.transcript_reads(60) + ReadSet([sim.segment_read("dup2", 100), sim.random_read()])
	return run_star(genome_sjdb, reads, "--outSAMtype BAM Unsorted --outSAMunmapped Within "
		"--outSAMattributes NH HI AS nM NM MD jM jI")


def test_bam_decoder_matches_samtools(genome_sjdb, sim):
	"""read_bam gives the same lines as `samtools view -h` (skipped without
	samtools on PATH)."""
	samtools = shutil.which("samtools")
	if samtools is None:
		pytest.skip("samtools not on PATH")
	run = _bam_run(genome_sjdb, sim)
	bam = run.path("Aligned.out.bam")
	text = subprocess.run([samtools, "view", "-h", str(bam)], stdout=subprocess.PIPE, check=True).stdout.decode()
	ref = sh.parse_sam(text)
	ours = sh.read_bam(bam)
	assert [r.to_line() for r in ours] == [r.to_line() for r in ref]
	assert ours.references == ref.references
	assert [h for h in ours.header if not h.startswith("@PG\tID:samtools")] == \
		[h for h in ref.header if not h.startswith("@PG\tID:samtools")]


def test_bam_decoder_matches_pysam(genome_sjdb, sim):
	"""read_bam agrees with pysam on every field (skipped without pysam)."""
	pysam = pytest.importorskip("pysam")
	run = _bam_run(genome_sjdb, sim)
	bam = run.path("Aligned.out.bam")
	ours = sh.read_bam(bam)
	with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as f:
		theirs = [sh.parse_sam_line(a.to_string()) for a in f]
	assert [r.to_line() for r in ours] == [r.to_line() for r in theirs]
