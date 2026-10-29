"""Helpers for the area-B tests (SAM/BAM output): a small second synthetic
genome with 2- to 6-copy segments and hidden transcripts with mixed-strand or
non-canonical introns, expected alignment scores, unmappable mates, FASTQ
files with hand-written name lines, and header comparisons."""

import gzip
import math
import re
import struct

import star_helpers as sh
from star_helpers import SimPair, SimRead, cigar_ops, motif_code

# The tiny genome. Segments of 150 bp give 100-bp reads with exactly that many
# loci, so one run covers NH 1-6 and the whole MAPQ ladder. The hidden
# transcripts (not in any GTF) carry: Tmix, a GT/AG and a CT/AC intron
# (inconsistent strands) with blocks of 34/32/34 bases in the read so that
# no one-junction alignment reaches the 0.66 score filter; Tnc, a
# non-canonical CATG intron; Tgtag and Tctac, one canonical intron each.
TINY_SEED = 424242
TINY_CONTIGS = (("chrA", 16000), ("chrB", 12000))
TINY_GENES = (
	dict(gene_id="X01", name="Mixed", contig="chrA", strand="+", transcripts={},
		hidden={"Tmix": [(1001, 1100), (1401, 1432), (1801, 1900)],
			"Tnc": [(3001, 3100), (3501, 3600)],
			"Tgtag": [(5001, 5100), (5501, 5600)],
			"Tctac": [(7001, 7100), (7501, 7600)]},
		motifs={(1101, 1400): "GTAG", (1433, 1800): "CTAC", (3101, 3500): "CATG",
			(5101, 5500): "GTAG", (7101, 7500): "CTAC"}),
)
TINY_SEGMENTS = (
	("dup2", 150, [("chrA", 9001, "+"), ("chrB", 1001, "-")], ()),
	("dup3", 150, [("chrA", 9501, "+"), ("chrB", 1501, "+"), ("chrB", 2001, "-")], ()),
	("dup4", 150, [("chrA", 10001, "+"), ("chrA", 10501, "-"), ("chrB", 3001, "+"), ("chrB", 3501, "+")], ()),
	("dup5", 150, [("chrA", 11001, "+"), ("chrA", 11501, "+"), ("chrB", 4001, "-"), ("chrB", 4501, "+"),
		("chrB", 5001, "+")], ()),
	("dup6", 150, [("chrA", 12001, "+"), ("chrA", 12501, "+"), ("chrA", 13001, "+"), ("chrB", 6001, "+"),
		("chrB", 6501, "-"), ("chrB", 7001, "+")], ()),
)

_TINY = {}


def tiny_synth(tmp_path_factory):
	"""The tiny SyntheticGenome, built and written once per process so that
	genome_factory (which caches by object) indexes it once per session."""
	if "synth" not in _TINY:
		g = sh.make_synthetic_genome(seed=TINY_SEED, contig_specs=TINY_CONTIGS, gene_specs=TINY_GENES,
			segment_specs=TINY_SEGMENTS)
		g.write(tmp_path_factory.mktemp("tiny_synth"))
		_TINY["synth"] = g
	return _TINY["synth"]


# STAR's default scores (source/parametersDefault).
_MOTIF_PENALTY = {0: -8, 3: -4, 4: -4, 5: -8, 6: -8}  # scoreGapNoncan, scoreGapGCAG, scoreGapATAC


def expected_score(synth, mates, annotated=frozenset()):
	"""AS of an alignment from its truth, with STAR's default scores worked
	out by hand: +1 per matching and -1 per mismatching aligned base;
	scoreDelOpen + scoreDelBase*n (-2-2n) per deletion and the same for
	insertions; per intron +sjdbScore (2) if annotated, else scoreGap (0)
	plus the motif penalty (non-canonical -8, GC/AG -4, AT/AC -8); plus
	ceil(log2(span)*scoreGenomicLengthLog2scale - 0.5) with span the bases
	from the leftmost to the rightmost aligned base of all mates; floored at
	0. `mates` is one SimRead or the SimReads of a pair; `annotated` holds
	(contig, start, end) of annotated introns."""
	if isinstance(mates, SimRead):
		mates = [mates]
	score = 0
	for r in mates:
		g = synth.contigs[r.contig]
		seq = r.sam_seq()
		rp, qp = r.pos, 0
		for n, op in cigar_ops(r.cigar):
			if op == "M":
				score += sum(1 if g[rp - 1 + i] == seq[qp + i] else -1 for i in range(n))
				rp += n
				qp += n
			elif op == "S":
				qp += n
			elif op == "I":
				score += -2 - 2 * n
				qp += n
			elif op == "D":
				score += -2 - 2 * n
				rp += n
			elif op == "N":
				if (r.contig, rp, rp + n - 1) in annotated:
					score += 2
				else:
					score += _MOTIF_PENALTY.get(motif_code(g[rp - 1:rp + 1], g[rp + n - 3:rp + n - 1]), 0)
				rp += n
	span = max(r.end for r in mates) - min(r.pos for r in mates) + 1
	score += int(math.ceil(math.log2(span) * -0.25 - 0.5))
	return max(0, score)


def n_read(name, length=100, qual="I"):
	"""A read of N bases: no seeds, so it never aligns."""
	return SimRead(name, "N" * length, qual * length, kind="unmappable")


def half_pair(name, mapped, mapped_mate=1):
	"""Pair whose mate `mapped_mate` is the SimRead `mapped` and whose other
	mate is all N."""
	other = n_read(name)
	r1, r2 = (mapped, other) if mapped_mate == 1 else (other, mapped)
	r1.mate, r2.mate = 1, 2
	r1.name = r2.name = name
	return SimPair(name, r1, r2, kind="half")


def unmappable_pair(name):
	r1, r2 = n_read(name), n_read(name)
	r1.mate, r2.mate = 1, 2
	return SimPair(name, r1, r2, kind="unmappable")


def write_fastq(path, entries):
	"""Write FASTQ from (name line without '@', seq, qual) tuples."""
	with open(str(path), "w") as f:
		for name, seq, qual in entries:
			f.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))
	return path


def stable_header(header):
	"""Header lines without the two that hold the command line (@PG ID:STAR
	CL:..., @CO user command line), which differ between runs."""
	return [h for h in header if not h.startswith("@PG\tID:STAR\t") and not h.startswith("@CO\tuser command line")]


def max_sort_memory(run):
	"""'Max memory needed for sorting' from Log.out."""
	return int(re.search(r"Max memory needed for sorting = (\d+)", run.log_out()).group(1))


def contiguous_by_name(records):
	"""True if every QNAME's records form one contiguous run."""
	seen = set()
	prev = None
	for r in records:
		if r.qname != prev:
			if r.qname in seen:
				return False
			seen.add(r.qname)
			prev = r.qname
	return True


def mate_truth(pair, mate):
	return pair.r1 if mate == 1 else pair.r2


def assert_fastx_order(names, expected, threads=None):
	"""Names in Unmapped.out.mate1/2 against the unmapped reads in input
	order. Each mapping thread appends its reads to the file when it
	finishes (ReadAlignChunk_processChunks.cpp), so the input order is fixed
	only with one thread; with more, only the multiset is compared.
	`threads` defaults to the session's --runThreadN."""
	threads = sh.default_threads() if threads is None else threads
	if threads == 1:
		assert list(names) == list(expected)
	else:
		assert sorted(names) == sorted(expected)


# The empty BGZF block that ends a BAM file (SAM specification 4.1.2).
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def reg2bin(beg, end):
	"""BAM bin of the 0-based half-open interval [beg, end), copied from the
	SAM specification (5.3)."""
	end -= 1
	if beg >> 14 == end >> 14:
		return ((1 << 15) - 1) // 7 + (beg >> 14)
	if beg >> 17 == end >> 17:
		return ((1 << 12) - 1) // 7 + (beg >> 17)
	if beg >> 20 == end >> 20:
		return ((1 << 9) - 1) // 7 + (beg >> 20)
	if beg >> 23 == end >> 23:
		return ((1 << 6) - 1) // 7 + (beg >> 23)
	if beg >> 26 == end >> 26:
		return ((1 << 3) - 1) // 7 + (beg >> 26)
	return 0


def bam_bins(raw):
	"""(qname, pos0, bin, reference length of the CIGAR) of every record of a
	BAM file's bytes; star_helpers.read_bam does not keep the bin field."""
	data = gzip.decompress(raw)
	(l_text,) = struct.unpack_from("<i", data, 4)
	off = 8 + l_text
	(n_ref,) = struct.unpack_from("<i", data, off)
	off += 4
	for _ in range(n_ref):
		(l_name,) = struct.unpack_from("<i", data, off)
		off += 8 + l_name
	out = []
	while off < len(data):
		(block, _ref, pos, l_qname, _mapq, bin_, n_cig) = struct.unpack_from("<iiiBBHH", data, off)
		p = off + 36
		qname = data[p:p + l_qname - 1].decode()
		cig = struct.unpack_from("<%dI" % n_cig, data, p + l_qname)
		ref_len = sum(c >> 4 for c in cig if "MIDNSHP=X"[c & 15] in "MDN=X")
		out.append((qname, pos, bin_, ref_len))
		off += 4 + block
	return out


_INT_FORMATS = {"c": "<b", "C": "<B", "s": "<h", "S": "<H", "i": "<i", "I": "<I"}
_AUX_SIZES = {"A": 1, "c": 1, "C": 1, "s": 2, "S": 2, "i": 4, "I": 4, "f": 4}


def bam_aux_types(raw):
	"""[(qname, [(tag, type code, value or None)])] for every record of a
	BAM file's bytes, with the type codes as stored (SAM specification
	4.2.4); star_helpers.read_bam reports every integer type as 'i'. Values
	are decoded for single integers only; B arrays get 'B' + subtype."""
	data = gzip.decompress(raw)
	(l_text,) = struct.unpack_from("<i", data, 4)
	off = 8 + l_text
	(n_ref,) = struct.unpack_from("<i", data, off)
	off += 4
	for _ in range(n_ref):
		(l_name,) = struct.unpack_from("<i", data, off)
		off += 8 + l_name
	out = []
	while off < len(data):
		(block, _ref, _pos, l_qname, _mapq, _bin, n_cig, _flag, l_seq) = struct.unpack_from("<iiiBBHHHI", data, off)
		p = off + 36
		qname = data[p:p + l_qname - 1].decode()
		a = p + l_qname + 4 * n_cig + (l_seq + 1) // 2 + l_seq
		end = off + 4 + block
		tags = []
		while a < end:
			tag, typ = data[a:a + 2].decode(), chr(data[a + 2])
			a += 3
			if typ in _INT_FORMATS:
				(val,) = struct.unpack_from(_INT_FORMATS[typ], data, a)
				tags.append((tag, typ, val))
				a += _AUX_SIZES[typ]
			elif typ in _AUX_SIZES:
				tags.append((tag, typ, None))
				a += _AUX_SIZES[typ]
			elif typ in "ZH":
				tags.append((tag, typ, None))
				a = data.index(b"\0", a) + 1
			elif typ == "B":
				sub = chr(data[a])
				(n,) = struct.unpack_from("<i", data, a + 1)
				tags.append((tag, "B" + sub, None))
				a += 5 + n * _AUX_SIZES[sub]
			else:
				raise ValueError("unknown BAM aux type %r in %s" % (typ, qname))
		out.append((qname, tags))
		off = end
	return out


def smallest_uint_type(value):
	"""Smallest unsigned BAM integer type code (C, S, I) that holds the
	non-negative `value`."""
	assert value >= 0
	return "C" if value <= 0xFF else "S" if value <= 0xFFFF else "I"
