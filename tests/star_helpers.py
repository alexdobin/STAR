"""Harness for the STAR unit tests.

Synthetic genomes with known annotation, a read simulator that records each
read's true alignment, a runner for the STAR binary, and parsers for STAR's
outputs (SAM, BAM, SJ.out.tab, ReadsPerGene.out.tab, Log.final.out, Solo.out,
Unmapped.out.mate*, Chimeric.out.junction). Only the standard library is used.

Coordinates are 1-based and inclusive throughout, as in SAM and SJ.out.tab.
"""

import gzip
import math
import os
import random
import re
import shlex
import struct
import subprocess
import tempfile
from collections import OrderedDict, namedtuple
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

GENOME_SEED = 20261024
DEFAULT_READ_LEN = 100
DEFAULT_SJDB_OVERHANG = DEFAULT_READ_LEN - 1
ILLUMINA_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
MAX_THREADS = 8
DEFAULT_TIMEOUT = 120

# --outSAMtype arguments for parametrizing a test over STAR's three writers
# (SAM text, unsorted BAM, coordinate-sorted BAM); run.alignments() reads
# whichever was written.
OUTPUT_FORMATS = OrderedDict([
	("sam", []),
	("bam", ["--outSAMtype", "BAM", "Unsorted"]),
	("sorted", ["--outSAMtype", "BAM", "SortedByCoordinate"]),
])

# --soloOutFileNames defaults (parametersDefault): output prefix, features,
# barcodes and matrix file names.
SOLO_OUT_FILE_NAMES = ("Solo.out/", "features.tsv", "barcodes.tsv", "matrix.mtx")

BUILD_HINT = (
	"STAR binary not found. Build it with `make -C source STAR` from the "
	"repository root, or point STAR_BIN at an existing binary."
)


# ---------------------------------------------------------------- sequences

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
	return seq.translate(_COMP)[::-1]


def hamming(a, b):
	return sum(1 for x, y in zip(a, b) if x != y)


# STAR's intron motif codes (SJ.out.tab column 5, jM attribute).
MOTIF_CODES = {
	("GT", "AG"): 1, ("CT", "AC"): 2,
	("GC", "AG"): 3, ("CT", "GC"): 4,
	("AT", "AC"): 5, ("GT", "AT"): 6,
}
# SJ.out.tab column 4: 0 undefined, 1 +, 2 -. Odd motif codes are on +.
MOTIF_STRAND = {0: 0, 1: 1, 2: 2, 3: 1, 4: 2, 5: 1, 6: 2}


def motif_code(left2, right2):
	"""STAR motif code (0-6) of an intron whose first two genomic bases are
	`left2` and last two are `right2`, both on the + strand."""
	return MOTIF_CODES.get((left2.upper(), right2.upper()), 0)


# -------------------------------------------------------------------- CIGAR

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_ops(cigar):
	"""'10M5N3M' -> [(10, 'M'), (5, 'N'), (3, 'M')]; '*' -> []."""
	if cigar in ("*", "", None):
		return []
	ops = [(int(n), op) for n, op in _CIGAR_RE.findall(cigar)]
	if "".join("%d%s" % o for o in ops) != cigar:
		raise ValueError("malformed CIGAR %r" % cigar)
	return ops


def ops_to_cigar(ops):
	"""Join (length, op) pairs, merging neighbours with the same op."""
	merged = []
	for n, op in ops:
		if n == 0:
			continue
		if merged and merged[-1][1] == op:
			merged[-1] = (merged[-1][0] + n, op)
		else:
			merged.append((n, op))
	return "".join("%d%s" % o for o in merged) or "*"


def cigar_ref_len(cigar):
	"""Reference bases spanned (M, D, N, =, X)."""
	return sum(n for n, op in cigar_ops(cigar) if op in "MDN=X")


def cigar_query_len(cigar):
	"""Read bases consumed (M, I, S, =, X)."""
	return sum(n for n, op in cigar_ops(cigar) if op in "MIS=X")


def cigar_introns(pos, cigar):
	"""(start, end) of every N operation, 1-based inclusive."""
	out = []
	p = pos
	for n, op in cigar_ops(cigar):
		if op == "N":
			out.append((p, p + n - 1))
		if op in "MDN=X":
			p += n
	return out


def cigar_blocks(pos, cigar):
	"""(start, end) of every aligned block (M, =, X), 1-based inclusive."""
	out = []
	p = pos
	for n, op in cigar_ops(cigar):
		if op in "M=X":
			out.append((p, p + n - 1))
		if op in "MDN=X":
			p += n
	return out


def cigar_softclips(cigar):
	"""(left, right) soft-clipped lengths."""
	ops = cigar_ops(cigar)
	left = ops[0][0] if ops and ops[0][1] == "S" else 0
	right = ops[-1][0] if len(ops) > 1 and ops[-1][1] == "S" else 0
	return left, right


# ------------------------------------------------------- synthetic genome

class Transcript(object):
	"""A transcript: exons are (start, end), 1-based inclusive, ascending."""

	def __init__(self, tx_id, gene_id, contig, strand, exons, annotated=True):
		self.tx_id = tx_id
		self.gene_id = gene_id
		self.contig = contig
		self.strand = strand
		self.exons = sorted(exons)
		self.annotated = annotated

	@property
	def start(self):
		return self.exons[0][0]

	@property
	def end(self):
		return self.exons[-1][1]

	@property
	def length(self):
		return sum(e - s + 1 for s, e in self.exons)

	@property
	def introns(self):
		return [(a[1] + 1, b[0] - 1) for a, b in zip(self.exons, self.exons[1:])]

	def genomic_positions(self):
		"""Genomic positions of the transcript's bases, ascending."""
		out = []
		for s, e in self.exons:
			out.extend(range(s, e + 1))
		return out

	def sense_positions(self):
		"""Genomic positions in transcript (5' to 3') order."""
		pos = self.genomic_positions()
		return pos if self.strand == "+" else pos[::-1]

	def __repr__(self):
		return "Transcript(%s %s:%s %s)" % (self.tx_id, self.contig, self.strand, self.exons)


class Gene(object):
	"""A gene: its transcripts (annotated and hidden), strand and biotype.
	annotated=False for genes left out of the GTF (G12)."""

	def __init__(self, gene_id, name, contig, strand, transcripts, biotype="protein_coding", annotated=True):
		self.gene_id = gene_id
		self.name = name
		self.contig = contig
		self.strand = strand
		self.transcripts = transcripts
		self.biotype = biotype
		self.annotated = annotated

	@property
	def start(self):
		return min(t.start for t in self.transcripts)

	@property
	def end(self):
		return max(t.end for t in self.transcripts)

	def exonic_positions(self):
		out = set()
		for t in self.transcripts:
			if t.annotated:
				out.update(t.genomic_positions())
		return out


Segment = namedtuple("Segment", "name length copies")
SegmentCopy = namedtuple("SegmentCopy", "contig start end strand")
Junction = namedtuple("Junction", "contig start end strand motif annotated")


# The default layout. Coordinates were picked by hand so tests can name them.
CONTIG_SPECS = (("chr1", 50000), ("chr2", 40000), ("chr3", 20000))

# Motifs are written in the transcript's sense: donor dinucleotide + acceptor
# dinucleotide. Introns not listed are GT/AG.
GENE_SPECS = (
	dict(gene_id="G01", name="GeneOne", contig="chr1", strand="+",
		transcripts={"T01a": [(1001, 1300), (2301, 2450), (3451, 3800)],
			"T01b": [(1001, 1300), (3451, 3800)]}),
	dict(gene_id="G02", name="GeneTwo", contig="chr1", strand="-",
		transcripts={"T02": [(6001, 6400), (7001, 7200), (8001, 8500)]}),
	dict(gene_id="G03", name="GeneThree", contig="chr1", strand="+",
		transcripts={"T03": [(11001, 11400), (11801, 12200)]},
		motifs={(11401, 11800): "GCAG"}),
	dict(gene_id="G04", name="GeneFour", contig="chr1", strand="+",
		transcripts={"T04": [(15001, 15300), (15901, 16300)]},
		motifs={(15301, 15900): "ATAC"}),
	dict(gene_id="G05", name="GeneFive", contig="chr1", strand="+",
		transcripts={"T05": [(20001, 20600), (21001, 21400)]}),
	dict(gene_id="G06", name="GeneSix", contig="chr1", strand="-",
		transcripts={"T06": [(21201, 21800), (22301, 22600)]}),
	dict(gene_id="G07", name="GeneSeven", contig="chr1", strand="-",
		transcripts={"T07": [(30001, 31000)]}),
	dict(gene_id="G08", name="GeneEight", contig="chr2", strand="-",
		transcripts={"T08": [(2001, 2300), (3001, 3150), (4001, 4200), (5001, 5400)]},
		motifs={(2301, 3000): "GCAG"}),
	dict(gene_id="G09", name="GeneNine", contig="chr2", strand="+",
		transcripts={"T09": [(12001, 12300), (17301, 17600)]}),
	dict(gene_id="G10", name="GeneTen", contig="chr2", strand="+",
		transcripts={"T10": [(30001, 30400), (30901, 31300), (31801, 32400)]},
		hidden={"T10skip": [(30001, 30400), (31801, 32400)]}),
	dict(gene_id="G11", name="GeneEleven", contig="chr3", strand="-",
		transcripts={"T11": [(12001, 12500), (13001, 13300)]},
		motifs={(12501, 13000): "ATAC"}),
	dict(gene_id="G12", name="GeneTwelve", contig="chr3", strand="+",
		transcripts={}, hidden={"T12nc": [(15001, 15200), (15801, 16000)]},
		motifs={(15201, 15800): "CATG"}),
)

# name, length, copies [(contig, start, strand)], substitutions in copies 2+
# (offsets within the segment, in the segment's own orientation).
SEGMENT_SPECS = (
	("dup2", 1500, [("chr1", 26001, "+"), ("chr2", 8001, "-")], ()),
	("dup3", 600, [("chr1", 33001, "+"), ("chr2", 20001, "+"), ("chr3", 2001, "+")], ()),
	("dupmm", 1000, [("chr1", 35001, "+"), ("chr3", 9001, "+")], tuple(range(50, 1000, 100))),
	("rep12", 300, [("chr1", 38001, "+"), ("chr1", 39001, "+"), ("chr1", 40001, "+"), ("chr1", 41001, "+"),
		("chr2", 23001, "+"), ("chr2", 24001, "+"), ("chr2", 25001, "+"), ("chr2", 26001, "+"),
		("chr3", 4001, "+"), ("chr3", 5001, "+"), ("chr3", 6001, "+"), ("chr3", 7001, "+")], ()),
)


class SyntheticGenome(object):
	"""Random chromosomes with a hand-placed annotation.

	Attributes: contigs (OrderedDict name -> sequence), genes (list of Gene,
	including the unannotated G12), transcripts (dict tx_id -> Transcript,
	including hidden transcripts that are not written to the GTF), segments
	(dict name -> Segment, repeated sequence with all copies listed),
	fasta_path, gtf_path (set by write()).
	"""

	def __init__(self, contigs, genes, segments, seed):
		self.contigs = contigs
		self.genes = genes
		self.segments = segments
		self.seed = seed
		self.transcripts = OrderedDict()
		for g in genes:
			for t in g.transcripts:
				self.transcripts[t.tx_id] = t
		self.fasta_path = None
		self.gtf_path = None

	@property
	def length(self):
		return sum(len(s) for s in self.contigs.values())

	def sa_index_nbases(self):
		"""min(14, floor(log2(GenomeLength)/2 - 1)), from STAR's manual."""
		return int(min(14, math.floor(math.log2(self.length) / 2 - 1)))

	def seq(self, contig, start, end):
		"""Forward-strand sequence of contig:start-end (1-based inclusive)."""
		return self.contigs[contig][start - 1:end]

	def gene(self, gene_id):
		for g in self.genes:
			if g.gene_id == gene_id:
				return g
		raise KeyError(gene_id)

	def annotated_genes(self):
		return [g for g in self.genes if g.annotated]

	def transcript_seq(self, tx_id):
		"""Transcript sequence in its sense (5' to 3') orientation."""
		t = self.transcripts[tx_id]
		s = "".join(self.seq(t.contig, a, b) for a, b in t.exons)
		return s if t.strand == "+" else revcomp(s)

	def junction(self, contig, start, end, strand, annotated):
		"""Junction for an intron, with its motif read from the sequence."""
		left = self.seq(contig, start, start + 1)
		right = self.seq(contig, end - 1, end)
		return Junction(contig, start, end, strand, motif_code(left, right), annotated)

	def junctions(self, annotated_only=True):
		"""Set of Junction(contig, start, end, strand, motif, annotated) for
		all introns of annotated (or all) transcripts."""
		out = set()
		annotated = set()
		for t in self.transcripts.values():
			if t.annotated:
				annotated.update((t.contig, s, e) for s, e in t.introns)
		for t in self.transcripts.values():
			if annotated_only and not t.annotated:
				continue
			for s, e in t.introns:
				out.add(self.junction(t.contig, s, e, t.strand, (t.contig, s, e) in annotated))
		return out

	def write(self, directory):
		"""Write genome.fa and annotation.gtf into `directory`; sets
		fasta_path and gtf_path. Returns self."""
		directory = Path(directory).resolve()
		directory.mkdir(parents=True, exist_ok=True)
		self.fasta_path = directory / "genome.fa"
		self.gtf_path = directory / "annotation.gtf"
		with open(self.fasta_path, "w") as f:
			for name, s in self.contigs.items():
				f.write(">%s synthetic\n" % name)
				for i in range(0, len(s), 60):
					f.write(s[i:i + 60] + "\n")
		with open(self.gtf_path, "w") as f:
			f.write(self.gtf_text())
		return self

	def gtf_text(self):
		lines = []
		for g in self.genes:
			if not g.annotated:
				continue
			txs = [t for t in g.transcripts if t.annotated]
			gattr = 'gene_id "%s"; gene_name "%s"; gene_biotype "%s";' % (g.gene_id, g.name, g.biotype)
			lines.append("\t".join([g.contig, "synthetic", "gene", str(g.start), str(g.end), ".", g.strand, ".", gattr]))
			for t in txs:
				tattr = 'gene_id "%s"; transcript_id "%s"; gene_name "%s"; gene_biotype "%s";' % (
					g.gene_id, t.tx_id, g.name, g.biotype)
				lines.append("\t".join([g.contig, "synthetic", "transcript", str(t.start), str(t.end), ".", t.strand, ".", tattr]))
				exons = t.exons if t.strand == "+" else t.exons[::-1]
				for i, (s, e) in enumerate(exons, 1):
					lines.append("\t".join([g.contig, "synthetic", "exon", str(s), str(e), ".", t.strand, ".",
						tattr + ' exon_number "%d";' % i]))
		return "\n".join(lines) + "\n"


def make_synthetic_genome(seed=GENOME_SEED, contig_specs=CONTIG_SPECS, gene_specs=GENE_SPECS,
		segment_specs=SEGMENT_SPECS):
	"""Build the default synthetic genome deterministically from `seed`."""
	rng = random.Random(seed)
	seqs = OrderedDict((name, [rng.choice("ACGT") for _ in range(n)]) for name, n in contig_specs)
	painted = {}

	def paint(contig, start, bases, what):
		for i, b in enumerate(bases):
			key = (contig, start + i)
			if key in painted and painted[key] != b:
				raise ValueError("layout conflict at %s:%d (%s)" % (contig, start + i, what))
			painted[key] = b
			seqs[contig][start + i - 1] = b

	segments = OrderedDict()
	for name, length, copies, subs in segment_specs:
		base = "".join(rng.choice("ACGT") for _ in range(length))
		seg_copies = []
		for k, (contig, start, strand) in enumerate(copies):
			s = list(base)
			if k > 0:
				for off in subs:
					s[off] = rng.choice([b for b in "ACGT" if b != s[off]])
			s = "".join(s)
			paint(contig, start, s if strand == "+" else revcomp(s), name)
			seg_copies.append(SegmentCopy(contig, start, start + length - 1, strand))
		segments[name] = Segment(name, length, seg_copies)

	genes = []
	for spec in gene_specs:
		txs = [Transcript(tid, spec["gene_id"], spec["contig"], spec["strand"], ex, True)
			for tid, ex in spec["transcripts"].items()]
		txs += [Transcript(tid, spec["gene_id"], spec["contig"], spec["strand"], ex, False)
			for tid, ex in spec.get("hidden", {}).items()]
		annotated = bool(spec["transcripts"])
		genes.append(Gene(spec["gene_id"], spec["name"], spec["contig"], spec["strand"], txs, annotated=annotated))
		motifs = spec.get("motifs", {})
		for t in txs:
			for s, e in t.introns:
				m = motifs.get((s, e))
				if m is None:
					# an exon-skipping intron takes the donor and acceptor of the
					# introns it spans
					donor = next((motifs[k] for k in motifs if (k[0] == s if t.strand == "+" else k[1] == e)), "GTAG")[:2]
					acc = next((motifs[k] for k in motifs if (k[1] == e if t.strand == "+" else k[0] == s)), "GTAG")[2:]
					m = donor + acc
				if t.strand == "+":
					left, right = m[:2], m[2:]
				else:
					left, right = revcomp(m[2:]), revcomp(m[:2])
				paint(t.contig, s, left, t.tx_id)
				paint(t.contig, e - 1, right, t.tx_id)

	# Make every junction position unambiguous: the last exon base differs from
	# the last intron base and the first exon base differs from the first intron
	# base, so the junction cannot slide (repeat length 0).
	for g in genes:
		for t in g.transcripts:
			c = seqs[t.contig]
			for s, e in t.introns:
				for exon_pos, intron_pos in ((s - 1, e), (e + 1, s)):
					forbidden = c[intron_pos - 1]
					if c[exon_pos - 1] == forbidden:
						if (t.contig, exon_pos) in painted:
							raise ValueError("cannot fix junction repeat at %s:%d" % (t.contig, exon_pos))
						c[exon_pos - 1] = next(b for b in "ACGT" if b != forbidden)
	contigs = OrderedDict((k, "".join(v)) for k, v in seqs.items())
	for g in genes:
		for t in g.transcripts:
			c = contigs[t.contig]
			for s, e in t.introns:
				assert c[s - 2] != c[e - 1] and c[e] != c[s - 1], (t.tx_id, s, e)
	return SyntheticGenome(contigs, genes, segments, seed)


# ------------------------------------------------------------ simulated reads

class SimRead(object):
	"""One simulated read and its true alignment.

	seq/qual are as sequenced (what goes into the FASTQ). The truth is in SAM
	terms: contig, strand ('+' when the read equals the forward genome),
	pos (leftmost aligned base, 1-based) and cigar (genome orientation, S for
	adapter or clipped bases, N for introns). contig is None for reads that
	should not align. mismatches lists the genomic positions of substitutions.
	extra holds kind-specific truth (tx_offset, loci, cb, umi, ...).
	"""

	__slots__ = ("name", "seq", "qual", "contig", "strand", "pos", "cigar", "mismatches",
		"kind", "tx", "gene", "mate", "extra")

	def __init__(self, name, seq, qual, contig=None, strand=None, pos=None, cigar=None,
			mismatches=(), kind="genomic", tx=None, gene=None, mate=None, extra=None):
		self.name = name
		self.seq = seq
		self.qual = qual
		self.contig = contig
		self.strand = strand
		self.pos = pos
		self.cigar = cigar
		self.mismatches = tuple(mismatches)
		self.kind = kind
		self.tx = tx
		self.gene = gene
		self.mate = mate
		self.extra = extra if extra is not None else {}

	@property
	def mappable(self):
		return self.contig is not None

	@property
	def end(self):
		return self.pos + cigar_ref_len(self.cigar) - 1

	@property
	def introns(self):
		return cigar_introns(self.pos, self.cigar) if self.cigar else []

	def truth(self):
		"""(contig, pos, cigar, strand), comparable with record_truth(rec)."""
		return (self.contig, self.pos, self.cigar, self.strand)

	def sam_seq(self):
		"""Sequence as it appears in SAM (reverse-complemented for '-')."""
		return self.seq if self.strand != "-" else revcomp(self.seq)

	def sam_qual(self):
		return self.qual if self.strand != "-" else self.qual[::-1]

	def __repr__(self):
		return "SimRead(%s %s:%s:%s %s %s)" % (self.name, self.contig, self.pos, self.strand, self.cigar, self.kind)


class SimPair(object):
	"""A read pair. r1 and r2 are in STAR's --readFilesIn order; for 10x
	reads r1 is the cDNA read and r2 the barcode read."""

	__slots__ = ("name", "r1", "r2", "insert", "kind", "extra")

	def __init__(self, name, r1, r2, insert=None, kind="pair", extra=None):
		self.name = name
		self.r1 = r1
		self.r2 = r2
		self.insert = insert
		self.kind = kind
		self.extra = extra if extra is not None else {}

	def __iter__(self):
		return iter((self.r1, self.r2))

	def __repr__(self):
		return "SimPair(%s %r %r)" % (self.name, self.r1, self.r2)


def _fastx_name(name, mate, style):
	"""Read name line: 'plain' name, 'slash' name/1, 'illumina' name 1:N:0:1
	(mate 1 for single-end reads)."""
	mate = mate or 1
	if style == "plain":
		return name
	if style == "slash":
		return "%s/%d" % (name, mate)
	if style == "illumina":
		return "%s %d:N:0:1" % (name, mate)
	raise ValueError("unknown name style %r" % style)


class ReadSet(object):
	"""Simulated reads: a list of SimRead (single-end) or SimPair (paired-end,
	10x). Iterable; `by_name[name]` gives the SimRead or SimPair."""

	def __init__(self, records, kind=None):
		self.records = list(records)
		if kind is None:
			kind = "pe" if self.records and isinstance(self.records[0], SimPair) else "se"
		self.kind = kind
		for r in self.records:
			if (kind == "se") != isinstance(r, SimRead):
				raise TypeError("ReadSet(kind=%r) cannot hold %r" % (kind, r))
		self.by_name = OrderedDict((r.name, r) for r in self.records)
		if len(self.by_name) != len(self.records):
			raise ValueError("duplicate read names in ReadSet")

	@property
	def paired(self):
		return self.kind != "se"

	def __iter__(self):
		return iter(self.records)

	def __len__(self):
		return len(self.records)

	def __add__(self, other):
		if self.kind != other.kind:
			raise TypeError("cannot add %s and %s read sets" % (self.kind, other.kind))
		return ReadSet(self.records + other.records, self.kind)

	def mates(self, i):
		"""SimReads of mate i (1 or 2); for single-end reads, mate 1 is all."""
		if not self.paired:
			if i != 1:
				raise ValueError("single-end ReadSet has only mate 1")
			return list(self.records)
		return [p.r1 if i == 1 else p.r2 for p in self.records]

	def write(self, directory, prefix="reads", fasta=False, gz=False, name_style="plain"):
		"""Write FASTQ (or FASTA) files and return their paths in
		--readFilesIn order (one for single-end, two otherwise)."""
		directory = Path(directory)
		directory.mkdir(parents=True, exist_ok=True)
		ext = ".fa" if fasta else ".fq"
		if gz:
			ext += ".gz"
		n = 2 if self.paired else 1
		paths = [directory / ("%s_%d%s" % (prefix, i + 1, ext)) for i in range(n)]
		for i, path in enumerate(paths):
			mate = (i + 1) if self.paired else None
			reads = self.mates(i + 1)
			lines = []
			for r in reads:
				name = _fastx_name(r.name, mate, name_style)
				if fasta:
					lines.append(">%s\n%s\n" % (name, r.seq))
				else:
					lines.append("@%s\n%s\n+\n%s\n" % (name, r.seq, r.qual))
			data = "".join(lines).encode()
			if gz:
				with gzip.open(path, "wb") as f:
					f.write(data)
			else:
				path.write_bytes(data)
		return paths


class ReadSimulator(object):
	"""Draws reads from a SyntheticGenome and records their truth.

	All randomness comes from `seed`, so the same calls give the same reads
	(also across Python versions). Base qualities are random Phred 20-40;
	const_qual='I' gives every base that quality character instead.
	"""

	def __init__(self, synth, seed=1, read_len=DEFAULT_READ_LEN, const_qual=None):
		self.g = synth
		self.rng = random.Random(seed)
		self.read_len = read_len
		self.const_qual = const_qual
		self._n = 0
		self._offsets = {}

	# -- basics
	def name(self, prefix="read"):
		self._n += 1
		return "%s%05d" % (prefix, self._n)

	def qual(self, n):
		if self.const_qual is not None:
			return self.const_qual * n
		return "".join(chr(33 + self.rng.randint(20, 40)) for _ in range(n))

	def random_seq(self, n):
		return "".join(self.rng.choice("ACGT") for _ in range(n))

	def _sub(self, b):
		return self.rng.choice([x for x in "ACGT" if x != b])

	def _build(self, contig, path, start, length, strand, mismatches=(), insertion=None,
			deletion=None, head="", tail="", n_offsets=()):
		"""Lay a read of `length` aligned read bases along `path` (ascending
		genomic positions) from index `start`. Offsets (mismatches, insertion,
		deletion, n_offsets) count read bases from the leftmost aligned base in
		genome orientation. Gaps in `path` become N, `deletion` becomes D.
		head/tail are extra bases at the read's 5'/3' end as sequenced."""
		ins_at, ins_n = insertion if insertion else (None, 0)
		del_at, del_n = deletion if deletion else (None, 0)
		n_ref = length - ins_n + del_n
		refs = path[start:start + n_ref]
		if start < 0 or len(refs) != n_ref:
			raise ValueError("read runs off the end of its template")
		cols = []  # (op, refpos)
		ri = 0
		off = 0
		ins_done = del_done = False
		while off < length:
			if not ins_done and off == ins_at:
				cols.extend(("I", None) for _ in range(ins_n))
				off += ins_n
				ins_done = True
				continue
			if not del_done and off == del_at:
				for _ in range(del_n):
					cols.append(("D", refs[ri]))
					ri += 1
				del_done = True
				continue
			cols.append(("M", refs[ri]))
			ri += 1
			off += 1
		ops = []
		fwd = []
		mm_genomic = []
		prev = None
		read_off = 0
		mismatches = set(mismatches)
		n_offsets = set(n_offsets)
		ctg = self.g.contigs[contig]
		for ci, (op, rp) in enumerate(cols):
			if op == "M":
				if prev is not None and rp != prev + 1:
					ops.append((rp - prev - 1, "N"))
				b = ctg[rp - 1]
				if read_off in mismatches:
					b = self._sub(b)
					mm_genomic.append(rp)
				if read_off in n_offsets:
					b = "N"
				fwd.append(b)
				ops.append((1, "M"))
				prev = rp
				read_off += 1
			elif op == "I":
				nxt = next(r for o, r in cols[ci:] if o == "M")
				flank = {ctg[nxt - 1]} if prev is None else {ctg[prev - 1], ctg[nxt - 1]}
				fwd.append(next(b for b in "ACGT" if b not in flank))
				ops.append((1, "I"))
				read_off += 1
			else:
				ops.append((1, "D"))
				prev = rp
		fwd = "".join(fwd)
		pos = refs[0]
		if strand == "+":
			left, right = head, tail
			seq = head + fwd + tail
		else:
			left, right = tail, head
			seq = head + revcomp(fwd) + tail
		cig = ops_to_cigar([(len(left), "S")] + ops + [(len(right), "S")])
		return seq, pos, cig, sorted(mm_genomic)

	def _make(self, name, seq, contig, strand, pos, cigar, mm, kind, tx=None, gene=None, extra=None, qual=None):
		return SimRead(name, seq, qual if qual is not None else self.qual(len(seq)), contig, strand, pos, cigar,
			mm, kind, tx, gene, extra=extra)

	# -- single reads
	def genomic_read(self, contig, pos, strand="+", length=None, mismatches=(), insertion=None,
			deletion=None, head="", tail="", n_offsets=(), name=None, kind="genomic"):
		"""Unspliced read whose leftmost aligned base is contig:pos.
		insertion/deletion are (read offset, n bases)."""
		length = length or self.read_len
		span = length + (deletion[1] if deletion else 0) + 10
		path = list(range(pos, min(pos + span, len(self.g.contigs[contig]) + 1)))
		seq, p, cig, mm = self._build(contig, path, 0, length, strand, mismatches, insertion, deletion,
			head, tail, n_offsets)
		return self._make(name or self.name(), seq, contig, strand, p, cig, mm, kind)

	def transcript_read(self, tx_id, offset, length=None, sense=True, mismatches=(), head="", tail="",
			n_offsets=(), name=None, kind=None):
		"""Read covering transcript bases [offset, offset+length) counted 5' to
		3' along the transcript; sense=False gives the antisense read. Spliced
		reads get N in their CIGAR. extra['tx_offset'] is `offset`."""
		t = self.g.transcripts[tx_id]
		length = length or self.read_len
		if offset < 0 or offset + length > t.length:
			raise ValueError("read [%d, %d) outside %s (length %d)" % (offset, offset + length, tx_id, t.length))
		path = t.genomic_positions()
		start = offset if t.strand == "+" else t.length - offset - length
		strand = t.strand if sense else ("-" if t.strand == "+" else "+")
		seq, p, cig, mm = self._build(t.contig, path, start, length, strand, mismatches, None, None,
			head, tail, n_offsets)
		if kind is None:
			kind = "spliced" if "N" in cig else "exonic"
		return self._make(name or self.name(), seq, t.contig, strand, p, cig, mm, kind, tx_id, t.gene_id,
			extra={"tx_offset": offset, "sense": sense})

	def segment_read(self, segment, offset, length=None, strand="+", copy=0, name=None):
		"""Read from repeated segment `segment` at `offset` in copy `copy`.
		Truth is that copy; extra['loci'] lists (contig, pos, strand) of every
		copy that the read matches, and extra['n_loci'] their number."""
		seg = self.g.segments[segment]
		length = length or self.read_len
		c = seg.copies[copy]
		if c.strand == "+":
			pos = c.start + offset
		else:
			pos = c.end - offset - length + 1
		read_strand = strand if c.strand == "+" else ("-" if strand == "+" else "+")
		r = self.genomic_read(c.contig, pos, read_strand, length, name=name, kind="multi")
		loci = []
		for k in seg.copies:
			if k.strand == "+":
				loci.append((k.contig, k.start + offset, strand))
			else:
				loci.append((k.contig, k.end - offset - length + 1, "-" if strand == "+" else "+"))
		r.extra["loci"] = loci
		r.extra["n_loci"] = len(loci)
		r.extra["segment"] = segment
		return r

	def unique_locus(self, length=None, avoid_genes=False, margin=50):
		"""(contig, pos) of a random window of `length` bases that overlaps no
		repeated segment (so a read there has one locus) and, with
		avoid_genes, no gene either."""
		length = length or self.read_len
		names = list(self.g.contigs)
		weights = [len(self.g.contigs[n]) for n in names]
		blocked = {}
		for seg in self.g.segments.values():
			for c in seg.copies:
				blocked.setdefault(c.contig, []).append((c.start, c.end))
		if avoid_genes:
			for gene in self.g.genes:
				blocked.setdefault(gene.contig, []).append((gene.start, gene.end))
		while True:
			contig = self.rng.choices(names, weights)[0]
			pos = self.rng.randint(margin, len(self.g.contigs[contig]) - length - margin)
			end = pos + length - 1
			if all(end + margin < a or pos - margin > b for a, b in blocked.get(contig, ())):
				return contig, pos

	def deletion_is_unambiguous(self, contig, start, n):
		"""True when deleting reference bases start..start+n-1 gives a read
		that no shifted deletion reproduces."""
		g = self.g.contigs[contig]
		return g[start - 2] != g[start + n - 2] and g[start - 1] != g[start + n - 1]

	def random_read(self, length=None, name=None, kind="random"):
		"""Random sequence that should not align."""
		length = length or self.read_len
		return self._make(name or self.name(), self.random_seq(length), None, None, None, None, (), kind)

	def chimeric_read(self, left, right, name=None):
		"""Read made of two genomic pieces, each (contig, pos, strand, length),
		joined in sequencing order: piece A is the donor, piece B the acceptor.
		The truth fields describe piece A. extra holds 'pieces', 'repeat'
		(left, right): how far the breakpoint could slide, 'motif': the
		Chimeric.out.junction junction type (1 GT/AG, 2 CT/AC, 0 other, read
		orientation), 'shift_motifs': {shift: junction type} for breakpoints
		moved by -8..8 read bases, and 'junction': the expected first 14 columns of
		Chimeric.out.junction (manual 6.4) for an unslidable breakpoint."""
		(ca, pa, sa, na), (cb, pb, sb, nb) = left, right
		ga, gb = self.g.contigs[ca], self.g.contigs[cb]
		k = 20
		if sa == "+":
			part_a = ga[pa - 1:pa - 1 + na]
			cont = ga[pa - 1 + na:pa - 1 + na + k]
		else:
			part_a = revcomp(ga[pa - 1:pa - 1 + na])
			cont = revcomp(ga[max(0, pa - 1 - k):pa - 1])
		if sb == "+":
			part_b = gb[pb - 1:pb - 1 + nb]
			prev = gb[max(0, pb - 1 - k):pb - 1]
		else:
			part_b = revcomp(gb[pb - 1:pb - 1 + nb])
			prev = revcomp(gb[pb - 1 + nb:pb - 1 + nb + k])
		rep_left = 0
		while rep_left < min(len(part_a), len(prev)) and part_a[-1 - rep_left] == prev[-1 - rep_left]:
			rep_left += 1
		rep_right = 0
		while rep_right < min(len(part_b), len(cont)) and part_b[rep_right] == cont[rep_right]:
			rep_right += 1
		motif = {"GTAG": 1, "CTAC": 2}.get(cont[:2] + prev[-2:], 0)
		# junction type if the breakpoint moved by s read bases (STAR may move it,
		# accepting mismatches, to reach GT/AG or CT/AC)
		a_ext, b_ext = part_a + cont, prev + part_b
		shift_motifs = {}
		for sh_ in range(-8, 9):
			donor = a_ext[na + sh_:na + sh_ + 2]
			acc = b_ext[len(prev) + sh_ - 2:len(prev) + sh_]
			shift_motifs[sh_] = {"GTAG": 1, "CTAC": 2}.get(donor + acc, 0)
		name = name or self.name()
		L = na + nb
		junction = OrderedDict([
			("chr_donorA", ca), ("brkpt_donorA", pa + na if sa == "+" else pa - 1), ("strand_donorA", sa),
			("chr_acceptorB", cb), ("brkpt_acceptorB", pb - 1 if sb == "+" else pb + nb), ("strand_acceptorB", sb),
			("junction_type", motif), ("repeat_left_lenA", rep_left), ("repeat_right_lenB", rep_right),
			("read_name", name),
			("start_alnA", pa), ("cigar_alnA", "%dM%dS" % (na, nb) if sa == "+" else "%dS%dM" % (nb, na)),
			("start_alnB", pb), ("cigar_alnB", "%dS%dM" % (na, nb) if sb == "+" else "%dM%dS" % (nb, na)),
		])
		seq = part_a + part_b
		return self._make(name, seq, ca, sa, pa, None, (), "chimeric",
			extra={"pieces": [left, right], "repeat": (rep_left, rep_right), "motif": motif,
				"shift_motifs": shift_motifs, "junction": junction, "length": L})

	def clean_chimeric_read(self, strands=("+", "+"), length=(50, 50), name=None):
		"""chimeric_read between unique, gene-free loci on different contigs
		whose breakpoint cannot slide (repeats 0) and has no GT/AG or CT/AC
		within 8 bases, so Chimeric.out.junction must match extra['junction']."""
		while True:
			ca, pa = self.unique_locus(length[0] + 10, avoid_genes=True)
			cb, pb = self.unique_locus(length[1] + 10, avoid_genes=True)
			if ca == cb:
				continue
			r = self.chimeric_read((ca, pa, strands[0], length[0]), (cb, pb, strands[1], length[1]), name=name)
			if r.extra["repeat"] == (0, 0) and not any(r.extra["shift_motifs"].values()):
				return r

	# -- batches
	def valid_offsets(self, tx_id, length=None, min_block=25, region=None):
		"""Transcript offsets whose read has every aligned block >= min_block.
		region=(lo, hi) restricts offsets to lo <= offset <= hi."""
		length = length or self.read_len
		key = (tx_id, length, min_block, region)
		if key in self._offsets:
			return self._offsets[key]
		t = self.g.transcripts[tx_id]
		bounds = []  # exon boundaries in transcript-sense offsets
		acc = 0
		lens = [e - s + 1 for s, e in t.exons]
		if t.strand == "-":
			lens = lens[::-1]
		for n in lens[:-1]:
			acc += n
			bounds.append(acc)
		lo, hi = 0, t.length - length
		if region:
			lo, hi = max(lo, region[0]), min(hi, region[1])
		out = []
		for o in range(lo, hi + 1):
			ok = True
			for b in bounds:
				if o < b < o + length and (b - o < min_block or o + length - b < min_block):
					ok = False
					break
			if ok:
				out.append(o)
		self._offsets[key] = out
		return out

	def transcript_reads(self, n, tx_ids=None, length=None, min_block=25, sense=None, spliced=None):
		"""n reads from annotated transcripts (all, or tx_ids). sense=None
		draws the orientation at random; spliced=True/False keeps only
		spliced/unspliced reads."""
		tx_ids = list(tx_ids or [t for t, x in self.g.transcripts.items() if x.annotated])
		length = length or self.read_len
		reads = []
		while len(reads) < n:
			tid = self.rng.choice(tx_ids)
			offs = self.valid_offsets(tid, length, min_block)
			o = self.rng.choice(offs)
			s = self.rng.random() < 0.5 if sense is None else sense
			r = self.transcript_read(tid, o, length, sense=s)
			if spliced is not None and (("N" in r.cigar) != spliced):
				continue
			reads.append(r)
		return ReadSet(reads)

	def transcript_pair(self, tx_id, offset, insert, length=None, forward=True, mismatches1=(),
			mismatches2=(), adapter=ILLUMINA_ADAPTER, name=None):
		"""Pair from fragment [offset, offset+insert) of a transcript (5' to 3').
		forward=True puts mate 1 on the transcript's sense strand. When insert
		< length the mates read into `adapter`."""
		length = length or self.read_len
		name = name or self.name("pair")
		L = min(length, insert)
		tail = (adapter * 3)[:length - L] if length > L else ""
		m5, m3 = (mismatches1, mismatches2) if forward else (mismatches2, mismatches1)
		five = self.transcript_read(tx_id, offset, L, sense=True, mismatches=m5, tail=tail, name=name)
		three = self.transcript_read(tx_id, offset + insert - L, L, sense=False, mismatches=m3, tail=tail, name=name)
		r1, r2 = (five, three) if forward else (three, five)
		r1.mate, r2.mate = 1, 2
		return SimPair(name, r1, r2, insert, "pair", extra={"tx": tx_id, "offset": offset})

	def genomic_pair(self, contig, pos, insert, strand="+", length=None, mismatches1=(), mismatches2=(),
			adapter=ILLUMINA_ADAPTER, name=None):
		"""Pair from fragment contig:pos..pos+insert-1. strand is mate 1's."""
		length = length or self.read_len
		name = name or self.name("pair")
		L = min(length, insert)
		tail = (adapter * 3)[:length - L] if length > L else ""
		ml, mr = (mismatches1, mismatches2) if strand == "+" else (mismatches2, mismatches1)
		left = self.genomic_read(contig, pos, "+", L, mismatches=ml, tail=tail, name=name)
		right = self.genomic_read(contig, pos + insert - L, "-", L, mismatches=mr, tail=tail, name=name)
		r1, r2 = (left, right) if strand == "+" else (right, left)
		r1.mate, r2.mate = 1, 2
		return SimPair(name, r1, r2, insert, "pair")

	def transcript_pairs(self, n, tx_ids=None, length=None, insert=(200, 300), min_block=25):
		"""n pairs from annotated transcripts with uniform insert sizes in
		[insert[0], insert[1]] and random orientation."""
		tx_ids = list(tx_ids or [t for t, x in self.g.transcripts.items()
			if x.annotated and x.length >= insert[1] + 10])
		length = length or self.read_len
		pairs = []
		while len(pairs) < n:
			tid = self.rng.choice(tx_ids)
			t = self.g.transcripts[tid]
			ins = self.rng.randint(insert[0], min(insert[1], t.length))
			valid1 = set(self.valid_offsets(tid, length, min_block))
			cands = [o for o in range(0, t.length - ins + 1)
				if o in valid1 and (o + ins - length) in valid1]
			if not cands:
				continue
			o = self.rng.choice(cands)
			pairs.append(self.transcript_pair(tid, o, ins, length, forward=self.rng.random() < 0.5))
		return ReadSet(pairs)

	def tenx_reads(self, whitelist, cells, cdna_len=90, cb_len=16, umi_len=12, reads_per_umi=1,
			three_prime=300, min_block=25, cb_errors=None):
		"""10x Chromium 3' (v3 geometry) reads.

		cells: {cell barcode: {gene_id: number of UMIs}}. Each UMI gives
		`reads_per_umi` cDNA reads (an int, or (lo, hi) drawn uniformly) from
		the last `three_prime` bases of the gene's first annotated transcript,
		in the sense orientation. UMIs within a cell are >= 3 mismatches
		apart. cb_errors: {cell barcode: number of reads whose barcode gets one
		substitution}, added on top. Returns ReadSet(kind='solo') whose r1 is
		the cDNA read and r2 is the CB+UMI read; truth for each pair is in
		pair.extra (cb, umi, gene, cb_read).
		"""
		pairs = []
		for cb in cells:
			if cb not in whitelist:
				raise ValueError("cell barcode %s not in whitelist" % cb)
		for cb, genes in cells.items():
			umis = []
			for gene_id, n_umi in genes.items():
				g = self.g.gene(gene_id)
				t = next(x for x in g.transcripts if x.annotated)
				offs = self.valid_offsets(t.tx_id, cdna_len, min_block,
					region=(t.length - three_prime, t.length - cdna_len))
				for _ in range(n_umi):
					while True:
						u = self.random_seq(umi_len)
						if all(hamming(u, v) >= 3 for v in umis) and len(set(u)) > 1:
							break
					umis.append(u)
					k = reads_per_umi if isinstance(reads_per_umi, int) else self.rng.randint(*reads_per_umi)
					for _ in range(k):
						pairs.append(self._tenx_pair(t, offs, cb, cb, u, cdna_len))
			for _ in range((cb_errors or {}).get(cb, 0)):
				gene_id = self.rng.choice(list(genes))
				t = next(x for x in self.g.gene(gene_id).transcripts if x.annotated)
				offs = self.valid_offsets(t.tx_id, cdna_len, min_block,
					region=(t.length - three_prime, t.length - cdna_len))
				i = self.rng.randrange(cb_len)
				bad = cb[:i] + self._sub(cb[i]) + cb[i + 1:]
				while True:
					u = self.random_seq(umi_len)
					if all(hamming(u, v) >= 3 for v in umis) and len(set(u)) > 1:
						break
				umis.append(u)
				pairs.append(self._tenx_pair(t, offs, cb, bad, u, cdna_len))
		return ReadSet(pairs, kind="solo")

	def _tenx_pair(self, t, offs, cb, cb_read, umi, cdna_len):
		name = self.name("cell")
		r1 = self.transcript_read(t.tx_id, self.rng.choice(offs), cdna_len, sense=True, name=name)
		bc = cb_read + umi
		r2 = SimRead(name, bc, self.qual(len(bc)), kind="barcode")
		r1.mate, r2.mate = 1, 2
		return SimPair(name, r1, r2, None, "solo",
			extra={"cb": cb, "umi": umi, "gene": t.gene_id, "cb_read": cb_read})


def make_whitelist(n, seed=7, length=16, min_dist=3):
	"""n random cell barcodes, pairwise >= min_dist mismatches apart."""
	rng = random.Random(seed)
	out = []
	while len(out) < n:
		b = "".join(rng.choice("ACGT") for _ in range(length))
		if all(hamming(b, x) >= min_dist for x in out):
			out.append(b)
	return out


def expected_solo_counts(readset):
	"""{(gene_id, cell barcode): number of distinct UMIs} from 10x truth."""
	umis = {}
	for p in readset:
		key = (p.extra["gene"], p.extra["cb"])
		umis.setdefault(key, set()).add(p.extra["umi"])
	return {k: len(v) for k, v in umis.items()}


# ------------------------------------------------------------- running STAR

def find_star():
	"""Path to the STAR binary from $STAR_BIN, else source/STAR; None if
	neither is an executable file."""
	cand = os.environ.get("STAR_BIN") or str(REPO_ROOT / "source" / "STAR")
	p = Path(cand).expanduser()
	if p.is_file() and os.access(str(p), os.X_OK):
		return p.resolve()
	return None


def default_threads():
	"""--runThreadN used when a test does not pass threads: $STAR_TEST_THREADS
	or 1."""
	n = int(os.environ.get("STAR_TEST_THREADS", "1"))
	if not 1 <= n <= MAX_THREADS:
		raise ValueError("STAR_TEST_THREADS must be between 1 and %d" % MAX_THREADS)
	return n


_RUN_ROOT = [None]


def set_run_root(path):
	"""Directory under which run_star creates one fresh directory per run
	(conftest points it at a pytest tmp directory)."""
	_RUN_ROOT[0] = Path(path).resolve()
	_RUN_ROOT[0].mkdir(parents=True, exist_ok=True)


def _fresh_dir(name):
	root = _RUN_ROOT[0]
	if root is None:
		root = Path(tempfile.mkdtemp(prefix="star_tests_"))
		_RUN_ROOT[0] = root
	safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name or "run")[:60]
	return Path(tempfile.mkdtemp(prefix=safe + "_", dir=str(root)))


def _split_args(extra_args):
	if extra_args is None:
		return []
	if isinstance(extra_args, str):
		return shlex.split(extra_args)
	return [str(a) for a in extra_args]


def _check_threads(threads, args):
	if "--runThreadN" in args:
		raise ValueError("pass threads=, not --runThreadN, so the thread count stays pinned")
	threads = default_threads() if threads is None else int(threads)
	if not 1 <= threads <= MAX_THREADS:
		raise ValueError("threads must be between 1 and %d" % MAX_THREADS)
	return threads


def _run(argv, cwd, timeout):
	return subprocess.run(argv, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
		timeout=timeout or DEFAULT_TIMEOUT)


class StarError(AssertionError):
	"""STAR exited non-zero (raised by run_star/build_genome with check=True)."""


def _failure(argv, proc, outdir):
	log = ""
	for name in ("Log.out",):
		p = Path(outdir) / name
		if p.exists():
			log = "\n".join(p.read_text(errors="replace").splitlines()[-15:])
	return StarError("STAR exited %d\nargv: %s\nstderr: %s\nstdout: %s\nLog.out tail:\n%s" % (
		proc.returncode, " ".join(shlex.quote(a) for a in argv),
		proc.stderr.decode(errors="replace")[-2000:], proc.stdout.decode(errors="replace")[-1000:], log))


class StarGenome(object):
	"""A generated STAR index. dir: index directory; synth: the
	SyntheticGenome; gtf: path of the GTF used (None without sjdb);
	sjdb_overhang: --sjdbOverhang (None without sjdb)."""

	def __init__(self, directory, synth, gtf=None, sjdb_overhang=None, args=()):
		self.dir = Path(directory)
		self.synth = synth
		self.gtf = gtf
		self.sjdb_overhang = sjdb_overhang
		self.args = list(args)

	def __repr__(self):
		return "StarGenome(%s, gtf=%s)" % (self.dir, bool(self.gtf))


def build_genome(synth, outdir, gtf=True, sjdb_overhang=DEFAULT_SJDB_OVERHANG, extra_args=(),
		threads=None, sa_index_nbases=None, timeout=None):
	"""Run genomeGenerate for `synth` (which must have been write()-n) into
	`outdir`. gtf=True uses synth.gtf_path; a path uses that GTF; False builds
	without annotations. Returns a StarGenome."""
	star = find_star()
	if star is None:
		raise StarError(BUILD_HINT)
	args = _split_args(extra_args)
	threads = _check_threads(threads, args)
	outdir = Path(outdir).resolve()
	outdir.mkdir(parents=True, exist_ok=True)
	n = sa_index_nbases if sa_index_nbases is not None else synth.sa_index_nbases()
	argv = [str(star), "--runMode", "genomeGenerate", "--runThreadN", str(threads),
		"--genomeDir", str(outdir), "--genomeFastaFiles", str(synth.fasta_path),
		"--genomeSAindexNbases", str(n), "--outFileNamePrefix", str(outdir) + "/"]
	gtf_path = None
	if gtf:
		gtf_path = Path(synth.gtf_path if gtf is True else gtf).resolve()
		argv += ["--sjdbGTFfile", str(gtf_path), "--sjdbOverhang", str(sjdb_overhang)]
	argv += args
	proc = _run(argv, outdir, timeout)
	if proc.returncode != 0:
		raise _failure(argv, proc, outdir)
	return StarGenome(outdir, synth, gtf_path, sjdb_overhang if gtf else None, args)


def run_star(genome, reads, extra_args=(), threads=None, outdir=None, name=None, fasta=False, gz=False,
		name_style="plain", check=True, timeout=None):
	"""Align `reads` with STAR in a fresh directory and return a StarRun.

	genome: StarGenome (or index directory path).
	reads: ReadSet; a list of ReadSets (joined with commas, one read group per
	set); a list of file paths given verbatim to --readFilesIn (one entry per
	mate, commas allowed); or None when extra_args carries --readFilesIn.
	extra_args: list of strings, or one string split with shlex.
	threads: --runThreadN; None means default_threads() ($STAR_TEST_THREADS,
	else 1). Never pass --runThreadN in extra_args.
	fasta/gz/name_style: how ReadSets are written; gz adds
	--readFilesCommand gzip -cd unless extra_args sets --readFilesCommand.
	check: raise StarError (an AssertionError) if STAR exits non-zero.
	"""
	star = find_star()
	if star is None:
		raise StarError(BUILD_HINT)
	args = _split_args(extra_args)
	threads = _check_threads(threads, args)
	out = Path(outdir).resolve() if outdir is not None else _fresh_dir(name)
	out.mkdir(parents=True, exist_ok=True)
	gdir = (genome.dir if isinstance(genome, StarGenome) else Path(genome)).resolve()
	argv = [str(star), "--runMode", "alignReads", "--runThreadN", str(threads), "--genomeDir", str(gdir),
		"--outFileNamePrefix", str(out) + "/"]
	sets = []
	if isinstance(reads, ReadSet):
		sets = [reads]
	elif isinstance(reads, (list, tuple)) and reads and all(isinstance(r, ReadSet) for r in reads):
		sets = list(reads)
	if sets:
		files = [s.write(out / "input", prefix="reads%d" % i, fasta=fasta, gz=gz, name_style=name_style)
			for i, s in enumerate(sets)]
		nf = {len(f) for f in files}
		if len(nf) != 1:
			raise ValueError("cannot mix single-end and paired read sets")
		argv += ["--readFilesIn"] + [",".join(str(f[i]) for f in files) for i in range(nf.pop())]
		if gz and "--readFilesCommand" not in args:
			argv += ["--readFilesCommand", "gzip", "-cd"]
	elif reads is not None:
		argv += ["--readFilesIn"] + [",".join(str(Path(x).resolve()) for x in str(r).split(",")) for r in reads]
	argv += args
	proc = _run(argv, out, timeout)
	run = StarRun(out, argv, proc, sets, genome)
	if check and proc.returncode != 0:
		raise _failure(argv, proc, out)
	return run


class StarRun(object):
	"""Result of run_star. outdir: run directory; argv; returncode;
	stdout_bytes/stdout/stderr; reads: the ReadSets aligned (truth);
	truth: dict read name -> SimRead or SimPair over all of them."""

	def __init__(self, outdir, argv, proc, readsets, genome):
		self.outdir = Path(outdir)
		self.argv = argv
		self.returncode = proc.returncode
		self.stdout_bytes = proc.stdout
		self.stderr = proc.stderr.decode(errors="replace")
		self.reads = readsets
		self.genome = genome
		self.truth = OrderedDict()
		for s in readsets:
			self.truth.update(s.by_name)

	@property
	def stdout(self):
		return self.stdout_bytes.decode(errors="replace")

	def path(self, name):
		return self.outdir / name

	def exists(self, name):
		return (self.outdir / name).exists()

	def main_alignment_path(self):
		for n in ("Aligned.out.sam", "Aligned.out.bam", "Aligned.sortedByCoord.out.bam"):
			if self.exists(n):
				return self.path(n)
		raise FileNotFoundError("no Aligned.* output in %s" % self.outdir)

	def alignments(self, name=None):
		"""AlignmentFile of the main output (Aligned.out.sam, Aligned.out.bam
		or Aligned.sortedByCoord.out.bam, first found) or of `name`."""
		p = self.path(name) if name else self.main_alignment_path()
		return read_alignments(p)

	def records(self, name=None):
		"""List of SamRecord from alignments(name)."""
		return self.alignments(name).records

	def transcriptome(self):
		"""AlignmentFile of Aligned.toTranscriptome.out.bam."""
		return read_alignments(self.path("Aligned.toTranscriptome.out.bam"))

	def log_final(self):
		"""Log.final.out as a dict (parse_log_final)."""
		return parse_log_final(self.path("Log.final.out"))

	def log_out(self):
		return self.path("Log.out").read_text(errors="replace")

	def sj(self, name="SJ.out.tab"):
		"""SJ.out.tab as a list of SJRecord, in file order."""
		return parse_sj_out(self.path(name))

	def gene_counts(self):
		"""ReadsPerGene.out.tab as {name: (unstranded, strand1, strand2)}."""
		return parse_reads_per_gene(self.path("ReadsPerGene.out.tab"))

	def unmapped(self, mate=1):
		"""Unmapped.out.mate<mate> as a list of FastxRecord."""
		return parse_fastx(self.path("Unmapped.out.mate%d" % mate))

	def chimeric(self):
		"""Chimeric.out.junction as a list of dicts (parse_chimeric_junction)."""
		return parse_chimeric_junction(self.path("Chimeric.out.junction"))

	def solo_out_names(self):
		"""The four --soloOutFileNames of this run (output prefix, features,
		barcodes, matrix), defaults from parametersDefault for those not
		given."""
		names = list(SOLO_OUT_FILE_NAMES)
		if "--soloOutFileNames" in self.argv:
			i = self.argv.index("--soloOutFileNames") + 1
			given = []
			while i < len(self.argv) and not self.argv[i].startswith("--"):
				given.append(self.argv[i])
				i += 1
			names[:len(given)] = given[:4]
		return names

	def solo_dir(self, feature="Gene"):
		"""Solo output directory of `feature` (Solo.out/<feature> by
		default; STAR prepends the first --soloOutFileNames value as is)."""
		return Path(str(self.outdir) + "/" + self.solo_out_names()[0] + feature)

	def solo(self, feature="Gene", matrix="raw", mtx=None):
		"""SoloMatrix from <solo_dir(feature)>/<matrix>/; the features,
		barcodes and (unless `mtx` is given) matrix file names follow
		--soloOutFileNames."""
		_, features, barcodes, default_mtx = self.solo_out_names()
		return read_solo_matrix(self.solo_dir(feature) / matrix, mtx or default_mtx, features, barcodes)

	def solo_features_stats(self, feature="Gene"):
		"""Solo.out/<feature>/Features.stats as {name: number}."""
		return parse_stats(self.solo_dir(feature) / "Features.stats")

	def solo_summary(self, feature="Gene"):
		"""Solo.out/<feature>/Summary.csv as {field: value}."""
		return parse_summary_csv(self.solo_dir(feature) / "Summary.csv")

	def solo_barcodes_stats(self):
		"""Solo.out/Barcodes.stats as {name: number}."""
		return parse_stats(self.solo_dir("Barcodes.stats"))


# ---------------------------------------------------------------- SAM / BAM

FLAG_PAIRED = 0x1
FLAG_PROPER = 0x2
FLAG_UNMAPPED = 0x4
FLAG_MATE_UNMAPPED = 0x8
FLAG_REVERSE = 0x10
FLAG_MATE_REVERSE = 0x20
FLAG_READ1 = 0x40
FLAG_READ2 = 0x80
FLAG_SECONDARY = 0x100
FLAG_QCFAIL = 0x200
FLAG_DUP = 0x400
FLAG_SUPPLEMENTARY = 0x800


class SamRecord(object):
	"""One alignment line. pos/pnext are 1-based (0 when unset), rname/rnext
	are names or '*'. tags: ordered dict tag -> value (int, float, str, or a
	list for B arrays); tag_types: tag -> 'i', 'f', 'A', 'Z', 'H' or 'B:<sub>'
	(BAM integer types are all reported as 'i', as in SAM text). For records
	read from BAM, bin is the bin field and bam_types maps each tag to its
	type code as stored ('c', 'C', 's', 'S', 'i', 'I', 'A', 'f', 'Z', 'H' or
	'B:<sub>'); both are None for SAM records."""

	__slots__ = ("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual",
		"tags", "tag_types", "bin", "bam_types")

	def __init__(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags, tag_types, bin=None,
			bam_types=None):
		self.qname = qname
		self.flag = flag
		self.rname = rname
		self.pos = pos
		self.mapq = mapq
		self.cigar = cigar
		self.rnext = rnext
		self.pnext = pnext
		self.tlen = tlen
		self.seq = seq
		self.qual = qual
		self.tags = tags
		self.tag_types = tag_types
		self.bin = bin
		self.bam_types = bam_types

	is_paired = property(lambda s: bool(s.flag & FLAG_PAIRED))
	is_proper_pair = property(lambda s: bool(s.flag & FLAG_PROPER))
	is_unmapped = property(lambda s: bool(s.flag & FLAG_UNMAPPED))
	mate_unmapped = property(lambda s: bool(s.flag & FLAG_MATE_UNMAPPED))
	is_reverse = property(lambda s: bool(s.flag & FLAG_REVERSE))
	mate_reverse = property(lambda s: bool(s.flag & FLAG_MATE_REVERSE))
	is_read1 = property(lambda s: bool(s.flag & FLAG_READ1))
	is_read2 = property(lambda s: bool(s.flag & FLAG_READ2))
	is_secondary = property(lambda s: bool(s.flag & FLAG_SECONDARY))
	is_supplementary = property(lambda s: bool(s.flag & FLAG_SUPPLEMENTARY))

	@property
	def strand(self):
		return "-" if self.is_reverse else "+"

	@property
	def mate(self):
		"""1 or 2 for paired records, None for single-end."""
		if self.is_read1:
			return 1
		if self.is_read2:
			return 2
		return None

	@property
	def end(self):
		"""Rightmost aligned reference base (1-based)."""
		return self.pos + cigar_ref_len(self.cigar) - 1

	@property
	def introns(self):
		return cigar_introns(self.pos, self.cigar)

	def tag(self, name, default=None):
		return self.tags.get(name, default)

	def to_line(self):
		fields = [self.qname, str(self.flag), self.rname, str(self.pos), str(self.mapq), self.cigar,
			self.rnext, str(self.pnext), str(self.tlen), self.seq, self.qual]
		for t, v in self.tags.items():
			ty = self.tag_types[t]
			if ty.startswith("B:"):
				fields.append("%s:B:%s,%s" % (t, ty[2:], ",".join(_fmt_num(x) for x in v)))
			else:
				fields.append("%s:%s:%s" % (t, ty, _fmt_num(v)))
		return "\t".join(fields)

	def __repr__(self):
		return "SamRecord(%s)" % self.to_line()


def _fmt_num(v):
	if isinstance(v, float):
		return "%g" % v
	return str(v)


def expected_nm_md(synth, contig, pos, cigar, sam_seq):
	"""(NM, MD) for an alignment of `sam_seq` (SAM orientation) at
	contig:pos with `cigar`, following the SAM specification; N in the read
	or genome counts as a mismatch, as in STAR."""
	g = synth.contigs[contig]
	md = []
	run = 0
	nm = 0
	rp = pos
	qp = 0
	for n, op in cigar_ops(cigar):
		if op in "M=X":
			for i in range(n):
				rb = g[rp - 1 + i]
				qb = sam_seq[qp + i]
				if rb != qb or rb == "N" or qb == "N":
					nm += 1
					md.append(str(run))
					md.append(rb)
					run = 0
				else:
					run += 1
			rp += n
			qp += n
		elif op == "I":
			nm += n
			qp += n
		elif op == "D":
			nm += n
			md.append(str(run))
			md.append("^" + g[rp - 1:rp - 1 + n])
			run = 0
			rp += n
		elif op == "N":
			rp += n
		elif op == "S":
			qp += n
	md.append(str(run))
	return nm, "".join(md)


def record_truth(rec):
	"""(rname, pos, cigar, strand) of a record, comparable with
	SimRead.truth()."""
	return (rec.rname, rec.pos, rec.cigar, rec.strand)


class AlignmentFile(object):
	"""header: list of header lines; references: list of (name, length);
	records: list of SamRecord in file order."""

	def __init__(self, header, references, records):
		self.header = header
		self.references = references
		self.records = records

	def __iter__(self):
		return iter(self.records)

	def __len__(self):
		return len(self.records)

	def header_lines(self, kind):
		"""Header lines of one type, e.g. '@RG'."""
		return [h for h in self.header if h.startswith(kind + "\t") or h == kind]

	def by_name(self):
		"""{qname: [records]} in file order."""
		out = OrderedDict()
		for r in self.records:
			out.setdefault(r.qname, []).append(r)
		return out

	def primary(self):
		"""Mapped records that are neither secondary nor supplementary."""
		return [r for r in self.records if not (r.flag & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY))]

	def mapped(self):
		return [r for r in self.records if not r.is_unmapped]


def _parse_tag(field):
	t, ty, val = field.split(":", 2)
	if ty == "i":
		return t, int(val), "i"
	if ty == "f":
		return t, float(val), "f"
	if ty == "B":
		sub, _, rest = val.partition(",")
		conv = float if sub == "f" else int
		return t, [conv(x) for x in rest.split(",")] if rest else [], "B:" + sub
	return t, val, ty


def parse_sam_line(line):
	"""One SAM text line -> SamRecord."""
	f = line.rstrip("\n").split("\t")
	if len(f) < 11:
		raise ValueError("SAM line with %d fields: %r" % (len(f), line))
	tags = OrderedDict()
	types = {}
	for x in f[11:]:
		t, v, ty = _parse_tag(x)
		tags[t] = v
		types[t] = ty
	return SamRecord(f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[5], f[6], int(f[7]), int(f[8]),
		f[9], f[10], tags, types)


def parse_sam(source):
	"""AlignmentFile from a SAM path, or from SAM text if `source` is a str
	containing newlines."""
	if isinstance(source, str) and "\n" in source:
		text = source
	else:
		text = Path(source).read_text()
	header, records, refs = [], [], []
	for line in text.splitlines():
		if not line:
			continue
		if line.startswith("@"):
			header.append(line)
			if line.startswith("@SQ\t"):
				d = dict(x.split(":", 1) for x in line.split("\t")[1:])
				refs.append((d["SN"], int(d["LN"])))
		else:
			records.append(parse_sam_line(line))
	return AlignmentFile(header, refs, records)


_SEQ_NIB = "=ACMGRSVTWYHKDBN"
_SEQ_BYTE = [_SEQ_NIB[b >> 4] + _SEQ_NIB[b & 15] for b in range(256)]
_BAM_INT = {"c": "<b", "C": "<B", "s": "<h", "S": "<H", "i": "<i", "I": "<I"}


def read_bam(source):
	"""AlignmentFile from a BAM path (or its bytes). BGZF is multi-member
	gzip, so gzip + struct decode it without pysam. Each record keeps the
	BAM bin field in `bin` and the stored tag type codes in `bam_types`."""
	raw = source if isinstance(source, (bytes, bytearray)) else Path(source).read_bytes()
	data = gzip.decompress(raw)
	if data[:4] != b"BAM\x01":
		raise ValueError("not a BAM file")
	(l_text,) = struct.unpack_from("<i", data, 4)
	text = data[8:8 + l_text].decode().rstrip("\x00")
	off = 8 + l_text
	(n_ref,) = struct.unpack_from("<i", data, off)
	off += 4
	refs = []
	for _ in range(n_ref):
		(l_name,) = struct.unpack_from("<i", data, off)
		off += 4
		name = data[off:off + l_name - 1].decode()
		off += l_name
		(l_ref,) = struct.unpack_from("<i", data, off)
		off += 4
		refs.append((name, l_ref))
	names = [r[0] for r in refs]
	records = []
	n = len(data)
	while off < n:
		(block,) = struct.unpack_from("<i", data, off)
		start = off + 4
		end = start + block
		(ref_id, pos, l_qname, mapq, bin_, n_cig, flag, l_seq, nref_id, npos, tlen) = struct.unpack_from(
			"<iiBBHHHiiii", data, start)
		p = start + 32
		qname = data[p:p + l_qname - 1].decode()
		p += l_qname
		cig = struct.unpack_from("<%dI" % n_cig, data, p)
		p += 4 * n_cig
		cigar = "".join("%d%s" % (c >> 4, "MIDNSHP=X"[c & 15]) for c in cig) or "*"
		nb = (l_seq + 1) // 2
		seq = "".join(_SEQ_BYTE[b] for b in data[p:p + nb])[:l_seq] or "*"
		p += nb
		qb = data[p:p + l_seq]
		p += l_seq
		qual = "*" if (l_seq == 0 or qb[0] == 0xFF) else "".join(chr(q + 33) for q in qb)
		tags = OrderedDict()
		types = {}
		raw_types = OrderedDict()
		while p < end:
			t = data[p:p + 2].decode()
			ty = chr(data[p + 2])
			raw_types[t] = ty
			p += 3
			if ty == "A":
				v = chr(data[p])
				p += 1
			elif ty in _BAM_INT:
				fmt = _BAM_INT[ty]
				(v,) = struct.unpack_from(fmt, data, p)
				p += struct.calcsize(fmt)
				ty = "i"
			elif ty == "f":
				(v,) = struct.unpack_from("<f", data, p)
				p += 4
			elif ty in "ZH":
				z = data.index(b"\x00", p)
				v = data[p:z].decode()
				p = z + 1
			elif ty == "B":
				sub = chr(data[p])
				(cnt,) = struct.unpack_from("<i", data, p + 1)
				p += 5
				fmt = "<%d%s" % (cnt, "f" if sub == "f" else _BAM_INT[sub][1])
				v = list(struct.unpack_from(fmt, data, p))
				p += struct.calcsize(fmt)
				ty = "B:" + sub
				raw_types[t] = ty
			else:
				raise ValueError("unknown BAM tag type %r" % ty)
			tags[t] = v
			types[t] = ty
		rname = names[ref_id] if ref_id >= 0 else "*"
		if nref_id < 0:
			rnext = "*"
		elif nref_id == ref_id:
			rnext = "="
		else:
			rnext = names[nref_id]
		records.append(SamRecord(qname, flag, rname, pos + 1, mapq, cigar, rnext, npos + 1, tlen, seq, qual,
			tags, types, bin_, raw_types))
		off = end
	return AlignmentFile(text.splitlines(), refs, records)


def read_alignments(path):
	"""parse_sam or read_bam by file extension."""
	path = Path(path)
	return read_bam(path) if path.suffix == ".bam" else parse_sam(path)


# --------------------------------------------------------- other outputs

SJRecord = namedtuple("SJRecord", "contig start end strand motif annotated unique multi max_overhang")


def parse_sj_out(path):
	"""SJ.out.tab -> list of SJRecord (all fields but contig are int)."""
	out = []
	for line in Path(path).read_text().splitlines():
		if line.strip():
			f = line.split("\t")
			out.append(SJRecord(f[0], *[int(x) for x in f[1:9]]))
	return out


def parse_reads_per_gene(path):
	"""ReadsPerGene.out.tab -> {name: (unstranded, strand1, strand2)} for
	the N_unmapped/N_multimapping/N_noFeature/N_ambiguous rows and every
	gene ID."""
	out = OrderedDict()
	for line in Path(path).read_text().splitlines():
		if line.strip():
			f = line.split("\t")
			out[f[0]] = tuple(int(x) for x in f[1:4])
	return out


def _num(v):
	v = v.strip()
	if v.endswith("%"):
		try:
			return float(v[:-1])
		except ValueError:
			return v
	try:
		return int(v)
	except ValueError:
		pass
	try:
		return float(v)
	except ValueError:
		return v


def parse_log_final(path):
	"""Log.final.out -> {field name: value}; integers become int, 'x%'
	becomes float x, other values stay str. Names are as printed, stripped,
	e.g. 'Number of input reads', 'Uniquely mapped reads %'."""
	out = OrderedDict()
	for line in Path(path).read_text().splitlines():
		if "|" in line:
			k, v = line.split("|", 1)
			out[k.strip()] = _num(v)
	return out


FastxRecord = namedtuple("FastxRecord", "name description seq qual")


def parse_fastx(path):
	"""FASTQ or FASTA -> list of FastxRecord(name, description, seq, qual);
	name is the first word without '@'/'>', qual is None for FASTA."""
	lines = Path(path).read_text().splitlines()
	out = []
	i = 0
	while i < len(lines):
		h = lines[i]
		if not h:
			i += 1
			continue
		name, _, desc = h[1:].partition(" ")
		if h[0] == "@":
			out.append(FastxRecord(name, desc, lines[i + 1], lines[i + 3]))
			i += 4
		elif h[0] == ">":
			seq = []
			i += 1
			while i < len(lines) and not lines[i].startswith(">"):
				seq.append(lines[i])
				i += 1
			out.append(FastxRecord(name, desc, "".join(seq), None))
		else:
			raise ValueError("not FASTA/FASTQ at line %d of %s" % (i + 1, path))
	return out


CHIMERIC_FIELDS = ("chr_donorA", "brkpt_donorA", "strand_donorA", "chr_acceptorB", "brkpt_acceptorB",
	"strand_acceptorB", "junction_type", "repeat_left_lenA", "repeat_right_lenB", "read_name",
	"start_alnA", "cigar_alnA", "start_alnB", "cigar_alnB")
_CHIM_INT = {"brkpt_donorA", "brkpt_acceptorB", "junction_type", "repeat_left_lenA", "repeat_right_lenB",
	"start_alnA", "start_alnB"}


def parse_chimeric_junction(path):
	"""Chimeric.out.junction -> list of dicts keyed by CHIMERIC_FIELDS
	(integers converted); columns past the 14th go to 'extra' (list of str).
	Header and '#' comment lines are skipped."""
	out = []
	for line in Path(path).read_text().splitlines():
		if not line.strip() or line.startswith("#") or line.startswith("chr_donorA"):
			continue
		f = line.split("\t")
		d = OrderedDict()
		for k, v in zip(CHIMERIC_FIELDS, f):
			d[k] = int(v) if k in _CHIM_INT else v
		d["extra"] = f[len(CHIMERIC_FIELDS):]
		out.append(d)
	return out


class SoloMatrix(object):
	"""A Solo.out matrix directory. features: list of tuples from
	features.tsv (gene_id, gene_name, type; for the SJ feature the
	SJ.out.tab columns); barcodes: list of str; keys: the feature key of
	each matrix row (solo_feature_keys); counts: {(feature key, barcode):
	value} for non-zero entries (value is an int, or a tuple when the .mtx
	has several value columns)."""

	def __init__(self, features, barcodes, counts, shape, keys=None):
		self.features = features
		self.barcodes = barcodes
		self.counts = counts
		self.shape = shape
		self.keys = keys if keys is not None else [f[0] for f in features]

	def get(self, feature_id, barcode, default=0):
		return self.counts.get((feature_id, barcode), default)

	def feature_ids(self):
		return [f[0] for f in self.features]

	def cell_totals(self):
		"""{barcode: sum of counts} over non-zero cells."""
		out = {}
		for (_, b), v in self.counts.items():
			out[b] = out.get(b, 0) + (v if not isinstance(v, tuple) else v[0])
		return out


def parse_mtx(path):
	"""MatrixMarket coordinate file -> (nrows, ncols, {(row, col): value})
	with 1-based indices."""
	rows = cols = None
	entries = {}
	for line in Path(path).read_text().splitlines():
		if not line or line.startswith("%"):
			continue
		f = line.split()
		if rows is None:
			rows, cols, n_entries = int(f[0]), int(f[1]), int(f[2])
			continue
		vals = tuple(_num(x) for x in f[2:])
		entries[(int(f[0]), int(f[1]))] = vals[0] if len(vals) == 1 else vals
	if rows is not None and len(entries) != n_entries:
		raise ValueError("%s: header says %d entries, found %d" % (path, n_entries, len(entries)))
	return rows, cols, entries


def solo_feature_keys(features):
	"""Key of each features.tsv row: column 1 (the gene ID) when column 1 is
	unique, else (contig, start, end) with integer coordinates when columns
	2 and 3 are integers (the SJ feature, whose features.tsv is SJ.out.tab),
	else the whole row."""
	ids = [f[0] for f in features]
	if len(set(ids)) == len(ids):
		return ids
	if all(len(f) >= 3 and f[1].isdigit() and f[2].isdigit() for f in features):
		keys = [(f[0], int(f[1]), int(f[2])) for f in features]
	else:
		keys = [tuple(f) for f in features]
	if len(set(keys)) != len(keys):
		raise ValueError("features.tsv rows are not unique")
	return keys


def read_solo_matrix(directory, mtx="matrix.mtx", features="features.tsv", barcodes="barcodes.tsv"):
	"""SoloMatrix from a Solo.out/<feature>/<raw|filtered> directory; the
	file names default to STAR's (--soloOutFileNames)."""
	directory = Path(directory)
	feats = [tuple(l.split("\t")) for l in (directory / features).read_text().splitlines() if l]
	bcs = [l for l in (directory / barcodes).read_text().splitlines() if l]
	nr, nc, entries = parse_mtx(directory / mtx)
	keys = solo_feature_keys(feats)
	counts = {(keys[r - 1], bcs[c - 1]): v for (r, c), v in entries.items()}
	return SoloMatrix(feats, bcs, counts, (nr, nc), keys)


def parse_stats(path):
	"""Whitespace 'name value' lines (Features.stats, Barcodes.stats) ->
	{name: number}."""
	out = OrderedDict()
	for line in Path(path).read_text().splitlines():
		f = line.split()
		if len(f) >= 2:
			out[f[0]] = _num(f[-1])
	return out


def parse_summary_csv(path):
	"""Solo Summary.csv -> {field: number or str}."""
	out = OrderedDict()
	for line in Path(path).read_text().splitlines():
		if "," in line:
			k, v = line.rsplit(",", 1)
			out[k] = _num(v)
	return out


# --------------------------------------------------------- comparisons

def record_multiset(records, fields=None, bam_types=False):
	"""Sorted list of record lines (or of tuples of `fields`) for comparing
	outputs whose order is not fixed, e.g. --runThreadN 1 vs 4. bam_types:
	append each BAM record's stored tag type codes to its line, so two BAM
	outputs differing only in integer tag types differ (leave it False when
	comparing SAM with BAM)."""
	if fields is None:
		if bam_types:
			return sorted(r.to_line() + "\t" + ",".join("%s:%s" % kv for kv in (r.bam_types or {}).items())
				for r in records)
		return sorted(r.to_line() for r in records)
	return sorted(tuple(getattr(r, f) for f in fields) for r in records)


def first_line_diff(a, b, limit=5):
	"""Human-readable summary of how two sorted line lists differ."""
	sa, sb = set(a), set(b)
	only_a = [x for x in a if x not in sb][:limit]
	only_b = [x for x in b if x not in sa][:limit]
	return "only in first (%d): %s\nonly in second (%d): %s" % (
		len(set(a) - sb), only_a, len(sb - sa), only_b)
