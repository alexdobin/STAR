"""Helpers for the area A tests (input, seeds, stitching, scoring, filters,
multimappers). Each helper states the STAR rule it models so that expected
values in the tests can be checked against the source by hand."""

import math

import star_helpers as sh
from star_helpers import SimRead

_OTHER = {"A": "C", "C": "G", "G": "T", "T": "A"}

UNMAPPED_FIELDS = (
	"Number of reads unmapped: too many mismatches",
	"Number of reads unmapped: too short",
	"Number of reads unmapped: other",
	"Number of reads mapped to too many loci",
)


def unmapped_counts(log):
	"""(too many mismatches, too short, other, too many loci) from a
	Log.final.out dict."""
	return tuple(log[f] for f in UNMAPPED_FIELDS)


def length_score(span, scale=-0.25):
	"""Genomic-length term of the alignment score (stitchWindowAligns.cpp):
	int(ceil(log2(span) * scoreGenomicLengthLog2scale - 0.5)), where span is
	the genomic distance from the first to the last aligned base + 1. STAR
	adds it only when the scale is non-zero."""
	if scale == 0:
		return 0
	return int(math.ceil(math.log2(span) * scale - 0.5))


def mismatched(seq):
	"""A sequence of the same length that differs from `seq` at every base,
	so a local extension over it never gains score."""
	return "".join(_OTHER[b] for b in seq.upper())


def make_read(name, seq, contig=None, pos=None, cigar=None, strand="+", qual=None, extra=None):
	"""SimRead with a hand-given truth (seq is as sequenced)."""
	return SimRead(name, seq, qual or "I" * len(seq), contig, strand if contig else None, pos, cigar,
		extra=extra or {})


def oriented(name, fwd, contig, pos, cigar, strand, extra=None):
	"""SimRead from its genome-orientation sequence `fwd`: sequenced as fwd
	for '+' and as its reverse complement for '-'."""
	seq = fwd if strand == "+" else sh.revcomp(fwd)
	return make_read(name, seq, contig, pos, cigar, strand, extra=extra)


def local_extension(read, genome, extend_to_end=False):
	"""Number of bases a seed end is extended over `read` against `genome`
	(both listed outward from the seed end), following extendAlign.cpp:
	match +1, mismatch -1, N in either scores 0; Local keeps the first
	position where the running score reaches a new strict maximum above 0;
	EndToEnd extends over everything. Mismatch limits are not modelled (the
	tests keep ends to a few mismatches)."""
	if extend_to_end:
		return len(read)
	score = best = ext = 0
	for i, (r, g) in enumerate(zip(read, genome)):
		if r == "N" or g == "N":
			continue
		score += 1 if r == g else -1
		if score > best:
			best, ext = score, i + 1
	return ext


def extension_score(read, genome, n):
	"""Score (+1 match, -1 mismatch, N 0) of the first n extended bases."""
	return sum(0 if "N" in (r, g) else (1 if r == g else -1) for r, g in zip(read[:n], genome[:n]))


def extension_past(contig_seq, pos, end, strand, extra):
	"""Bases local_extension() adds to an alignment of genome pos..end
	(1-based) over `extra` bases sequenced after its 3' end."""
	if strand == "+":
		return local_extension(extra, contig_seq[end:end + len(extra)])
	comp = sh.revcomp(extra)[::-1]
	return local_extension(comp, contig_seq[max(0, pos - 1 - len(extra)):pos - 1][::-1])


def end_clips(fwd, ref, mm_left, mm_right, left_to_end=False, right_to_end=False):
	"""(left clip, right clip) of an unspliced read whose genome-orientation
	sequence `fwd` aligns to `ref` (same length) with mismatches only at
	offsets in mm_left (near the left end) and mm_right (near the right end).
	Seeds stop at the innermost end mismatches, so the seed covers
	[max(mm_left)+1, min(mm_right)-1] and each end is extended outward by
	local_extension()."""
	lo = max(mm_left) + 1 if mm_left else 0
	hi = min(mm_right) if mm_right else len(fwd)
	left = lo - local_extension(fwd[:lo][::-1], ref[:lo][::-1], left_to_end)
	right = (len(fwd) - hi) - local_extension(fwd[hi:], ref[hi:], right_to_end)
	return left, right


# STAR's junction motifs on the + strand (donor dinucleotide, acceptor
# dinucleotide) -> (canonSJ code, penalty option); stitchAlignToTranscript.cpp.
MOTIF_PENALTY_KEY = {
	("GT", "AG"): (1, None), ("CT", "AC"): (2, None),
	("GC", "AG"): (3, "GCAG"), ("CT", "GC"): (4, "GCAG"),
	("AT", "AC"): (5, "ATAC"), ("GT", "AT"): (6, "ATAC"),
}


def motif_and_penalty(donor, acceptor, pens):
	"""(canonSJ code, penalty) for an intron starting with `donor` and ending
	with `acceptor` (+ strand); pens = {'Noncan':, 'GCAG':, 'ATAC':}."""
	code, key = MOTIF_PENALTY_KEY.get((donor, acceptor), (0, "Noncan"))
	return code, (0 if key is None else pens[key])


DEFAULT_PENS = {"Noncan": -8, "GCAG": -4, "ATAC": -8}


def junction_placement(contig_seq, read_fwd, n_donor, d, a, shift=1, pens=DEFAULT_PENS):
	"""Junction position STAR picks when stitching a read that joins genome
	bases ..d (1-based, last donor base) to a.. (first acceptor base), with
	n_donor read bases before the join, seeds ending exactly at the join and
	no read gap. Follows the scan in stitchAlignToTranscript.cpp: move the
	junction left while matches lost minus scoreStitchSJshift stay >= 0,
	then scan right over the acceptor piece for the highest (matches gained -
	matches lost + motif penalty), keeping the leftmost on ties; for a
	non-canonical winner flush left over repeats. Returns (jR, motif code,
	mismatches): jR is the junction shift in read bases (0 = at the join,
	-k = k bases to the left), mismatches are read bases that then disagree
	with the genome."""
	G = contig_seq
	R = read_fwd
	rA = n_donor - 1          # read index of the last donor base
	gA = d - 1                # 0-based genome index of the last donor base
	gB1 = a - 2               # 0-based index of the last intron base
	rBend = len(R) - 1
	s1 = 0
	j = 1
	while True:
		j -= 1
		if R[rA + j] != G[gB1 + j] and R[rA + j] == G[gA + j]:
			s1 -= 1
		if not (s1 + shift >= 0 and n_donor + j > 1):
			break
	best = None
	s1 = 0
	while True:
		if R[rA + j] == G[gA + j] and R[rA + j] != G[gB1 + j]:
			s1 += 1
		if R[rA + j] != G[gA + j] and R[rA + j] == G[gB1 + j]:
			s1 -= 1
		code, pen = motif_and_penalty(G[gA + j + 1:gA + j + 3], G[gB1 + j - 1:gB1 + j + 1], pens)
		if best is None or s1 + pen > best[0]:
			best = (s1 + pen, j, code)
		j += 1
		if not j < rBend - rA:
			break
	_, jr, code = best
	if code == 0:
		while G[gA + jr] == G[gB1 + jr]:
			jr -= 1
	mm = 0
	if jr > 0:
		mm = sum(R[rA + i] != G[gA + i] for i in range(1, jr + 1))
	elif jr < 0:
		mm = sum(R[rA + i] != G[gB1 + i] for i in range(jr + 1, 1))
	return jr, code, mm


class MT19937(object):
	"""std::mt19937 seeded with seed(s) (init_genrand), 32-bit outputs."""

	def __init__(self, seed):
		self.mt = [0] * 624
		self.mt[0] = seed & 0xFFFFFFFF
		for i in range(1, 624):
			self.mt[i] = (1812433253 * (self.mt[i - 1] ^ (self.mt[i - 1] >> 30)) + i) & 0xFFFFFFFF
		self.i = 624

	def next32(self):
		if self.i >= 624:
			mt = self.mt
			for k in range(624):
				y = (mt[k] & 0x80000000) | (mt[(k + 1) % 624] & 0x7FFFFFFF)
				mt[k] = mt[(k + 397) % 624] ^ (y >> 1) ^ (0x9908B0DF if y & 1 else 0)
			self.i = 0
		y = self.mt[self.i]
		self.i += 1
		y ^= y >> 11
		y ^= (y << 7) & 0x9D2C5680
		y ^= (y << 15) & 0xEFC60000
		y ^= y >> 18
		return y

	def uniform01(self):
		"""std::uniform_real_distribution<double>(0, 1) as libstdc++ draws it
		from a 32-bit engine: generate_canonical takes two outputs,
		(g1 + g2 * 2^32) / 2^64."""
		g1 = self.next32()
		g2 = self.next32()
		return (g1 + g2 * 4294967296.0) / 18446744073709551616.0


def random_multimapper_order(order, n_best, rng):
	"""The order --outMultimapperOrder Random gives a read's alignments
	(ReadAlign_multMapSelect.cpp): best-scoring alignments first, each group
	shuffled with Fisher-Yates using rand1 = int(u * itr + 0.5)."""
	out = list(order)
	for itr in range(n_best - 1, 0, -1):
		r = int(rng.uniform01() * itr + 0.5)
		out[itr], out[r] = out[r], out[itr]
	for itr in range(len(out) - n_best - 1, 0, -1):
		r = int(rng.uniform01() * itr + 0.5)
		out[n_best + itr], out[n_best + r] = out[n_best + r], out[n_best + itr]
	return out


def tx_junction_read(sim, tx_id, k, left, length=100, sense=True, name=None):
	"""Read over intron k (0-based, genome order) of transcript tx_id with
	`left` of its bases in the exon left of the intron (genome orientation)."""
	t = sim.g.transcripts[tx_id]
	b = sum(e - s + 1 for s, e in t.exons[:k + 1])
	a = b - left
	off = a if t.strand == "+" else t.length - a - length
	return sim.transcript_read(tx_id, off, length, sense=sense, name=name)


def n_head_read(name, synth, contig, pos, k, strand, length=100):
	"""Read whose genome-orientation sequence is k Ns followed by genome
	bases pos..pos+length-k-1: Ns never seed and score 0, so the alignment
	is exactly kS(length-k)M at pos."""
	fwd = "N" * k + synth.seq(contig, pos, pos + length - k - 1)
	return oriented(name, fwd, contig, pos, ("%dS%dM" % (k, length - k)) if k else "%dM" % length, strand)


def unspliced_fallback(contig_seq, fwd, pos1, left, pos2):
	"""Unspliced alignments of a two-block spliced read (genome orientation
	`fwd`; `left` bases at pos1, the rest at pos2) when its junction is not
	allowed: each block locally extended into the intron. Returns the list
	[(score, pos, cigar)] for both blocks, best first; score = aligned
	matches - mismatches + length_score(span)."""
	n = len(fwd)
	right = n - left
	rest = fwd[left:]
	gen = contig_seq[pos1 - 1 + left:pos1 - 1 + n]
	e1 = local_extension(rest, gen)
	a = (left + extension_score(rest, gen, e1) + length_score(left + e1), pos1,
		sh.ops_to_cigar([(left + e1, "M"), (right - e1, "S")]))
	head = fwd[:left][::-1]
	gen = contig_seq[max(0, pos2 - 1 - left):pos2 - 1][::-1]
	e2 = local_extension(head, gen)
	b = (right + extension_score(head, gen, e2) + length_score(right + e2), pos2 - e2,
		sh.ops_to_cigar([(left - e2, "S"), (right + e2, "M")]))
	return sorted([a, b], reverse=True)
