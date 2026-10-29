"""Helpers for the area-D tests (STARsolo and chimeric detection).

Read builders that the shared harness does not have (10x pairs around an
arbitrary cDNA read, chimeric reads whose segments cannot be extended past
the breakpoint, chimeric pairs), and reference implementations of the
documented STARsolo rules (UMI collapsing, multi-gene UMI distribution, cell
calling) that the tests compare STAR's matrices with.
"""

import itertools
import zlib
from collections import OrderedDict

import star_helpers as sh
from star_helpers import ReadSet, SimPair, SimRead

CB_LEN = 16
UMI_LEN = 12
TENX_GEOMETRY = ["--soloCBstart", "1", "--soloCBlen", str(CB_LEN), "--soloUMIstart", "17", "--soloUMIlen", str(UMI_LEN)]


def solo_args(whitelist, *extra):
	"""CB_UMI_Simple in 10x 3' v3 geometry (16 bp CB at 1, 12 bp UMI at 17)."""
	return ["--soloType", "CB_UMI_Simple", "--soloCBwhitelist", str(whitelist)] + TENX_GEOMETRY + [str(x) for x in extra]


def module_sim(synth, module, name, **kw):
	"""A ReadSimulator for a module-scoped fixture, seeded like conftest.sim."""
	seed = zlib.crc32(("%s::%s" % (module, name)).encode())
	return sh.ReadSimulator(synth, seed=seed, **kw)


# Genes added to the synthetic annotation for the area-D genome (1-based,
# inclusive, single transcript each): GN1/GN2 nest inside G09's intron
# (chr2:12301-17300) on the same and the opposite strand; GR1/GR2 cover the
# two copies of the exact repeat dup2 (chr1:26001-27500 +, chr2:8001-9500
# reverse-complemented) plus 200 unique bases each, so dup2 reads align to
# two loci that belong to different genes.
EXTRA_GENES = (
	("GN1", "GeneNestedSense", "chr2", "+", [(13001, 13200), (13401, 13600)]),
	("GN2", "GeneNestedAnti", "chr2", "-", [(16001, 16400)]),
	("GR1", "GeneRepeatOne", "chr1", "+", [(25801, 27500)]),
	("GR2", "GeneRepeatTwo", "chr2", "-", [(8001, 9700)]),
)


def extra_genes_gtf():
	"""GTF lines for EXTRA_GENES in the format of SyntheticGenome.gtf_text."""
	lines = []
	for gid, name, contig, strand, exons in EXTRA_GENES:
		ga = 'gene_id "%s"; gene_name "%s"; gene_biotype "protein_coding";' % (gid, name)
		ta = 'gene_id "%s"; transcript_id "%s"; gene_name "%s"; gene_biotype "protein_coding";' % (gid, "T" + gid, name)
		s, e = exons[0][0], exons[-1][1]
		lines.append("\t".join([contig, "synthetic", "gene", str(s), str(e), ".", strand, ".", ga]))
		lines.append("\t".join([contig, "synthetic", "transcript", str(s), str(e), ".", strand, ".", ta]))
		ordered = exons if strand == "+" else exons[::-1]
		for i, (a, b) in enumerate(ordered, 1):
			lines.append("\t".join([contig, "synthetic", "exon", str(a), str(b), ".", strand, ".",
				ta + ' exon_number "%d";' % i]))
	return "\n".join(lines) + "\n"


def area_d_genome(genome_factory, synth):
	"""STAR index of the synthetic genome annotated with its GTF plus
	EXTRA_GENES (built once per session: the GTF path is fixed)."""
	from pathlib import Path
	gtf = Path(synth.gtf_path).parent / "annotation_area_d.gtf"
	text = synth.gtf_text() + extra_genes_gtf()
	if not gtf.exists() or gtf.read_text() != text:
		gtf.write_text(text)
	return genome_factory(gtf=str(gtf))


def write_whitelist(path, barcodes):
	path.write_text("\n".join(barcodes) + "\n")
	return path


def substitute(seq, offset, base=None):
	"""seq with the base at `offset` replaced (by `base`, else by the next
	base in ACGT order)."""
	if base is None:
		base = "ACGT"["ACGT".index(seq[offset]) + 1 & 3]
	return seq[:offset] + base + seq[offset + 1:]


class UmiPool(object):
	"""Distinct UMIs, pairwise >= 3 mismatches apart and never homopolymers,
	drawn from the simulator's generator."""

	def __init__(self, sim, length=UMI_LEN):
		self.sim = sim
		self.length = length
		self.used = []

	def take(self):
		while True:
			u = self.sim.random_seq(self.length)
			if len(set(u)) > 1 and all(sh.hamming(u, v) >= 3 for v in self.used):
				self.used.append(u)
				return u


def tenx_pair(sim, cdna, cb, umi, cb_read=None, umi_read=None, bc_qual=None, gene=None):
	"""10x pair around an existing cDNA SimRead (renamed to the pair name).
	cb/umi are the truth, cb_read/umi_read what is sequenced (default: the
	truth). extra holds cb, umi, cb_read, umi_read and gene."""
	name = cdna.name
	cb_read = cb if cb_read is None else cb_read
	umi_read = umi if umi_read is None else umi_read
	bc = cb_read + umi_read
	cdna.mate = 1
	r2 = SimRead(name, bc, bc_qual if bc_qual is not None else sim.qual(len(bc)), kind="barcode")
	r2.mate = 2
	return SimPair(name, cdna, r2, None, "solo",
		extra={"cb": cb, "umi": umi, "cb_read": cb_read, "umi_read": umi_read, "gene": gene})


def exonic_read(sim, tx_id, length=90, sense=True, from_end=150):
	"""Unspliced read inside the last exon of tx_id (in transcript order),
	`from_end` bases before the transcript's 3' end."""
	t = sim.g.transcripts[tx_id]
	return sim.transcript_read(tx_id, t.length - from_end, length, sense=sense, name=sim.name("cell"))


def genomic_read(sim, contig, pos, strand, length=90):
	return sim.genomic_read(contig, pos, strand, length, name=sim.name("cell"))


# ------------------------------------------------------------- UMI collapsing

def _one_mm(a, b):
	return sh.hamming(a, b) == 1


def dedup_exact(counts):
	"""Exact: distinct UMI sequences (parametersDefault soloUMIdedup)."""
	return len(counts)


def dedup_nodedup(counts):
	"""NoDedup: every read is counted."""
	return sum(counts.values())


def dedup_1mm_all(counts):
	"""1MM_All: UMIs joined by a chain of 1-mismatch neighbours are counted
	once, i.e. the number of connected components of the 1MM graph."""
	umis = list(counts)
	parent = {u: u for u in umis}

	def find(u):
		while parent[u] != u:
			u = parent[u]
		return u

	for a, b in itertools.combinations(umis, 2):
		if _one_mm(a, b):
			parent[find(a)] = find(b)
	return len({find(u) for u in umis})


def _directional_once(order, counts, add):
	cluster = {}
	for i, b in enumerate(order):
		cluster[b] = b
		for a in order[:i]:
			if _one_mm(a, b) and counts[a] >= 2 * counts[b] + add:
				cluster[b] = cluster[a]
				break
	return len(set(cluster.values()))


def dedup_directional(counts, add):
	"""UMI-tools "directional" (Smith et al. 2017): taking UMIs from the
	highest count down, a UMI b joins the cluster of the first higher-ranked
	1-mismatch neighbour a with count(a) >= 2*count(b) + add (add=-1 for
	1MM_Directional_UMItools, 0 for the stricter 1MM_Directional). Returns
	the number of clusters; raises if the answer depends on how UMIs with
	equal counts are ordered, so test data must avoid such ties."""
	results = set()
	for rev in (False, True):
		order = sorted(counts, key=lambda u: (-counts[u], u), reverse=False)
		if rev:
			order = sorted(counts, key=lambda u: (-counts[u], tuple(-ord(c) for c in u)))
		results.add(_directional_once(order, counts, add))
	if len(results) != 1:
		raise ValueError("directional result depends on tie order: %r" % counts)
	return results.pop()


def cr_corrections(counts):
	"""CellRanger 2-4 UMI correction: each UMI is replaced by its 1-mismatch
	neighbour with the largest (count, sequence), if that is larger than its
	own (count, sequence); one pass, not transitive."""
	out = {}
	for u in counts:
		best = (counts[u], u)
		for v in counts:
			if _one_mm(u, v) and (counts[v], v) > best:
				best = (counts[v], v)
		out[u] = best[1]
	return out


def dedup_1mm_cr(counts):
	"""1MM_CR: number of distinct UMIs after the CellRanger correction."""
	return len(set(cr_corrections(counts).values()))


DEDUP = OrderedDict([
	("1MM_All", dedup_1mm_all),
	("1MM_Directional_UMItools", lambda c: dedup_directional(c, -1)),
	("1MM_Directional", lambda c: dedup_directional(c, 0)),
	("Exact", dedup_exact),
	("NoDedup", dedup_nodedup),
	("1MM_CR", dedup_1mm_cr),
])


def umi_read_counts(readset, key=lambda p: (p.extra["gene"], p.extra["cb"])):
	"""{key: {sequenced UMI: number of reads}} from 10x truth."""
	out = OrderedDict()
	for p in readset:
		d = out.setdefault(key(p), OrderedDict())
		d[p.extra["umi_read"]] = d.get(p.extra["umi_read"], 0) + 1
	return out


# --------------------------------------------------------- multi-gene UMIs

def distribute_multigene(unique, multi, method):
	"""Fractional counts that one cell's multi-gene UMIs add to each gene,
	following docs/STARsolo.md "Multi-gene reads".

	unique: {gene: unique-gene UMI count}; multi: list of gene tuples, one
	per multi-gene UMI. Returns {gene: added count} for every gene that
	appears in a multi-gene UMI.
	- Uniform: each UMI adds 1/N to each of its N genes.
	- PropUnique: proportional to the unique counts, uniform if all are 0.
	- Rescue: proportional to unique + Uniform counts (first EM step).
	- EM: maximum likelihood by expectation-maximisation, started from
	  unique + Uniform, iterated until no gene changes by >= 0.01 (STAR's
	  stopping rule, SoloFeature_collapseUMIall.cpp; counts < 0.01 are
	  zeroed before each step).
	"""
	genes = []
	for gs in multi:
		for g in gs:
			if g not in genes:
				genes.append(g)
	uni = {g: 0.0 for g in genes}
	for gs in multi:
		for g in gs:
			uni[g] += 1.0 / len(gs)
	if method == "Uniform":
		return uni
	u = {g: float(unique.get(g, 0)) for g in genes}
	out = {g: 0.0 for g in genes}
	if method == "PropUnique":
		for gs in multi:
			norm = sum(u[g] for g in gs)
			for g in gs:
				out[g] += 1.0 / len(gs) if norm == 0 else u[g] / norm
		return out
	if method == "Rescue":
		for gs in multi:
			norm = sum(uni[g] + u[g] for g in gs)
			for g in gs:
				out[g] += (uni[g] + u[g]) / norm
		return out
	if method == "EM":
		old = {g: uni[g] + u[g] for g in genes}
		for _ in range(102):
			for g in genes:
				if old[g] < 0.01:
					old[g] = 0.0
			new = dict(u)
			for gs in multi:
				norm = sum(old[g] for g in gs)
				for g in gs:
					new[g] += old[g] / norm
			change = max(abs(new[g] - old[g]) for g in genes)
			old = new
			if change < 0.01:
				break
		return {g: old[g] - u[g] for g in genes}
	raise ValueError(method)


# -------------------------------------------------------------- cell calling

def knee_filter(totals, n_expected, max_percentile, max_min_ratio):
	"""CellRanger 2.2 cell calling (docs/STARsolo.md "Knee filtering"): the
	robust maximum is the UMI total of the cell ranked
	round(n_expected * (1 - max_percentile)) (0-based, capped at the last
	cell); cells with at least max(1, round(robust max / max_min_ratio))
	UMIs are kept. Avoid .5 in both roundings (Python rounds half to even)."""
	s = sorted((t for t in totals.values() if t > 0), reverse=True)
	idx = min(len(s) - 1, int(round(n_expected * (1 - max_percentile))))
	cutoff = max(1, int(round(s[idx] / float(max_min_ratio))))
	return {cb for cb, t in totals.items() if t >= cutoff}


# ------------------------------------------------------------------ chimeras

def segment_locus(sim, segment, length, copy=0, margin=30):
	"""(contig, pos) of a window of `length` bases inside copy `copy` of a
	repeated segment, at least `margin` bases from its ends."""
	c = sim.g.segments[segment].copies[copy]
	return c.contig, sim.rng.randint(c.start + margin, c.end - margin - length + 1)


def other_copy(sim, segment, contig, pos, strand, length):
	"""(contig, pos, strand) of the same sequence in each other copy of an
	exact repeated segment, for a piece at contig:pos on `strand` in one copy."""
	seg = sim.g.segments[segment]
	src = next(c for c in seg.copies if c.contig == contig and c.start <= pos <= c.end)
	off = pos - src.start if src.strand == "+" else src.end - (pos + length - 1)
	out = []
	for c in seg.copies:
		if c is src:
			continue
		p = c.start + off if c.strand == "+" else c.end - off - length + 1
		s = strand if c.strand == src.strand else ("-" if strand == "+" else "+")
		out.append((c.contig, p, s))
	return out


def piece_context(g, piece, flank):
	"""(sequence before, piece, sequence after) of piece=(contig, pos,
	strand, n) in the piece's own orientation, `flank` bases each side."""
	contig, pos, strand, n = piece
	s = g.contigs[contig]
	left = s[max(0, pos - 1 - flank):pos - 1]
	seq = s[pos - 1:pos - 1 + n]
	right = s[pos - 1 + n:pos - 1 + n + flank]
	if strand == "+":
		return left, seq, right
	return sh.revcomp(right), sh.revcomp(seq), sh.revcomp(left)


def junction_is_clean(sim, a, b, motif=0, flank=6):
	"""True when the breakpoint between consecutive pieces a and b cannot
	move: repeats (0, 0), junction type `motif` there and no GT/AG or CT/AC
	within 8 bases otherwise, and the `flank` bases of each piece's genomic
	continuation all differ from the read, so neither segment's alignment
	extends past the breakpoint."""
	r = sim.chimeric_read(a, b, name="probe")
	if r.extra["repeat"] != (0, 0) or r.extra["motif"] != motif:
		return False
	if any(v for k, v in r.extra["shift_motifs"].items() if k != 0):
		return False
	_, seq_a, after_a = piece_context(sim.g, a, flank)
	before_b, seq_b, _ = piece_context(sim.g, b, flank)
	return (len(after_a) == len(before_b) == flank
		and all(x != y for x, y in zip(after_a, seq_b[:flank]))
		and all(x != y for x, y in zip(before_b[::-1], seq_a[::-1][:flank])))


def _motif_shift(g, piece, side, dinuc, span):
	"""Shift (0..span) that makes the bases after (side='after') or before
	(side='before') the piece read `dinuc` in its orientation, or None."""
	contig, pos, strand, n = piece
	for d in range(span + 1):
		p = (contig, pos + d, strand, n)
		before, _, after = piece_context(g, p, 2)
		if (after if side == "after" else before) == dinuc:
			return d
	return None


def strict_chimeric_read(sim, strands=("+", "+"), length=(50, 50), same_contig=False, min_dist=5000,
		flank=6, motif=0, name=None, locus_a=None, locus_b=None, max_tries=20000):
	"""Two-piece chimeric read (sim.chimeric_read) whose breakpoint is clean
	(junction_is_clean), so Chimeric.out.junction must equal
	extra['junction'] and each segment's block is exactly its piece. Pieces
	are unique and gene-free unless locus_a/locus_b (callables returning
	(contig, pos)) place them elsewhere; same_contig puts both on one contig
	at least min_dist apart. motif 1 (GT/AG) or 2 (CT/AC) asks for that
	junction type."""
	dinuc = {1: ("GT", "AG"), 2: ("CT", "AC")}.get(motif)
	span = 300 if dinuc else 0
	for _ in range(max_tries):
		if locus_a:
			ca, pa = locus_a()
		else:
			ca, pa = sim.unique_locus(length[0] + 2 * flank + span, avoid_genes=True)
			pa += flank
		if locus_b:
			cb, pb = locus_b()
		else:
			cb, pb = sim.unique_locus(length[1] + 2 * flank + span, avoid_genes=True)
			pb += flank
		if same_contig != (ca == cb):
			continue
		if same_contig and abs(pa - pb) < min_dist:
			continue
		a = (ca, pa, strands[0], length[0])
		b = (cb, pb, strands[1], length[1])
		if dinuc:
			da = _motif_shift(sim.g, a, "after", dinuc[0], span)
			db = _motif_shift(sim.g, b, "before", dinuc[1], span)
			if da is None or db is None:
				continue
			a = (ca, pa + da, a[2], a[3])
			b = (cb, pb + db, b[2], b[3])
		if junction_is_clean(sim, a, b, motif, flank):
			return sim.chimeric_read(a, b, name=name)
	raise RuntimeError("no chimeric read found")


def multi_piece_read(sim, pieces, name=None, flank=6):
	"""Read joining several genomic pieces (contig, pos, strand, n) in
	sequencing order, or None unless every consecutive breakpoint is clean
	(non-canonical). extra['pieces'] lists them."""
	for a, b in zip(pieces, pieces[1:]):
		if not junction_is_clean(sim, a, b, 0, flank):
			return None
	seq = "".join(piece_context(sim.g, p, 0)[1] for p in pieces)
	c, p, s, n = pieces[0]
	return SimRead(name or sim.name(), seq, sim.qual(len(seq)), c, s, p, None, kind="chimeric",
		extra={"pieces": list(pieces)})


def junction_columns(read, i, motif=0):
	"""Expected first 14 Chimeric.out.junction columns (manual 6.4) for the
	chimera joining pieces i (donor) and i+1 (acceptor) of a multi-piece read:
	breakpoints are the first intron base on each side, CIGARs are on the +
	strand with every base outside the segment soft-clipped."""
	pieces = read.extra["pieces"]
	L = sum(p[3] for p in pieces)
	starts = [sum(p[3] for p in pieces[:k]) for k in range(len(pieces))]

	def cigar(k):
		s, e, strand = starts[k], starts[k] + pieces[k][3], pieces[k][2]
		ops = [(s, "S"), (e - s, "M"), (L - e, "S")]
		return sh.ops_to_cigar(ops if strand == "+" else ops[::-1])

	(ca, pa, sa, na), (cb, pb, sb, nb) = pieces[i], pieces[i + 1]
	return OrderedDict([
		("chr_donorA", ca), ("brkpt_donorA", pa + na if sa == "+" else pa - 1), ("strand_donorA", sa),
		("chr_acceptorB", cb), ("brkpt_acceptorB", pb - 1 if sb == "+" else pb + nb), ("strand_acceptorB", sb),
		("junction_type", motif), ("repeat_left_lenA", 0), ("repeat_right_lenB", 0), ("read_name", read.name),
		("start_alnA", pa), ("cigar_alnA", cigar(i)), ("start_alnB", pb), ("cigar_alnB", cigar(i + 1)),
	])


def shifted_junction(read, shift):
	"""First 14 Chimeric.out.junction columns for a two-piece read
	(sim.chimeric_read) whose breakpoint is moved by `shift` read bases
	(> 0: piece A longer, B shorter), with the repeat lengths measured from
	the moved breakpoint: left = repeat[0] + shift, right = repeat[1] - shift."""
	(ca, pa, sa, na), (cb, pb, sb, nb) = read.extra["pieces"]
	rl, rr = read.extra["repeat"]
	la, lb = na + shift, nb - shift
	if sa == "+":
		start_a, brk_a, cig_a = pa, pa + la, "%dM%dS" % (la, lb)
	else:
		start_a, brk_a, cig_a = pa - shift, pa - shift - 1, "%dS%dM" % (lb, la)
	if sb == "+":
		start_b, brk_b, cig_b = pb + shift, pb + shift - 1, "%dS%dM" % (la, lb)
	else:
		start_b, brk_b, cig_b = pb, pb + nb - shift, "%dM%dS" % (lb, la)
	return OrderedDict([
		("chr_donorA", ca), ("brkpt_donorA", brk_a), ("strand_donorA", sa),
		("chr_acceptorB", cb), ("brkpt_acceptorB", brk_b), ("strand_acceptorB", sb),
		("junction_type", 0), ("repeat_left_lenA", rl + shift), ("repeat_right_lenB", rr - shift),
		("read_name", read.name), ("start_alnA", start_a), ("cigar_alnA", cig_a),
		("start_alnB", start_b), ("cigar_alnB", cig_b),
	])


def chimeric_pair(sim, a, b, length=100, name=None):
	"""Pair whose mate 1 is the `length` bases of genome at a=(contig, pos,
	strand) read on that strand and whose mate 2 is the reverse complement
	of the `length` bases at b, so the chimeric junction lies between the
	mates. extra['junction'] holds the expected first 14 columns of
	Chimeric.out.junction (manual 6.4, encompassing junction: type -1, the
	base after mate 1's 3' end and the base before mate 2's 5' end, both on
	their own strands)."""
	(ca, pa, sa), (cb, pb, sb) = a, b
	name = name or sim.name("pair")
	r1 = sim.genomic_read(ca, pa, sa, length, name=name)
	piece_b = sim.genomic_read(cb, pb, sb, length, name=name)
	r2 = SimRead(name, sh.revcomp(piece_b.seq), piece_b.qual[::-1], cb, "-" if sb == "+" else "+", pb,
		"%dM" % length, kind="genomic")
	r1.mate, r2.mate = 1, 2
	junction = OrderedDict([
		("chr_donorA", ca), ("brkpt_donorA", pa + length if sa == "+" else pa - 1), ("strand_donorA", sa),
		("chr_acceptorB", cb), ("brkpt_acceptorB", pb - 1 if sb == "+" else pb + length), ("strand_acceptorB", sb),
		("junction_type", -1), ("repeat_left_lenA", 0), ("repeat_right_lenB", 0), ("read_name", name),
	])
	return SimPair(name, r1, r2, None, "chimeric_pair", extra={"junction": junction, "pieces": [a, b]})
