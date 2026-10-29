"""Helpers for the area C tests (SJ.out.tab, sjdb insertion, 2-pass mapping,
GeneCounts, TranscriptomeSAM, Log.final.out).

Coordinates are 1-based and inclusive, as in star_helpers.
"""

from collections import Counter

import star_helpers as sh

# Gap penalties by motif code with STAR's defaults (parametersDefault:
# scoreGap 0, scoreGapNoncan -8, scoreGapGCAG -4, scoreGapATAC -8).
GAP_PENALTY = {0: -8, 1: 0, 2: 0, 3: -4, 4: -4, 5: -8, 6: -8}


def spliced_read(sim, contig, blocks, strand="+", mismatches=(), name=None, kind="spliced"):
	"""Read whose aligned bases are the genomic blocks [(start, end), ...]
	(ascending), sequenced from `strand`. Gaps between blocks become N in the
	truth CIGAR; mismatches are read offsets from the leftmost base."""
	path = []
	for s, e in blocks:
		path.extend(range(s, e + 1))
	seq, pos, cig, mm = sim._build(contig, path, 0, len(path), strand, mismatches)
	return sim._make(name or sim.name(), seq, contig, strand, pos, cig, mm, kind)


def junction_read(sim, contig, intron, left, right, strand="+", name=None):
	"""Read crossing `intron` (start, end) with `left` exonic bases before
	it and `right` after it."""
	s, e = intron
	return spliced_read(sim, contig, [(s - left, s - 1), (e + 1, e + right)], strand, name=name)


def motif_of(synth, contig, start, end):
	"""STAR motif code (0-6) of the intron contig:start-end."""
	return sh.motif_code(synth.seq(contig, start, start + 1), synth.seq(contig, end - 1, end))


def junction_is_clean(synth, contig, start, end, left=60, right=60, max_mm=1):
	"""True when, for a read with `left` and `right` bases around the intron,
	every other placement of a junction of the same length with at most
	`max_mm` mismatches scores strictly worse than the true one (each
	mismatch costs 2, plus the motif penalty), so the junction can neither
	slide nor move to a better motif. STAR only searches shifts that cost
	little (parametersDefault: scoreStitchSJshift 1, "maximum score reduction
	while searching for SJ boundaries"), hence the default max_mm=1."""
	g = synth.contigs[contig]
	read = g[start - 1 - left:start - 1] + g[end:end + right]
	gap = end - start + 1
	true_score = GAP_PENALTY[motif_of(synth, contig, start, end)]
	for k in range(1, left + right):
		if k == left:
			continue
		ref = g[start - 1 - left:start - 1 - left + k] + g[start - 1 - left + k + gap:start - 1 - left + gap + left + right]
		mm = sum(1 for a, b in zip(read, ref) if a != b)
		if mm > max_mm:
			continue
		s2 = start - left + k
		alt = GAP_PENALTY[motif_of(synth, contig, s2, s2 + gap - 1)] - 2 * mm
		if alt >= true_score:
			return False
	return True


def extension_length(read_part, ref_part):
	"""Bases a local alignment extends over `ref_part` (match +1, mismatch
	-1): the shortest prefix reaching the highest positive running score, 0
	when no prefix scores above 0. Both sequences start at the base next to
	the aligned block and run away from it."""
	best = best_len = score = 0
	for i, (a, b) in enumerate(zip(read_part, ref_part)):
		score += 1 if a == b else -1
		if score > best:
			best, best_len = score, i + 1
	return best_len


def segment_intron(synth, segment, motif=("GT", "AG"), length=(200, 600), flank=60, avoid=()):
	"""(a, b): 0-based offsets, in the segment's own orientation, of the first
	and last intron base of an intron inside `segment` whose first and last
	two bases are `motif` (segment orientation), with `flank` bases of
	segment on each side, that is clean (junction_is_clean) in every copy.
	The first such intron by donor offset that is not in `avoid`."""
	seg = synth.segments[segment]
	c0 = seg.copies[0]
	assert c0.strand == "+"
	s = synth.seq(c0.contig, c0.start, c0.end)
	for a in range(flank, seg.length - flank - length[0]):
		if s[a:a + 2] != motif[0]:
			continue
		for b in range(a + length[0] - 1, min(a + length[1], seg.length - flank)):
			if s[b - 1:b + 1] != motif[1] or (a, b) in avoid:
				continue
			if all(junction_is_clean(synth, c, i0, i1, flank, flank, max_mm=3)
					for c, i0, i1 in segment_intron_loci(synth, segment, a, b)):
				return a, b
	raise ValueError("no clean intron in %s" % segment)


def segment_intron_loci(synth, segment, a, b):
	"""[(contig, start, end)] of segment-offset intron a..b in every copy."""
	out = []
	for c in synth.segments[segment].copies:
		if c.strand == "+":
			out.append((c.contig, c.start + a, c.start + b))
		else:
			out.append((c.contig, c.end - b, c.end - a))
	return out


def segment_junction_read(sim, segment, a, b, left, right, copy=0, strand="+", name=None):
	"""Read from copy `copy` of `segment` crossing the segment intron a..b
	with `left`/`right` bases before/after it in segment orientation. Truth
	is the read's copy; extra['loci'] lists (contig, pos, cigar, strand) of
	the same alignment in every copy."""
	seg = sim.g.segments[segment]
	loci = []
	for c in seg.copies:
		if c.strand == "+":
			blocks = [(c.start + a - left, c.start + a - 1), (c.start + b + 1, c.start + b + right)]
			st = strand
		else:
			blocks = [(c.end - b - right, c.end - b - 1), (c.end - a + 1, c.end - a + left)]
			st = "-" if strand == "+" else "+"
		loci.append((c.contig, blocks, st))
	contig, blocks, st = loci[copy]
	r = spliced_read(sim, contig, blocks, st, name=name, kind="multi")
	r.extra["loci"] = []
	for contig, blocks, st in loci:
		pos = blocks[0][0]
		cig = "%dM%dN%dM" % (blocks[0][1] - blocks[0][0] + 1, blocks[1][0] - blocks[0][1] - 1, blocks[1][1] - blocks[1][0] + 1)
		r.extra["loci"].append((contig, pos, cig, st))
	r.extra["n_loci"] = len(loci)
	return r


def expected_sj(alignments, synth, annotated=()):
	"""SJRecord rows expected from (contig, pos, cigar, unique) alignments,
	one tuple per alignment (a multimapper contributes one per locus, with
	unique False). annotated: set of (contig, start, end) in the sjdb.
	Column 9 is the largest min(left block, right block) over the
	alignments crossing the junction; rows are in genome order."""
	uniq = Counter()
	multi = Counter()
	overhang = {}
	for contig, pos, cigar, unique in alignments:
		blocks = sh.cigar_blocks(pos, cigar)
		introns = [x for x in sh.cigar_introns(pos, cigar)]
		for i, (s, e) in enumerate(introns):
			key = (contig, s, e)
			(uniq if unique else multi)[key] += 1
			left = blocks[i][1] - blocks[i][0] + 1
			right = blocks[i + 1][1] - blocks[i + 1][0] + 1
			overhang[key] = max(overhang.get(key, 0), min(left, right))
	order = list(synth.contigs)
	rows = []
	for key in sorted(overhang, key=lambda k: (order.index(k[0]), k[1], k[2])):
		contig, s, e = key
		m = motif_of(synth, contig, s, e)
		rows.append(sh.SJRecord(contig, s, e, sh.MOTIF_STRAND[m], m, int(key in annotated), uniq[key],
			multi[key], overhang[key]))
	return rows


def truth_alignments(reads):
	"""(contig, pos, cigar, unique) for SimReads: every locus of a
	multimapper built by segment_junction_read, else the read's truth."""
	out = []
	for r in reads:
		loci = r.extra.get("loci")
		if loci and isinstance(loci[0][2], str) and "N" in loci[0][2]:
			out.extend((c, p, cig, False) for c, p, cig, _ in loci)
		else:
			out.append((r.contig, r.pos, r.cigar, True))
	return out


def annotated_junctions(synth):
	"""{(contig, start, end)} of the introns in the GTF."""
	return {(j.contig, j.start, j.end) for j in synth.junctions()}


def write_gtf(path, lines):
	"""Write GTF lines (lists of 9 fields or strings) to `path`."""
	with open(str(path), "w") as f:
		for line in lines:
			f.write((line if isinstance(line, str) else "\t".join(str(x) for x in line)) + "\n")
	return path


def tx_offset_of(synth, tx_id, genomic_pos):
	"""0-based offset, 5' to 3' along transcript `tx_id`, of a genomic
	position inside one of its exons."""
	t = synth.transcripts[tx_id]
	return t.sense_positions().index(genomic_pos)


def fmt2(x):
	"""A value as STAR's Log.final.out prints it (fixed, 2 decimals), back
	to float for comparison with parse_log_final."""
	return float("%.2f" % x)


def bgzf_block_types(path):
	"""Deflate BTYPE (0 stored, 1 fixed, 2 dynamic Huffman) of the first
	deflate block of every BGZF member with a non-empty payload."""
	data = open(str(path), "rb").read()
	out = []
	off = 0
	while off < len(data):
		assert data[off:off + 2] == b"\x1f\x8b", "not a gzip member at %d" % off
		xlen = data[off + 10] | (data[off + 11] << 8)
		bsize = None
		p = off + 12
		while p < off + 12 + xlen:
			si1, si2 = data[p], data[p + 1]
			slen = data[p + 2] | (data[p + 3] << 8)
			if si1 == 66 and si2 == 67:
				bsize = data[p + 4] | (data[p + 5] << 8)
			p += 4 + slen
		assert bsize is not None, "gzip member without a BGZF BC field"
		payload = data[off + 12 + xlen:off + bsize + 1 - 8]
		isize = int.from_bytes(data[off + bsize + 1 - 4:off + bsize + 1], "little")
		if isize > 0:
			out.append((payload[0] >> 1) & 3)
		off += bsize + 1
	return out
