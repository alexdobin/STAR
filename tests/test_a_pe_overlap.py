"""Area A, merging of overlapping mates (ReadAlign_peOverlapMergeMap.cpp):
--peOverlapNbasesMin and --peOverlapMMp."""

import helpers_a as ha
import star_helpers as sh
from star_helpers import ReadSet, record_truth, run_star

TXS = ("T01a", "T02", "T03", "T05", "T08", "T09", "T10", "T06")


def _short_overhang_pairs(sim, kinds=("R", "L"), mm_left=()):
	"""Pairs from 150-base fragments over each intron (no annotation used),
	so the mates overlap by 50 bases. Kind 'R': the right mate (genome order)
	starts 3 bases before the junction (3 + 97 bases), the left mate is
	split 53 + 47. Kind 'L': the left mate ends 3 bases after the junction
	(97 + 3), the right mate is split 47 + 53. mm_left: substitutions at
	these genome-orientation offsets of the left mate. extra['short'] names
	the mate (1 or 2) with the 3-base overhang."""
	pairs = []
	for tx in TXS:
		t = sim.g.transcripts[tx]
		lens = [e - s + 1 for s, e in t.exons]
		for k in range(len(lens) - 1):
			b = sum(lens[:k + 1])
			for kind in kinds:
				a = b - 53 if kind == "R" else b - 97
				if a < b - lens[k] or a + 150 > b + lens[k + 1]:
					continue
				off = a if t.strand == "+" else t.length - a - 150
				forward = sim.rng.random() < 0.5
				# the left mate (genome order) is mate 1 when mate 1 is on the
				# transcript's + strand
				left_is_1 = forward == (t.strand == "+")
				mm1, mm2 = (mm_left, ()) if left_is_1 else ((), mm_left)
				p = sim.transcript_pair(tx, off, 150, forward=forward, mismatches1=mm1, mismatches2=mm2,
					name="%s_%d_%s" % (tx, k, kind))
				p.extra["short"] = (2 if left_is_1 else 1) if kind == "R" else (1 if left_is_1 else 2)
				pairs.append(p)
	return ReadSet(pairs)


def _unmerged(synth, mate):
	"""(pos, cigar) of the 3-base-overhang mate aligned on its own: the
	97-base block extended into the intron by the local rule, the rest
	soft-clipped."""
	(s1, e1), (s2, e2) = sh.cigar_blocks(mate.pos, mate.cigar)
	_, pos, cig = ha.unspliced_fallback(synth.contigs[mate.contig], mate.sam_seq(), s1, e1 - s1 + 1, s2)[0]
	return pos, cig


def _check(aln, pairs, synth, merged):
	by_name = aln.by_name()
	assert sorted(by_name) == sorted(pairs.by_name)
	unspliced = 0
	for name, recs in by_name.items():
		p = pairs.by_name[name]
		got = {r.mate: r for r in recs}
		assert sorted(got) == [1, 2]
		for i, truth in ((1, p.r1), (2, p.r2)):
			rec = got[i]
			if merged or i != p.extra.get("short"):
				assert record_truth(rec) == truth.truth(), (name, i)
			else:
				pos, cig = _unmerged(synth, truth)
				assert (rec.rname, rec.pos, rec.cigar, rec.strand) == (truth.contig, pos, cig, truth.strand), name
				unspliced += "N" not in cig
	return unspliced


def test_overlap_merging_rescues_short_overhangs(genome_plain, sim):
	"""Without merging (default --peOverlapNbasesMin 0) the 3-base overhang
	of one mate is below what seeding and --alignSJoverhangMin can splice,
	so that mate aligns unspliced (its 97-base block plus local extension);
	the other mate is spliced. With --peOverlapNbasesMin 10 the mates are
	merged into one 150-base read over their 50-base overlap, the merged
	read is spliced with long blocks on both sides, and both mates get their
	true spliced CIGARs. Non-overlapping pairs (insert 300) are unchanged."""
	pairs = _short_overhang_pairs(sim)
	controls = ReadSet([sim.genomic_pair(*sim.unique_locus(400), insert=300, strand="+-"[i % 2], name="ctrl%d" % i)
		for i in range(4)])
	assert len(pairs) >= 20
	aln = run_star(genome_plain, pairs + controls).alignments()
	assert _check(sh.AlignmentFile(aln.header, aln.references, [r for r in aln if not r.qname.startswith("ctrl")]),
		pairs, sim.g, merged=False) == len(pairs)
	aln = run_star(genome_plain, pairs + controls, "--peOverlapNbasesMin 10").alignments()
	_check(sh.AlignmentFile(aln.header, aln.references, [r for r in aln if not r.qname.startswith("ctrl")]),
		pairs, sim.g, merged=True)
	ctrl = sh.AlignmentFile(aln.header, aln.references, [r for r in aln if r.qname.startswith("ctrl")]).by_name()
	assert sorted(ctrl) == sorted(controls.by_name)
	for name, recs in ctrl.items():
		p = controls.by_name[name]
		assert sorted(record_truth(r) for r in recs) == sorted([p.r1.truth(), p.r2.truth()])


def test_overlap_mismatch_proportion(genome_plain, sim):
	"""The mates are merged only if their overlap has at most
	--peOverlapMMp mismatches per matched base (localSearchNisMM). The left
	mate carries substitutions at offsets 70 and 85, inside the 50-base
	overlap: 2/48 = 0.042 exceeds the default 0.01, so no merge happens and
	the right mate keeps its unspliced alignment; with --peOverlapMMp 0.05
	the mates merge and both get their true CIGARs (the left mate with nM
	2)."""
	pairs = _short_overhang_pairs(sim, kinds=("R",), mm_left=(70, 85))
	assert len(pairs) >= 10
	aln = run_star(genome_plain, pairs, "--peOverlapNbasesMin 10").alignments()
	_check(aln, pairs, sim.g, merged=False)
	aln = run_star(genome_plain, pairs, "--peOverlapNbasesMin 10 --peOverlapMMp 0.05").alignments()
	_check(aln, pairs, sim.g, merged=True)
	assert all(r.tag("nM") == 2 for r in aln)
