"""Area C: every count field of Log.final.out, worked out by hand from the
simulated read set.

Manual 5.1: statistics are per read; most are collected over unique mappers;
mismatch/indel rates are the number of mismatches/indels in unique mappers
divided by their mapped bases. Values print with two decimals.
"""

import star_helpers as sh
from star_helpers import ReadSet, run_star

import helpers_c as hc


def _unique_genomic(sim, n, length=None, strand=None, **kw):
	"""n unspliced reads at unique loci outside genes."""
	out = []
	for i in range(n):
		c, p = sim.unique_locus(length=(length or 100) + 10, avoid_genes=True)
		out.append(sim.genomic_read(c, p, strand or "+-"[i % 2], length=length, **kw))
	return out


def _deletion_read(sim, n):
	while True:
		c, p = sim.unique_locus(length=120, avoid_genes=True)
		if sim.deletion_is_unambiguous(c, p + 50, n):
			return sim.genomic_read(c, p, "+", deletion=(50, n))


def test_log_final_every_field(genome_sjdb, sim):
	"""86 reads on the GTF index, default parameters:
	- 30 unspliced 100-base and 10 unspliced 76-base unique reads;
	- 9 reads over GTF introns: GT/AG class 4 (T01a, T05 GT/AG; T02, T06
	  CT/AC), GC/AG class 2 (T03 GC/AG, T08 CT/GC), AT/AC class 3 (T04
	  AT/AC, 2 x T11 GT/AT); one 200-base read over both T01a introns; 2
	  reads over the unannotated T10skip intron (GT/AG) and 2 over the
	  unannotated non-canonical T12nc intron;
	- 5 reads with 1 mismatch and 3 with 2 (11 mismatches); 3 reads with a
	  2-base and 1 with a 3-base deletion (9 bases, 4 deletions); 2 reads
	  with a 3-base and 1 with a 1-base insertion (7 bases, 3 insertions);
	- 6 dup3 (3 loci) and 2 dup2 (2 loci) reads: multi-mapping;
	- 4 rep12 reads: 12 loci > 10, too many loci;
	- 3 reads of 55 genomic + 45 random bases: best alignment < 0.66 of the
	  read, unmapped too short; 2 all-N reads: no seeds, unmapped other.
	By hand: input 86, bases 75 x 100 + 10 x 76 + 200 = 8460, average
	8460 / 86 = 98.4 -> 98; unique 69 (80.23%); mapped bases of unique
	reads 6753 (insertion reads map 97, 97 and 99 bases), average 97.87;
	splices 15 = GT/AG 4 + 2 + 2 = 8, GC/AG 2, AT/AC 3, non-canonical 2,
	annotated 9 + 2 = 11; mismatch rate 11 / 6753 = 0.16%; deletion rate
	9 / 6753 = 0.13%, average 2.25; insertion rate 7 / 6753 = 0.10%,
	average 2.33; multi 8 (9.30%); too many loci 4 (4.65%); too many
	mismatches 0; too short 3 (3.49%); other 2 (2.33%); chimeric 0."""
	g = sim.g
	reads = _unique_genomic(sim, 30) + _unique_genomic(sim, 10, length=76)
	spliced = [("chr1", (1301, 2300), 50, 50, "+"), ("chr1", (20601, 21000), 40, 60, "-"),
		("chr1", (6401, 7000), 30, 70, "+"), ("chr1", (21801, 22300), 60, 40, "-"),
		("chr1", (11401, 11800), 50, 50, "+"), ("chr2", (2301, 3000), 45, 55, "-"),
		("chr1", (15301, 15900), 50, 50, "+"), ("chr3", (12501, 13000), 35, 65, "+"), ("chr3", (12501, 13000), 65, 35, "-"),
		("chr2", (30401, 31800), 50, 50, "+"), ("chr2", (30401, 31800), 30, 70, "-"),
		("chr3", (15201, 15800), 50, 50, "+"), ("chr3", (15201, 15800), 40, 60, "-")]
	reads += [hc.junction_read(sim, c, i, l, r, s) for c, i, l, r, s in spliced]
	assert [hc.motif_of(g, c, i[0], i[1]) for c, i, _, _, _ in spliced] == [1, 1, 2, 2, 3, 4, 5, 6, 6, 1, 1, 0, 0]
	reads.append(sim.transcript_read("T01a", 280, length=200))
	reads += _unique_genomic(sim, 5, mismatches=(40,)) + _unique_genomic(sim, 3, mismatches=(30, 70))
	reads += [_deletion_read(sim, 2) for _ in range(3)] + [_deletion_read(sim, 3)]
	for n in (3, 3, 1):
		c, p = sim.unique_locus(length=110, avoid_genes=True)
		reads.append(sim.genomic_read(c, p, "-", insertion=(50, n)))
	unique = list(reads)
	reads += [sim.segment_read("dup3", o) for o in (0, 90, 180, 270, 360, 450)]
	reads += [sim.segment_read("dup2", o) for o in (200, 900)]
	reads += [sim.segment_read("rep12", o) for o in (10, 70, 130, 190)]
	for _ in range(3):
		c, p = sim.unique_locus(avoid_genes=True)
		base = sim.genomic_read(c, p, "+", length=55)
		reads.append(sh.SimRead(sim.name(), base.seq + sim.random_seq(45), sim.qual(100), kind="random"))
	reads += [sh.SimRead(sim.name(), "N" * 100, sim.qual(100), kind="random") for _ in range(2)]
	assert len(unique) == 69 and len(reads) == 86
	total_bases = sum(len(r.seq) for r in reads)
	mapped = sum(n for r in unique for n, op in sh.cigar_ops(r.cigar) if op == "M")
	assert (total_bases, mapped) == (8460, 6753)
	run = run_star(genome_sjdb, ReadSet(reads))
	recs = run.alignments().by_name()
	for r in unique:
		assert [sh.record_truth(x) for x in recs[r.name]] == [r.truth()], r.name
	log = run.log_final()
	pct = lambda n: hc.fmt2(100.0 * n / 86)
	expected = {
		"Number of input reads": 86,
		"Average input read length": 98,
		"Uniquely mapped reads number": 69,
		"Uniquely mapped reads %": pct(69),
		"Average mapped length": hc.fmt2(6753 / 69.0),
		"Number of splices: Total": 15,
		"Number of splices: Annotated (sjdb)": 11,
		"Number of splices: GT/AG": 8,
		"Number of splices: GC/AG": 2,
		"Number of splices: AT/AC": 3,
		"Number of splices: Non-canonical": 2,
		"Mismatch rate per base, %": hc.fmt2(100.0 * 11 / 6753),
		"Deletion rate per base": hc.fmt2(100.0 * 9 / 6753),
		"Deletion average length": 2.25,
		"Insertion rate per base": hc.fmt2(100.0 * 7 / 6753),
		"Insertion average length": hc.fmt2(7 / 3.0),
		"Number of reads mapped to multiple loci": 8,
		"% of reads mapped to multiple loci": pct(8),
		"Number of reads mapped to too many loci": 4,
		"% of reads mapped to too many loci": pct(4),
		"Number of reads unmapped: too many mismatches": 0,
		"% of reads unmapped: too many mismatches": 0.0,
		"Number of reads unmapped: too short": 3,
		"% of reads unmapped: too short": pct(3),
		"Number of reads unmapped: other": 2,
		"% of reads unmapped: other": pct(2),
		"Number of chimeric reads": 0,
		"% of chimeric reads": 0.0,
	}
	assert {k: log[k] for k in expected} == expected
	assert (expected["Uniquely mapped reads %"], expected["Mismatch rate per base, %"],
		expected["Insertion average length"]) == (80.23, 0.16, 2.33)


def test_log_final_too_many_mismatches(genome_plain, sim):
	"""--alignEndsType EndToEnd forbids soft-clipping (parametersDefault),
	so a read whose first and last bases mismatch has 2 mismatches; with
	--outFilterMismatchNmax 1 its alignment is rejected and the read is
	counted as "unmapped: too many mismatches" (manual 5.2.2, uT attribute: type
	2). Of 8 reads: 3 such reads, 4 exact and 1 with one mismatch in the
	middle (5 unique, 62.50%); too many mismatches 3 (37.50%); mismatch
	rate 1 / 500 bases = 0.20%."""
	bad = _unique_genomic(sim, 3, mismatches=(0, 99))
	good = _unique_genomic(sim, 4) + _unique_genomic(sim, 1, mismatches=(50,))
	run = run_star(genome_plain, ReadSet(bad + good),
		"--alignEndsType EndToEnd --outFilterMismatchNmax 1 --outSAMunmapped Within")
	recs = run.alignments().by_name()
	for r in good:
		assert [sh.record_truth(x) for x in recs[r.name]] == [r.truth()]
	for r in bad:
		assert [(x.is_unmapped, x.tag("uT")) for x in recs[r.name]] == [(True, "2")]
	log = run.log_final()
	assert log["Number of input reads"] == 8
	assert log["Uniquely mapped reads number"] == 5
	assert log["Uniquely mapped reads %"] == 62.5
	assert log["Number of reads unmapped: too many mismatches"] == 3
	assert log["% of reads unmapped: too many mismatches"] == 37.5
	assert log["Number of reads unmapped: too short"] == 0
	assert log["Number of reads unmapped: other"] == 0
	assert log["Mismatch rate per base, %"] == 0.2


def test_log_final_chimeric(genome_plain, sim):
	"""With --chimSegmentMin 20 chimeric detection runs in addition to
	normal mapping (manual 6): 3 reads made of two 50-base pieces from
	different contigs are counted as chimeric (3 of 8 reads, 37.50%), and
	since their best linear alignment covers 50 of 100 bases (< 0.66) they
	are also unmapped "too short"; the 5 ordinary reads are unique."""
	chim = [sim.clean_chimeric_read() for _ in range(3)]
	plain = _unique_genomic(sim, 5)
	run = run_star(genome_plain, ReadSet(chim + plain), "--chimSegmentMin 20")
	assert sorted(j["read_name"] for j in run.chimeric()) == sorted(r.name for r in chim)
	log = run.log_final()
	assert log["Number of input reads"] == 8
	assert log["Uniquely mapped reads number"] == 5
	assert log["Number of chimeric reads"] == 3
	assert log["% of chimeric reads"] == 37.5
	assert log["Number of reads unmapped: too short"] == 3
