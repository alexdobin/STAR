"""Fixtures for the STAR unit tests; tests/README.md documents the API."""

import os
import zlib

import pytest

import star_helpers as sh


def pytest_report_header(config):
	star = sh.find_star()
	return ["STAR binary: %s" % (star or "not found"), "STAR_TEST_THREADS: %s" % sh.default_threads()]


@pytest.fixture(scope="session", autouse=True)
def star_bin(tmp_path_factory):
	"""Path of the STAR binary. Skips the whole session when it is missing."""
	star = sh.find_star()
	if star is None:
		pytest.skip(sh.BUILD_HINT)
	sh.set_run_root(tmp_path_factory.mktemp("runs"))
	return star


@pytest.fixture(scope="session")
def star_threads():
	"""The default --runThreadN for this session ($STAR_TEST_THREADS or 1)."""
	return sh.default_threads()


@pytest.fixture(scope="session")
def synth_genome(tmp_path_factory):
	"""The default SyntheticGenome, with genome.fa and annotation.gtf written."""
	g = sh.make_synthetic_genome()
	g.write(tmp_path_factory.mktemp("synth"))
	return g


@pytest.fixture(scope="session")
def genome_factory(synth_genome, tmp_path_factory):
	"""make(gtf=True, sjdb_overhang=99, extra_args=(), synth=None) -> StarGenome,
	built once per distinct argument set for the session."""
	cache = {}

	def make(gtf=True, sjdb_overhang=sh.DEFAULT_SJDB_OVERHANG, extra_args=(), synth=None):
		synth = synth or synth_genome
		key = (id(synth), str(gtf), sjdb_overhang if gtf else None, tuple(extra_args))
		if key not in cache:
			d = tmp_path_factory.mktemp("genome")
			cache[key] = sh.build_genome(synth, d, gtf=gtf, sjdb_overhang=sjdb_overhang, extra_args=extra_args)
		return cache[key]

	return make


@pytest.fixture(scope="session")
def genome_plain(genome_factory):
	"""STAR index of the synthetic genome without annotations."""
	return genome_factory(gtf=False)


@pytest.fixture(scope="session")
def genome_sjdb(genome_factory):
	"""STAR index with the synthetic GTF inserted (--sjdbOverhang 99)."""
	return genome_factory(gtf=True)


def stable_seed(*parts):
	return zlib.crc32("::".join(str(p) for p in parts).encode())


@pytest.fixture
def sim(synth_genome, request):
	"""ReadSimulator seeded from the test's module and name, so each test
	draws the same reads however the suite is invoked or ordered."""
	module = os.path.splitext(os.path.basename(request.node.module.__file__))[0]
	return sh.ReadSimulator(synth_genome, seed=stable_seed(module, request.node.name))


@pytest.fixture(scope="session")
def tenx_whitelist(tmp_path_factory):
	"""(barcodes, path): 24 cell barcodes of 16 bp written one per line."""
	barcodes = sh.make_whitelist(24, seed=11)
	path = tmp_path_factory.mktemp("whitelist") / "whitelist.txt"
	path.write_text("\n".join(barcodes) + "\n")
	return barcodes, path
