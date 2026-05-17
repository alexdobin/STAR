# STAR Genome Insert

`--runMode genomeInsert` adds named FASTA records to an existing STAR genome index. It is intended for small sequence additions such as spike-ins, controls, plasmid/transgene records, or other named sequences that should become normal reference contigs for subsequent alignment runs.

The default output mode writes a complete updated index. Overlay modes instead write a small directory that points to the prebuilt base index and inserted FASTA/GTF inputs, so the base index does not need to be copied or rebuilt.

## Usage

```bash
STAR \
  --runMode genomeInsert \
  --runThreadN 16 \
  --genomeDir /path/to/prebuilt/index \
  --genomeFastaFiles customer_transgenes.fa \
  --genomeInsertOutDir /path/to/updated/index
```

`--genomeDir` is read as the base index and is not modified. With the default `--genomeInsertOutMode Full`, `--genomeInsertOutDir` receives a complete updated index containing `Genome`, `SA`, `SAindex`, chromosome metadata, `genomeParameters.txt`, and supported annotation sidecar files.

Multiple inserted FASTA files can be supplied with `--genomeFastaFiles`:

```bash
STAR \
  --runMode genomeInsert \
  --genomeDir /path/to/prebuilt/index \
  --genomeFastaFiles gfp.fa gst.fa \
  --genomeInsertOutDir /path/to/updated/index
```

Annotations for inserted sequences can be supplied as an insert-only GTF:

```bash
STAR \
  --runMode genomeInsert \
  --genomeDir /path/to/prebuilt/index \
  --genomeFastaFiles customer_transgenes.fa \
  --sjdbGTFfile customer_transgenes.gtf \
  --genomeInsertOutDir /path/to/updated/index
```

The inserted GTF must contain annotations only for sequences supplied in the same `--genomeFastaFiles` command. Existing annotations from the base index are preserved automatically. If the inserted GTF contains exons on a base genome chromosome, STAR exits with an input error to avoid duplicating base gene models.

## Output Modes

`--genomeInsertOutMode Full` writes a complete updated index. Use it when the updated reference should be self-contained and reusable without access to the original base index or inserted FASTA/GTF files.

`--genomeInsertOutMode Overlay` writes only `genomeInsertOverlay.tsv`. The overlay directory can be passed as `--genomeDir` for later alignment; STAR loads the base index listed in the manifest and inserts the extra FASTA/GTF records at alignment startup.

```bash
STAR \
  --runMode genomeInsert \
  --genomeDir /path/to/prebuilt/index \
  --genomeFastaFiles customer_transgenes.fa \
  --sjdbGTFfile customer_transgenes.gtf \
  --genomeInsertOutMode Overlay \
  --genomeInsertOutDir /path/to/customer/overlay

STAR \
  --genomeDir /path/to/customer/overlay \
  --readFilesIn reads.fq
```

`--genomeInsertOutMode Delta` writes an overlay plus `genomeInsertDelta.bin`, a cached suffix insertion plan for the inserted FASTA. Later alignments use the cached plan to avoid repeating suffix search and sorting before replaying the in-memory insertion.

```bash
STAR \
  --runMode genomeInsert \
  --genomeDir /path/to/prebuilt/index \
  --genomeFastaFiles customer_transgenes.fa \
  --sjdbGTFfile customer_transgenes.gtf \
  --genomeInsertOutMode Delta \
  --genomeInsertOutDir /path/to/customer/delta

STAR \
  --genomeDir /path/to/customer/delta \
  --readFilesIn reads.fq
```

The delta file is validated against the loaded base index dimensions and inserted sequence content before it is used. Delta mode still constructs the expanded genome index in memory at alignment startup; it does not write or require a persistent `Genome`, `SA`, or `SAindex` in the overlay directory.

## Validation Expectations

An updated index should be validated against a full rebuild from the base FASTA/GTF files plus the inserted FASTA/GTF files. The following outputs are expected to match byte-for-byte when the base build inputs and parameters are identical:

- `Genome`
- `SAindex`
- `chrName.txt`
- `chrStart.txt`
- `chrLength.txt`
- `chrNameLength.txt`
- annotation sidecar files such as `sjdbInfo.txt`, `sjdbList.out.tab`, `exonInfo.tab`, `geneInfo.tab`, and `transcriptInfo.tab`

`SA` may not be byte-identical to a full rebuild because suffixes with equal ordering keys can be stored in a different but equivalent order. Validate `SA` behavior by aligning reads that target the inserted sequences and comparing alignment records against a full rebuild.

## Limitations

- `genomeInsert` adds new named FASTA sequences. It does not edit existing chromosome sequence.
- `genomeInsert` does not update in place; it writes a complete new index, overlay, or delta overlay to `--genomeInsertOutDir`.
- `--sjdbGTFfile` is supported only for annotations on inserted sequences. It is not a way to edit or replace existing base genome annotations.
- Overlay and delta directories depend on the base genome index and inserted FASTA/GTF paths recorded in `genomeInsertOverlay.tsv`.
- `--sjdbFileChrStartEnd` and `--twopassMode` are not supported with `--runMode genomeInsert`.

## Test and Benchmark Helpers

Fast regression test:

```bash
make -C source STAR
extras/tests/scripts/testGenomeInsert.sh
```

Benchmark harness:

```bash
BASE_GENOME_DIR=/path/to/prebuilt/index \
BASE_FASTA_FILES="/path/reference.fa /path/ercc.fa" \
INSERT_FASTA_FILES="/path/gfp.fa /path/gst.fa" \
INSERT_GTF_FILE=/path/inserted_sequences.gtf \
SJDB_GTF_FILE=/path/base_plus_inserted.gtf \
SJDB_OVERHANG=93 \
THREADS=96 \
extras/tests/scripts/benchmarkGenomeInsert.sh
```

The benchmark harness records command lines, `time -v` metrics, iostat snapshots when available, and byte-level validation for files expected to match a full rebuild.
